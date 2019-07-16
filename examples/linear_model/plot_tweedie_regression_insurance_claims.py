"""
======================================
Tweedie regression on insurance claims
======================================

This example illustrates the use of Poisson, Gamma and Tweedie regression
on the French Motor Third-Party Liability Claims dataset, and is inspired
by an R tutorial [1].

Insurance claims data consist of the number of claims and the total claim
amount. Often, the final goal is to predict the expected value, i.e. the mean,
of the total claim amount. There are several possibilities to do that, two of
which are:

1. Model the number of claims with a Poisson distribution, the average
   claim amount as a Gamma distribution and multiply the predictions of both in
   order to get the total claim amount.
2. Model total claim amount directly, typically with a Tweedie distribution of
   Tweedie power :math:`p \\in (1, 2)`.

In this example we will illustrate both approaches. We start by defining a few
helper functions for loading the data and visualizing results.


.. [1]  A. Noll, R. Salzmann and M.V. Wuthrich, Case Study: French Motor
    Third-Party Liability Claims (November 8, 2018).
    `doi:10.2139/ssrn.3164764 <http://dx.doi.org/10.2139/ssrn.3164764>`_

"""
print(__doc__)

# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#          Roman Yurchak <rth.yurchak@gmail.com>
# License: BSD 3 clause
from functools import partial

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.compose import ColumnTransformer
from sklearn.linear_model import GeneralizedLinearRegressor
from sklearn.linear_model._glm.distribution import TweedieDistribution
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder
from sklearn.preprocessing import StandardScaler, KBinsDiscretizer

from sklearn.metrics import mean_absolute_error, mean_squared_error


def load_mtpl2(n_samples=100000):
    """Fetch the French Motor Third-Party Liability Claims dataset.

    Parameters
    ----------
    n_samples: int, default=100000
      number of samples to select (for faster run time).
    """

    # Note: this should use the OpenML DataFrame fetcher in the future
    df_freq = pd.read_csv(
        "https://www.openml.org/data/get_csv/20649148/freMTPL2freq.csv",
        dtype={"IDpol": np.int},
        index_col=0,
    )

    df_sev = pd.read_csv(
        "https://www.openml.org/data/get_csv/20649149/freMTPL2sev.arff",
        index_col=0,
    )

    # sum ClaimAmount over identical IDs
    df_sev = df_sev.groupby(level=0).sum()

    df = df_freq.join(df_sev, how="left")
    df["ClaimAmount"].fillna(0, inplace=True)

    # unquote string fields
    for column_name in df.columns[df.dtypes.values == np.object]:
        df[column_name] = df[column_name].str.strip("'")
    return df.iloc[:n_samples]


def plot_obs_pred(df, feature, weight, observed, predicted, y_label=None,
                  title=None, ax=None):
    """Plot observed and predicted - aggregated per feature level.

    Parameters
    ----------
    df : DataFrame with at least three columns named feature, weight and
         observed
    feature: str
        a column name of df for the feature to be plotted
    weight : str
        column name of df with the values of weights or exposure
    observed : str
        a column name of df with the observed target
    predicted : frame
        a dataframe, with the same index as df, with the predicted target
    """
    # aggregate observed and predicted variables by feature level
    df_ = df.loc[:, [feature, weight]].copy()
    df_["observed"] = df[observed] * df[weight]
    df_["predicted"] = predicted * df[weight]
    df_ = (
        df_.groupby([feature])[weight, "observed", "predicted"]
        .sum()
        .assign(observed=lambda x: x["observed"] / x[weight])
        .assign(predicted=lambda x: x["predicted"] / x[weight])
    )

    ax = df_.loc[:, ["observed", "predicted"]].plot(style=".", ax=ax)
    y_max = df_.loc[:, ["observed", "predicted"]].values.max() * 0.8
    ax.fill_between(
        df_.index,
        0,
        y_max * df_[weight] / df_[weight].values.max(),
        color="g",
        alpha=0.1,
    )
    ax.set(
        ylabel=y_label if y_label is not None else None,
        title=title if title is not None else "Train: Observed vs Predicted",
    )


##############################################################################
#
# 1. Loading datasets and pre-processing
# --------------------------------------
#
# We construct the freMTPL2 dataset by joining the freMTPL2freq table,
# containing the number of claims (``ClaimNb``), with the freMTPL2sev table,
# containing the claim amount (``ClaimAmount``) for the same policy ids
# (``IDpol``).

df = load_mtpl2(n_samples=100000)

# Note: filter out claims with zero amount, as the severity model
# requires a strictly positive target values.
df.loc[(df.ClaimAmount == 0) & (df.ClaimNb >= 1), "ClaimNb"] = 0

# Correct for unreasonable observations (that might be data error)
# and a few exceptionally large claim amounts
df["ClaimNb"] = df["ClaimNb"].clip(upper=4)
df["Exposure"] = df["Exposure"].clip(upper=1)
df["ClaimAmount"] = df["ClaimAmount"].clip(upper=200000)

column_trans = ColumnTransformer(
    [
        ("Veh_Driv_Age", KBinsDiscretizer(n_bins=10), ["VehAge", "DrivAge"]),
        (
            "Veh_Brand_Gas_Region",
            OneHotEncoder(),
            ["VehBrand", "VehPower", "VehGas", "Region", "Area"],
        ),
        ("BonusMalus", "passthrough", ["BonusMalus"]),
        (
            "Density_log",
            make_pipeline(
                FunctionTransformer(np.log, validate=False), StandardScaler()
            ),
            ["Density"],
        ),
    ],
    remainder="drop",
)
X = column_trans.fit_transform(df)


df["Frequency"] = df.ClaimNb / df.Exposure
df["AvgClaimAmount"] = df.ClaimAmount / np.fmax(df.ClaimNb, 1)

print(df[df.ClaimAmount > 0].head())

##############################################################################
#
# 2. Frequency model -- Poisson distribution
# -------------------------------------------
#
# The number of claims (``ClaimNb``) is a positive integer that can be modeled
# as a Poisson distribution. It is then assumed to be the number of discrete
# events occuring with a constant rate in a given time interval (``Exposure``).
# Here we model the frequency ``y = ClaimNb / Exposure``,
# which is still a (scaled) Poisson distribution.
#
# A very important property of the Poisson distribution is its mean-variance
# relation: The variance is proportional to the mean.

df_train, df_test, X_train, X_test = train_test_split(df, X, random_state=2)

# Some of the features are colinear, we use a weak penalization to avoid
# numerical issues.
glm_freq = GeneralizedLinearRegressor(family="poisson", alpha=1e-2)
glm_freq.fit(X_train, df_train.Frequency, sample_weight=df_train.Exposure)


def mean_deviance(estimator, y, y_pred, weights):
    if hasattr(estimator, "_family_instance"):
        return estimator._family_instance.deviance(y, y_pred, weights) / len(y)
    else:
        return np.nan


def score_estimator(
    estimator, X_train, X_test, df_train, df_test, target, weights
):
    res = []

    for subset_label, X, df in [
        ("train", X_train, df_train),
        ("test", X_test, df_test),
    ]:
        y, _weights = df[target], df[weights]

        for score_label, metric in [
            ("D² explained", None),
            ("mean deviance", partial(mean_deviance, estimator)),
            ("mean abs. error", mean_absolute_error),
            ("mean squared error", mean_squared_error),
        ]:
            if estimator.__class__.__name__ == "ClaimProdEstimator":
                # ClaimProdEstimator is the product of frequency and severity
                # models, denormalized by the exposure values.
                # It does not fully follow the scikit-learn API and we
                # must handle it separately.
                y_pred = estimator.predict(X, exposure=df.Exposure.values)
            else:
                y_pred = estimator.predict(X)
            if metric is None:
                if not hasattr(estimator, "score"):
                    continue
                score = estimator.score(X, y, _weights)
            else:
                score = metric(y, y_pred, _weights)

            res.append(
                {"subset": subset_label, "metric": score_label, "score": score}
            )

    res = (
        pd.DataFrame(res)
        .set_index(["metric", "subset"])
        .score.unstack(-1)
        .round(3)
    )
    return res


scores = score_estimator(
    glm_freq,
    X_train,
    X_test,
    df_train,
    df_test,
    target="Frequency",
    weights="Exposure",
)
print(scores)

##############################################################################
#
# We can visually compare observed and predicted values, aggregated by
# the drivers age (``DrivAge``), vehicle age (``VehAge``) and the insurance
# bonus/malus (``BonusMalus``).

fig, ax = plt.subplots(2, 2, figsize=(16, 8))
fig.subplots_adjust(hspace=0.3, wspace=0.2)

plot_obs_pred(
    df=df_train,
    feature="DrivAge",
    weight="Exposure",
    observed="Frequency",
    predicted=glm_freq.predict(X_train),
    y_label="Claim Frequency",
    title="train data",
    ax=ax[0, 0],
)

plot_obs_pred(
    df=df_test,
    feature="DrivAge",
    weight="Exposure",
    observed="Frequency",
    predicted=glm_freq.predict(X_test),
    y_label="Claim Frequency",
    title="test data",
    ax=ax[0, 1],
)

plot_obs_pred(
    df=df_test,
    feature="VehAge",
    weight="Exposure",
    observed="Frequency",
    predicted=glm_freq.predict(X_test),
    y_label="Claim Frequency",
    title="test data",
    ax=ax[1, 0],
)

plot_obs_pred(
    df=df_test,
    feature="BonusMalus",
    weight="Exposure",
    observed="Frequency",
    predicted=glm_freq.predict(X_test),
    y_label="Claim Frequency",
    title="test data",
    ax=ax[1, 1],
)


##############################################################################
#
# 3. Severity model -  Gamma Distribution
# ---------------------------------------
# The mean claim amount or severity (`AvgClaimAmount`) can be empirically
# shown to follow approximately a Gamma distribution. We fit a GLM model for
# the severity with the same features as the frequency model.
#
# Note:
#
# - We filter out ``ClaimAmount == 0`` as the Gamma distribution has support
#   on :math:`(0, \infty)`, not :math:`[0, \infty)`.
# - We use ``ClaimNb`` as sample weights.

mask_train = df_train["ClaimAmount"] > 0
mask_test = df_test["ClaimAmount"] > 0

glm_sev = GeneralizedLinearRegressor(family="gamma")

glm_sev.fit(
    X_train[mask_train.values],
    df_train.loc[mask_train, "AvgClaimAmount"],
    sample_weight=df_train.loc[mask_train, "ClaimNb"],
)


scores = score_estimator(
    glm_sev,
    X_train[mask_train.values],
    X_test[mask_test.values],
    df_train[mask_train],
    df_test[mask_test],
    target="AvgClaimAmount",
    weights="ClaimNb",
)
print(scores)

##############################################################################
#
# Note that the resulting model is the average claim amount per claim. As such,
# it is conditional on having at least one claim, and cannot be used to predict
# the average claim amount per policy in general.

print(
    "Mean AvgClaim Amount per policy:              %.2f "
    % df_train.AvgClaimAmount.mean()
)
print(
    "Mean AvgClaim Amount | NbClaim > 0:           %.2f"
    % df_train.AvgClaimAmount[df_train.AvgClaimAmount > 0].mean()
)
print(
    "Predicted Mean AvgClaim Amount | NbClaim > 0: %.2f"
    % glm_sev.predict(X_train).mean()
)


##############################################################################
#
# We can visually compare observed and predicted values, aggregated for
# the drivers age (``DrivAge``).

fig, ax = plt.subplots(1, 2, figsize=(16, 4))

# plot DivAge
plot_obs_pred(
    df=df_train.loc[mask_train],
    feature="DrivAge",
    weight="Exposure",
    observed="AvgClaimAmount",
    predicted=glm_sev.predict(X_train[mask_train.values]),
    y_label="Average Claim Severity",
    title="train data",
    ax=ax[0],
)

plot_obs_pred(
    df=df_test.loc[mask_test],
    feature="DrivAge",
    weight="Exposure",
    observed="AvgClaimAmount",
    predicted=glm_sev.predict(X_test[mask_test.values]),
    y_label="Average Claim Severity",
    title="test data",
    ax=ax[1],
)


##############################################################################
#
# 4. Total Claims Amount -- Compound Poisson distribution
# -------------------------------------------------------
#
# As mentionned in the introduction, the total claim amount can be modeled
# either as the product of the frequency model by the severity model,


class ClaimProdEstimator:
    """Total claim amount estimator.

    Computed as the product of the frequency model by the serverity model,
    denormalized by exposure. Use Tweedie deviance with `p=1.5`.
    """

    def __init__(self, est_freq, est_sev):
        self.est_freq = est_freq
        self.est_sev = est_sev
        self._family_instance = TweedieDistribution(power=1.5)

    def predict(self, X, exposure):
        """Predict the total claim amount.

        The predict method is not compatible with the scikit-learn API.
        """
        return exposure * self.est_freq.predict(X) * self.est_sev.predict(X)

    def score(self, X, y, sample_weight=None):
        """Compute D², the percentage of deviance explained."""
        mu = self.predict(X, exposure=sample_weight)
        dev = self._family_instance.deviance(y, mu, weights=sample_weight)
        y_mean = np.average(y, weights=sample_weight)
        dev_null = self._family_instance.deviance(y, y_mean,
                                                  weights=sample_weight)
        return 1. - dev / dev_null


est_prod = ClaimProdEstimator(glm_freq, glm_sev)

scores = score_estimator(
    est_prod,
    X_train,
    X_test,
    df_train,
    df_test,
    target="ClaimAmount",
    weights="Exposure",
)
print(scores)


##############################################################################
#
# or as a unique Compound Poisson model, also corresponding to a Tweedie model
# with a power :math:`p \in (1, 2)`. We determine the optimal hyperparameter
# ``p`` with a grid search,

from sklearn.model_selection import GridSearchCV

# this takes a while
params = {
    "family": [
        TweedieDistribution(power=power) for power in np.linspace(1, 2, 8)
    ]
}

glm_total = GridSearchCV(
    GeneralizedLinearRegressor(), cv=3, param_grid=params, n_jobs=-1
)
glm_total.fit(
    X_train, df_train["ClaimAmount"], sample_weight=df_train["Exposure"]
)


print(
    "Best hyperparameters: power=%.2f\n"
    % glm_total.best_estimator_.family.power
)

scores = score_estimator(
    glm_total.best_estimator_,
    X_train,
    X_test,
    df_train,
    df_test,
    target="ClaimAmount",
    weights="Exposure",
)
print(scores)

##############################################################################
#
# In this example, the mean absolute error is lower for the Compound Poisson
# model than when using separate models for frequency and severity.
#
# We can additionally validate these models by comparing observed and predicted
# total claim amount over the test and train subsets. We see that in our case
# the frequency-severity model underestimates the total claim amount, whereas
# the Tweedie model overestimates.

res = []
for subset_label, X, df in [
    ("train", X_train, df_train),
    ("test", X_test, df_test),
]:
    res.append(
        {
            "subset": subset_label,
            "observed": df.ClaimAmount.values.sum(),
            "predicted, frequency*severity model": np.sum(
                est_prod.predict(X, exposure=df.Exposure.values)
            ),
            "predicted, tweedie, p=%.2f"
            % glm_total.best_estimator_.family.power: np.sum(
                glm_total.best_estimator_.predict(X)
            ),
        }
    )

print(pd.DataFrame(res).set_index("subset").T)
