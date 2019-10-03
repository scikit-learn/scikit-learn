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
   claim amount per claim, also known as severity, as a Gamma distribution and
   multiply the predictions of both in order to get the total claim amount.
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

from sklearn.datasets import fetch_openml
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import PoissonRegressor, GammaRegressor
from sklearn.linear_model import TweedieRegressor
from sklearn.metrics import mean_tweedie_deviance
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
      number of samples to select (for faster run time). Full dataset has
      678013 samples.
    """

    # freMTPL2freq dataset from https://www.openml.org/d/41214
    df_freq = fetch_openml(data_id=41214, as_frame=True)['data']
    df_freq['IDpol'] = df_freq['IDpol'].astype(np.int)
    df_freq.set_index('IDpol', inplace=True)

    # freMTPL2sev dataset from https://www.openml.org/d/41215
    df_sev = fetch_openml(data_id=41215, as_frame=True)['data']

    # sum ClaimAmount over identical IDs
    df_sev = df_sev.groupby('IDpol').sum()

    df = df_freq.join(df_sev, how="left")
    df["ClaimAmount"].fillna(0, inplace=True)

    # unquote string fields
    for column_name in df.columns[df.dtypes.values == np.object]:
        df[column_name] = df[column_name].str.strip("'")
    return df.iloc[:n_samples]


def plot_obs_pred(df, feature, weight, observed, predicted, y_label=None,
                  title=None, ax=None, fill_legend=False):
    """Plot observed and predicted - aggregated per feature level.

    Parameters
    ----------
    df : DataFrame
        input data
    feature: str
        a column name of df for the feature to be plotted
    weight : str
        column name of df with the values of weights or exposure
    observed : str
        a column name of df with the observed target
    predicted : frame
        a dataframe, with the same index as df, with the predicted target
    fill_legend : bool, default=False
        whether to show fill_between legend
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
    p2 = ax.fill_between(
        df_.index,
        0,
        y_max * df_[weight] / df_[weight].values.max(),
        color="g",
        alpha=0.1,
    )
    if fill_legend:
        ax.legend([p2], ["{} distribution".format(feature)])
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

df = load_mtpl2(n_samples=60000)

# Note: filter out claims with zero amount, as the severity model
# requires strictly positive target values.
df.loc[(df["ClaimAmount"] == 0) & (df["ClaimNb"] >= 1), "ClaimNb"] = 0

# Correct for unreasonable observations (that might be data error)
# and a few exceptionally large claim amounts
df["ClaimNb"] = df["ClaimNb"].clip(upper=4)
df["Exposure"] = df["Exposure"].clip(upper=1)
df["ClaimAmount"] = df["ClaimAmount"].clip(upper=200000)

log_scale_transformer = make_pipeline(
    FunctionTransformer(np.log, validate=False),
    StandardScaler()
)

column_trans = ColumnTransformer(
    [
        ("binned_numeric", KBinsDiscretizer(n_bins=10), ["VehAge", "DrivAge"]),
        ("onehot_categorical", OneHotEncoder(),
         ["VehBrand", "VehPower", "VehGas", "Region", "Area"]),
        ("passthrough_numeric", "passthrough", ["BonusMalus"]),
        ("log_scaled_numeric", log_scale_transformer, ["Density"]),
    ],
    remainder="drop",
)
X = column_trans.fit_transform(df)


df["Frequency"] = df["ClaimNb"] / df["Exposure"]
df["AvgClaimAmount"] = df["ClaimAmount"] / np.fmax(df["ClaimNb"], 1)

print(df[df.ClaimAmount > 0].head())

##############################################################################
#
# 2. Frequency model -- Poisson distribution
# -------------------------------------------
#
# The number of claims (``ClaimNb``) is a positive integer that can be modeled
# as a Poisson distribution. It is then assumed to be the number of discrete
# events occuring with a constant rate in a given time interval
# (``Exposure``, in units of years). Here we model the frequency
# ``y = ClaimNb / Exposure``, which is still a (scaled) Poisson distribution,
# and use ``Exposure`` as `sample_weight`.

df_train, df_test, X_train, X_test = train_test_split(df, X, random_state=0)

# Some of the features are colinear, we use a weak penalization to avoid
# numerical issues.
glm_freq = PoissonRegressor(alpha=1e-2)
glm_freq.fit(X_train, df_train.Frequency, sample_weight=df_train.Exposure)


def score_estimator(
    estimator, X_train, X_test, df_train, df_test, target, weights
):
    """Evaluate an estimator on train and test sets with different metrics"""
    res = []

    for subset_label, X, df in [
        ("train", X_train, df_train),
        ("test", X_test, df_test),
    ]:
        y, _weights = df[target], df[weights]

        for score_label, metric in [
            ("D² explained", None),
            ("mean deviance", mean_tweedie_deviance),
            ("mean abs. error", mean_absolute_error),
            ("mean squared error", mean_squared_error),
        ]:
            if isinstance(estimator, tuple) and len(estimator) == 2:
                # Score the model consisting of the product of frequency and
                # severity models, denormalized by the exposure values.
                est_freq, est_sev = estimator
                y_pred = (df.Exposure.values * est_freq.predict(X) *
                          est_sev.predict(X))
                power = 1.5
            else:
                y_pred = estimator.predict(X)
                power = getattr(getattr(estimator, "_family_instance"),
                                "power")

            if score_label == "mean deviance":
                metric = partial(mean_tweedie_deviance, power=power)

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
        .round(2)
        .loc[:, ['train', 'test']]
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

fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(16, 8))
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
    fill_legend=True
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
    fill_legend=True
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
    fill_legend=True
)


##############################################################################
#
# According to the observed data, the frequency of accidents is higher for
# drivers younger than 30 years old, and it positively correlated with the
# `BonusMalus` variable. Our model is able to mostly correctly model
# this behaviour.
#
# 3. Severity model -  Gamma distribution
# ---------------------------------------
# The mean claim amount or severity (`AvgClaimAmount`) can be empirically
# shown to follow approximately a Gamma distribution. We fit a GLM model for
# the severity with the same features as the frequency model.
#
# Note:
#
# - We filter out ``ClaimAmount == 0`` as the Gamma distribution has support
#   on :math:`(0, \infty)`, not :math:`[0, \infty)`.
# - We use ``ClaimNb`` as `sample_weight`.

mask_train = df_train["ClaimAmount"] > 0
mask_test = df_test["ClaimAmount"] > 0

glm_sev = GammaRegressor()

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
# Here, the scores for the test data call for caution as they are significantly
# worse than for the training data indicating an overfit.
# Note that the resulting model is the average claim amount per claim. As such,
# it is conditional on having at least one claim, and cannot be used to predict
# the average claim amount per policy in general.

print("Mean AvgClaim Amount per policy:              %.2f "
      % df_train["AvgClaimAmount"].mean())
print("Mean AvgClaim Amount | NbClaim > 0:           %.2f"
      % df_train["AvgClaimAmount"][df_train["AvgClaimAmount"] > 0].mean())
print("Predicted Mean AvgClaim Amount | NbClaim > 0: %.2f"
      % glm_sev.predict(X_train).mean())


##############################################################################
#
# We can visually compare observed and predicted values, aggregated for
# the drivers age (``DrivAge``).

fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(16, 4))

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
    fill_legend=True
)


##############################################################################
#
# Overall, the drivers age (``DrivAge``) has a weak impact on the claim
# severity, both in observed and predicted data.
#
# 4. Total claim amount -- Compound Poisson distribution
# -------------------------------------------------------
#
# As mentioned in the introduction, the total claim amount can be modeled
# either as the product of the frequency model by the severity model,
# denormalized by exposure. In the following code sample, the
# ``score_estimator`` is extended to score such a model. The mean deviance
# is computed assuming a Tweedie distribution with ``power=1.5`` to be
# comparable with the model from the following section,


scores = score_estimator(
    (glm_freq, glm_sev),
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
# Indeed, an alternative approach for modeling the total loss is with a unique
# Compound Poisson model, also corresponding to a Tweedie model
# with a power :math:`p \in (1, 2)`. We determine the optimal hyperparameter
# ``p`` with a grid search,

from sklearn.model_selection import GridSearchCV

# exclude upper bound as power>=2 does not support y=0.
params = {"power": np.linspace(1 + 1e-4, 2 - 1e-4, 8)}


# this takes a while
glm_total = GridSearchCV(
    TweedieRegressor(tol=1e-3, max_iter=500), cv=3,
    param_grid=params, n_jobs=-1
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
            "observed": df["ClaimAmount"].values.sum(),
            "predicted, frequency*severity model": np.sum(
                df["Exposure"].values*glm_freq.predict(X)*glm_sev.predict(X)
            ),
            "predicted, tweedie, power=%.2f"
            % glm_total.best_estimator_.family.power: np.sum(
                glm_total.best_estimator_.predict(X)
            ),
        }
    )

print(pd.DataFrame(res).set_index("subset").T)
plt.plot()
