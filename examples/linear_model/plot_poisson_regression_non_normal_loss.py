"""
======================================
Poisson regression and non-normal loss
======================================

This example illustrates the use of log-linear Poisson regression
on the French Motor Third-Party Liability Claims dataset [1] and compares
it with models learned with least squared error. The goal is to predict the
number of insurance claims (or frequency) following car accidents for a
policyholder given historical data over a population of policyholders.

We start by defining a few helper functions for loading the data and
visualizing results.


.. [1]  A. Noll, R. Salzmann and M.V. Wuthrich, Case Study: French Motor
    Third-Party Liability Claims (November 8, 2018).
    `doi:10.2139/ssrn.3164764 <http://dx.doi.org/10.2139/ssrn.3164764>`_

"""
print(__doc__)

# Authors: Christian Lorentzen <lorentzen.ch@gmail.com>
#          Roman Yurchak <rth.yurchak@gmail.com>
# License: BSD 3 clause
import warnings

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from sklearn.datasets import fetch_openml
from sklearn.dummy import DummyRegressor
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge, PoissonRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import FunctionTransformer, OneHotEncoder
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import StandardScaler, KBinsDiscretizer
from sklearn.ensemble import RandomForestRegressor
from sklearn.utils import gen_even_slices

from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.metrics import mean_poisson_deviance


def load_mtpl2(n_samples=100000):
    """Fetcher for French Motor Third-Party Liability Claims dataset

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


##############################################################################
#
# 1. Loading datasets and pre-processing
# --------------------------------------
#
# We construct the freMTPL2 dataset by joining the freMTPL2freq table,
# containing the number of claims (``ClaimNb``) with the freMTPL2sev table
# containing the claim amount (``ClaimAmount``) for the same policy ids
# (``IDpol``).

df = load_mtpl2(n_samples=50000)

# Note: filter out claims with zero amount, as the severity model
# requires strictly positive target values.
df.loc[(df.ClaimAmount == 0) & (df.ClaimNb >= 1), "ClaimNb"] = 0

# correct for unreasonable observations (that might be data error)
df["ClaimNb"] = df["ClaimNb"].clip(upper=4)
df["Exposure"] = df["Exposure"].clip(upper=1)

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

##############################################################################
#
# The number of claims (``ClaimNb``) is a positive integer that can be modeled
# as a Poisson distribution. It is then assumed to be the number of discrete
# events occurring with a constant rate in a given time interval
# (``Exposure``, in units of years). Here we model the frequency
# ``y = ClaimNb / Exposure``, which is still a (scaled) Poisson distribution,
# and use ``Exposure`` as `sample_weight`.

df["Frequency"] = df.ClaimNb / df.Exposure

print(
   pd.cut(df.Frequency, [-1e-6, 1e-6, 1, 2, 3, 4, 5]).value_counts()
)

print("Average Frequency = {}"
      .format(np.average(df.Frequency, weights=df.Exposure)))

##############################################################################
#
# It worth noting that 96 % of policyholders have zero claims, and if we were
# to convert this problem into a binary classification task, it would be
# significantly imbalanced.
#
# To evaluate the pertinence of the used metrics, we will consider as a
# baseline an estimator that returns the mean of the training sample.

df_train, df_test = train_test_split(df, random_state=0)

dummy = make_pipeline(
    column_trans,
    DummyRegressor(strategy='mean')
)
dummy.fit(df_train, df_train["Frequency"],
          dummyregressor__sample_weight=df_train["Exposure"])


def score_estimator(estimator, df_test):
    """Score an estimatr on the test set"""

    y_pred = estimator.predict(df_test)

    print("MSE: %.3f" % mean_squared_error(
              df_test["Frequency"], y_pred,
              df_test["Exposure"]))
    print("MAE: %.3f" % mean_absolute_error(
              df_test["Frequency"], y_pred,
              df_test["Exposure"]))

    # ignore negative predictions, as they are invalid for
    # the Poisson deviance
    mask = y_pred > 0
    if (~mask).any():
        warnings.warn("estimator yields negative predictions for {} samples "
                      "out of {}. These will be ignored while computing the "
                      "poisson deviance".format((~mask).sum(), mask.shape[0]))

    print("mean Poisson deviance: %.3f" % mean_poisson_deviance(
            df_test["Frequency"][mask],
            y_pred[mask],
            df_test["Exposure"][mask]))


print("Constant mean frequency evaluation:")
score_estimator(dummy, df_test)

##############################################################################
#
# We start by modeling the target variable with the least squares linear
# regression model,

ridge = make_pipeline(
    column_trans,
    Ridge(alpha=1.0)
)
ridge.fit(df_train, df_train["Frequency"],
          ridge__sample_weight=df_train["Exposure"])

##############################################################################
#
# The Poisson deviance cannot be computed on negative values predicted by the
# model. For models that do return a few negative predictions
# (e.g. :class:`linear_model.Ridge`) we ignore the corresponding samples,
# meaning that the obtained Poisson deviance is approximate. An alternative
# apporach could be to use class:`compose.TransformedTargetRegressor`
# meta-estimator to map ``y_pred`` to strictly positive domain.

print("Ridge evaluation:")
score_estimator(ridge, df_test)

##############################################################################
#
# Next we fit the Poisson regressor on the target variable,

poisson = make_pipeline(
    column_trans,
    PoissonRegressor(alpha=1/df_train.shape[0], max_iter=1000)
)
poisson.fit(df_train, df_train["Frequency"],
            poissonregressor__sample_weight=df_train["Exposure"])

print("PoissonRegressor evaluation:")
score_estimator(poisson, df_test)

##############################################################################
#
# Finally, we will consider a non-linear model, namely a random forest. Random
# forests do not require the categorical data to be one-hot encoded, instead
# we encode each category label with an arbirtrary integer using
# :class:`preprocessing.OrdinalEncoder` to make the model faster to train (the
# same information is encoded with a small number of features than with
# one-hot encoding).

rf = make_pipeline(
    ColumnTransformer(
        [
            (
                "Veh_Brand_Gas_Region", OrdinalEncoder(),
                ["VehBrand", "VehPower", "VehGas", "Region", "Area"],
            ),
            (
                "Continious", "passthrough",
                ["VehAge", "DrivAge", "BonusMalus", "Density"]
            ),
        ],
        remainder="drop",
    ),
    RandomForestRegressor(min_weight_fraction_leaf=0.01, n_jobs=2)
)
rf.fit(df_train, df_train["Frequency"].values,
       randomforestregressor__sample_weight=df_train["Exposure"].values)


print("RandomForestRegressor evaluation:")
score_estimator(rf, df_test)


##############################################################################
#
# The random forest model also minimizes the conditional least square error.
# However because of a higher predictive power it also results in a smaller
# Poisson deviance than the Poisson regression model.
#
# Not that Evaluating models with a single train / test split is prone to
# random fluctuations. We can verify that we would also get equivalent
# conclusions with cross-validated performance metrics.
#
# The qualitative difference between these models can also be visualized by
# comparing the histogram of observed target values with that of predicted
# values,


fig, axes = plt.subplots(1, 4, figsize=(16, 3))
fig.subplots_adjust(bottom=0.2)
n_bins = 20
df_train["Frequency"].hist(bins=np.linspace(-1, 10, n_bins), ax=axes[0])

axes[0].set_title("Data")
axes[0].set_yscale('log')
axes[0].set_xlabel("y (observed Frequency)")

for idx, model in enumerate([ridge, poisson, rf]):
    y_pred = model.predict(df_train)

    pd.Series(y_pred).hist(bins=np.linspace(-1, 4, n_bins), ax=axes[idx+1])
    axes[idx + 1].set(
        title=model[-1].__class__.__name__,
        yscale='log',
        xlabel="y_pred (predicted expected Frequency)"
    )

##############################################################################
#
# The experimental data presents a long tail distribution for ``y``. In all
# models we predict the mean expected value, so we will have necessairily
# fewer extreme values. Additionally normal distribution used in ``Ridge`` and
# ``RandomForestRegressor`` has a constant variance, while for the Poisson
# distribution used in ``PoissonRegressor``, the variance is proportional to
# the mean predicted value.
#
# Thus, among the considered estimators, ``PoissonRegressor`` is better suited
# for modeling the long tail distribution of the data as compared to the
# ``Ridge`` and ``RandomForestRegressor`` estimators.
#
# To ensure that estimators yield reasonable predictions for different
# policyholder types, we can bin test samples according to `y_pred` returned
# by each model. Then for each bin, compare the mean predicted `y_pred`, with
# the mean observed target:


def _mean_frequency_by_risk_group(y_true, y_pred, sample_weight=None,
                                  n_bins=100):
    """Compare predictions and observations for bins ordered by y_pred

    We order the samples by ``y_pred`` and split it in bins.
    In each bin the observed mean is compared with the predicted
    mean.

    Parameters
    ----------
    y_true: array-like of shape (n_samples,)
        Ground truth (correct) target values.
    y_pred: array-like of shape (n_samples,)
        Estimated target values.
    sample_weight : array-like of shape (n_samples,)
        Sample weights.
    n_bins: int
        number of bins to use

    Returns
    -------
    bin_centers: ndarray of shape (n_bins,)
        bin centers
    y_true_bin: ndarray of shape (n_bins,)
        average y_pred for each bin
    y_pred_bin: ndarray of shape (n_bins,)
        average y_pred for each bin
    """
    idx_sort = np.argsort(y_pred)
    bin_centers = np.arange(0, 1, 1/n_bins) + 0.5/n_bins
    y_pred_bin = np.zeros(n_bins)
    y_true_bin = np.zeros(n_bins)

    for n, sl in enumerate(gen_even_slices(len(y_true), n_bins)):
        weights = sample_weight[idx_sort][sl]
        y_pred_bin[n] = np.average(
               y_pred[idx_sort][sl], weights=weights
        )
        y_true_bin[n] = np.average(
            y_true[idx_sort][sl],
            weights=weights
        )
    return bin_centers, y_true_bin, y_pred_bin


fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 3.2))
plt.subplots_adjust(wspace=0.3)

for axi, (label, model) in zip(ax, [
        ('Ridge', ridge),
        ('PoissonRegressor', poisson),
        ('Random Forest', rf)
]):
    y_pred = model.predict(df_test)

    q, y_true_seg, y_pred_seg = _mean_frequency_by_risk_group(
        df_test["Frequency"].values,
        y_pred,
        sample_weight=df_test["Exposure"].values,
        n_bins=5)

    axi.plot(q, y_pred_seg, marker='o', linestyle="-", label="predictions")
    axi.plot(q, y_true_seg, marker='x', linestyle="--", label="observations")
    axi.set_xlim(0, 1.0)
    axi.set_ylim(0, 0.3)
    axi.set(
        title=label,
        xlabel='Fraction of samples sorted by y_pred',
        ylabel='Mean Frequency (y_pred)'

    )
    axi.legend()


##############################################################################
#
# On the above figure, ``PoissonRegressor`` is the model which presents the
# best consistency between predicted and observed targets, both for low and
# high target values.
#
# The ridge regression model tends to predict very low expected frequencies
# that do not match the data.
#
# The random forest regression model also tends to exaggerate low predicted
# frequencies although to a lower extent than ridge. It also tends to
# exaggerate high frequencies on the other hand.
#
# However for some business applications we are not necessarily interested in
# the the ability of the model in predicting the expected frequency value but
# instead in predicting which customer profiles are the riskiest and which are
# the safest. In this case the model evaluation would cast the problem as a
# ranking problem rather than a regression problem.
#
# To compare the 3 models under this light on, one can plot the fraction
# of cumulated number of claims vs the fraction of cumulated of exposure
# for test samples ordered by the model predictions, from riskiest to safest
# according to each model:


def _cumulated_claims(y_true, y_pred, exposure):
    idx_sort = np.argsort(y_pred)[::-1]  # from riskiest to safest
    sorted_exposure = exposure[idx_sort]
    sorted_frequencies = y_true[idx_sort]
    cumulated_exposure = np.cumsum(sorted_exposure)
    cumulated_exposure /= cumulated_exposure[-1]
    cumulated_claims = np.cumsum(sorted_exposure * sorted_frequencies)
    cumulated_claims /= cumulated_claims[-1]
    return cumulated_exposure, cumulated_claims


fig, ax = plt.subplots(figsize=(8, 8))
plt.subplots_adjust(wspace=0.3)

for (label, model) in [
        ('Ridge', ridge),
        ('PoissonRegressor', poisson),
        ('Random Forest', rf)
]:
    y_pred = model.predict(df_test)
    cum_exposure, cum_claims = _cumulated_claims(
        df_test["Frequency"].values,
        y_pred,
        df_test["Exposure"].values)
    ax.plot(cum_exposure, cum_claims, linestyle="-", label=label)

# Oracle model: y_pred == y_test
cum_exposure, cum_claims = _cumulated_claims(
    df_test["Frequency"].values,
    df_test["Frequency"].values,
    df_test["Exposure"].values)
ax.plot(cum_exposure, cum_claims, linestyle="-.", color="gray", label="Oracle")

# Random Baseline
ax.plot([0, 1], [0, 1], linestyle="--", color="black", label="Random baseline")
ax.set(
    title="Cumulated claims by model",
    xlabel='Fraction of cumulated exposure (from riskiest to safest)',
    ylabel='Fraction of cumulated number of claims'
)
ax.legend()

##############################################################################
#
# This plot reveals that the random forest model is almost uniformly the best
# at sorting customers by risk profiles even if the absolute value of the
# predicted expected frequencies are less well calibrated than for the linear
# Poisson model.
#
#
# All three models are significantly better than chance but also very far from
# making perfect predictions.
#
# This last point is expected due to the nature of the problem: the occurence
# of accidents is mostly dominated by environmental causes that are not
# captured in the columns of the dataset.

plt.show()
