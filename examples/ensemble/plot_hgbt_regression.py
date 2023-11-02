"""
=============================================================
Usecase of advanced features in Histogram Boosting Regression
=============================================================

:ref:`histogram_based_gradient_boosting` (HGBT) models may be the most useful
supervised learning models in scikit-learn. They are based on a modern gradient
boosting implementation comparable to LightGBM and XGBoost. As such, HGBT models
are more feature rich than -and often outperforms- alternative models like
random forests, especially when the number of samples is larger than tens of
thousands of samples (see
:ref:`sphx_glr_auto_examples_ensemble_plot_forest_hist_grad_boosting_comparison.py`).

The top usability features of HGBT models are:

1. :ref:`categorical_support_gbdt` (see
  :ref:`sphx_glr_auto_examples_ensemble_plot_gradient_boosting_categorical.py`).
1. Early stopping.
1. :ref:`nan_support_hgbt`, which avoids the need for an imputer.
1. Support for several losses such as the :ref:`Quantile loss
   <quantile_support_hgbdt>`.
1. :ref:`monotonic_cst_gbdt`.
1. :ref:`interaction_cst_hgbt`.

This example aims at showcasing points 2-5 in a real life setting.
"""


# Author: Arturo Amor <david-arturo.amor-quiroz@inria.fr>
# License: BSD 3 clause

# %%
# Preparing the data
# ==================
# The `electricity dataset <http://www.openml.org/d/151>`_ consists of data
# collected from the Australian New South Wales Electricity Market. In this
# market, prices are not fixed and are affected by supply and demand. They are
# set every five minutes. Electricity transfers to/from the neighboring state of
# Victoria were done to alleviate fluctuations.
#
# The dataset (originally named ELEC2) contains 45,312 instances dated from 7
# May 1996 to 5 December 1998. Each example of the dataset refers to a period of
# 30 minutes, i.e. there are 48 instances for each time period of one day. Each
# example on the dataset has 5 fields: the day of week, the time stamp, the New
# South Wales electricity demand, the Victoria electricity demand. It is
# originally a classification task, but here we use it as a regression where the
# target is the scheduled electricity transfer between states.

from sklearn.datasets import fetch_openml

electricity = fetch_openml(
    name="electricity", version=1, as_frame=True, parser="pandas"
)
df = electricity.frame
X = df.drop(columns=["transfer", "class"])
y = df["transfer"]
X

# %%
# Let us explore the hourly electricity transfer over different days of the week:

import matplotlib.pyplot as plt
import seaborn as sns

colors = sns.color_palette("colorblind")

fig, ax = plt.subplots(figsize=(15, 10))
pointplot = sns.pointplot(x=df["period"], y=df["transfer"], hue=df["day"], ax=ax)
handles, lables = ax.get_legend_handles_labels()
ax.set(
    title="Hourly energy transfer for different days of the week",
    xticks=[i * 2 for i in range(24)],
    xticklabels=list(range(24)),
    xlabel="Time of the day",
    ylabel="Normalized energy transfer",
)
_ = ax.legend(handles, ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"])

# %%
# Notice energy transfer increases systematically during weekends.
#
# Effect of number of trees and early stopping
# ============================================
# For the sake of illustrating the effect of the (maximum) number of trees, we
# train a :class:`~sklearn.ensemble.HistGradientBoostingRegressor` over the
# daily electricity transfer using the whole dataset. Then we visualize its
# predictions depending on the `max_iter` parameter.

from sklearn.ensemble import HistGradientBoostingRegressor

max_iter_list = [10, 50]

fig, ax = plt.subplots(figsize=(10, 5))
average_week_demand = df.groupby(["day", "period"], observed=False)["transfer"].mean()
average_week_demand.plot(color=colors[0], label="training data", linewidth=2, ax=ax)

for idx, max_iter in enumerate(max_iter_list):
    hgbt = HistGradientBoostingRegressor(max_iter=max_iter)
    hgbt.fit(X, y)
    y_pred = hgbt.predict(X)
    prediction_df = df.copy()
    prediction_df["y_pred"] = y_pred
    average_pred = prediction_df.groupby(["day", "period"], observed=False)[
        "y_pred"
    ].mean()
    average_pred.plot(
        color=colors[idx + 1], label=f"max_iter={max_iter}", linewidth=2, ax=ax
    )
ax.set(
    title="Average daily energy transfer during the week",
    xticks=[(i + 0.2) * 48 for i in range(7)],
    xticklabels=["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"],
    xlabel="Time of the week",
    ylabel="Normalized energy transfer",
)
_ = ax.legend()

# %%
# With just a few iterations, HGBT models can achieve convergence (see
# :ref:`sphx_glr_auto_examples_ensemble_plot_forest_hist_grad_boosting_comparison.py`).
#
# Instead of relying solely on `max_iter` to determine when to stop, the HGBT
# implementations in scikit-learn support early stopping. With it, the model
# uses a fraction of the training data as a validation set
# (`validation_fraction`) and stops training if the validation score does not
# improve (or degrades) after `n_iter_no_change` iterations up to a certain
# `tol`.
#
# Notice that there is a trade-off between `learning_rate` and `max_iter`:
# Generally, smaller learning rates require more iterations to converge to the
# minimum loss, while larger learning rates might converge faster but are at
# risk of overfitting.
#
# Indeed, a good practice is to tune the learning rate along with any other
# hyperparameters, fit the HBGT on the training set with a large enough value
# for `max_iter` and determine the best `max_iter` via early stopping and some
# explicit `validation_fraction`.

common_params = {
    "max_iter": 1_000,
    "learning_rate": 0.3,
    "validation_fraction": 0.2,
    "random_state": 42,
    "scoring": "neg_root_mean_squared_error",
}

hgbt = HistGradientBoostingRegressor(early_stopping=True, **common_params)
hgbt.fit(X, y)
plt.plot(-hgbt.validation_score_)
plt.xlabel("number of iterations")
plt.ylabel("root mean squared error")
_ = plt.title(f"Loss of hgbt with early stopping (n_iter={hgbt.n_iter_})")

# %%
# We can then overwrite the value for `max_iter` to a reasonable value and avoid
# the extra computational cost of the inner validation. In this case, rounding
# up the number of iterations to 600 may account for variability of the training
# set:

common_params["max_iter"] = 600
common_params["early_stopping"] = False
hgbt = HistGradientBoostingRegressor(**common_params)

# %%
# .. note:: The inner validation done during early stopping is not optimal for
#    time series with the implementation as of scikit-learn v1.3.
#
# Support for missing values
# ==========================
# HGBT models have native support of missing values. During training, the tree
# grower decides where samples with missing values should go (left or right
# child) at each split, based on the potential gain. When predicting, these
# samples are sent to either child accordingly. If a feature had no missing
# values during training, samples with missing values for that feature are sent
# to the child with the most samples.
#
# The present example shows how HGBT regressions deal with values missing
# completely at random (MCAR), i.e. the missingness does not depend on the
# observed data or the unobserved data. We can simulate such scenario by
# randomly replacing values from randomly selected features with `Nan` values.

import numpy as np

from sklearn.model_selection import TimeSeriesSplit, cross_validate

np.random.seed(42)

ts_cv = TimeSeriesSplit(n_splits=5, gap=48, max_train_size=10000, test_size=1000)
train_0, test_0 = next(ts_cv.split(df))
last_days = slice(-192, None)
total_cells = X.shape[0] * X.shape[1]
missing_fraction_list = [0, 0.01, 0.03]

fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(y.iloc[test_0].values[last_days], label="Actual transfer")

for missing_fraction in missing_fraction_list:
    num_missing_cells = int(total_cells * missing_fraction)
    row_indices = np.random.choice(X.shape[0], num_missing_cells, replace=True)
    col_indices = np.random.choice(X.shape[1], num_missing_cells, replace=True)
    X = df.drop(columns=["transfer", "class"])
    X.iloc[row_indices, col_indices] = np.nan

    hgbt.fit(X.iloc[train_0], y.iloc[train_0])
    hgbt_predictions = hgbt.predict(X.iloc[test_0])
    cv_results = cross_validate(
        hgbt,
        X,
        y,
        cv=ts_cv,
        scoring="neg_root_mean_squared_error",
    )
    rmse = -cv_results["test_score"]
    ax.plot(
        hgbt_predictions[last_days],
        label=(
            f"missing_fraction={missing_fraction}, RMSE={rmse.mean():.2f} +/-"
            f" {rmse.std():.2f}"
        ),
        alpha=0.5,
    )
ax.set(
    title="Daily energy transfer predictions on data with MCAR values",
    xticks=[(i + 0.25) * 48 for i in range(4)],
    xticklabels=["Tue", "Wed", "Thu", "Fri"],
    xlabel="Time of the week",
    ylabel="Normalized energy transfer",
)
_ = ax.legend()

# %%
# Support for quantile loss
# =========================
#
# The quantile loss in regression enables a view of the variability or
# uncertainty of the target variable. For instance, predicting the 5th and 95th
# percentiles can provide a 90% prediction interval, i.e. the range within which
# we expect the true value to fall with 90% probability.

from sklearn.metrics import make_scorer, mean_pinball_loss

quantiles = [0.95, 0.05]
predictions = []
X = df.drop(columns=["transfer", "class"])  # reset X

fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(y.iloc[test_0].values[last_days], label="Actual transfer")

for quantile in quantiles:
    hgbt_quantile = HistGradientBoostingRegressor(
        loss="quantile", quantile=quantile, **common_params
    )
    hgbt_quantile.fit(X.iloc[train_0], y.iloc[train_0])
    hgbt_predictions = hgbt_quantile.predict(X.iloc[test_0])

    predictions.append(hgbt_predictions)
    cv_results = cross_validate(
        hgbt_quantile,
        X,
        y,
        cv=ts_cv,
        scoring=make_scorer(mean_pinball_loss, alpha=quantile),
    )
    score = cv_results["test_score"]
    ax.plot(
        hgbt_predictions[last_days],
        label=(
            f"quantile={quantile}, pinball loss={score.mean():.3f} +/-"
            f" {score.std():.3f}"
        ),
        alpha=0.5,
    )

ax.fill_between(
    range(len(predictions[0][last_days])),
    predictions[0][last_days],
    predictions[1][last_days],
    color=colors[0],
    alpha=0.1,
)
ax.set(
    title="Daily energy transfer predictions with quantile loss",
    xticks=[(i + 0.25) * 48 for i in range(4)],
    xticklabels=["Tue", "Wed", "Thu", "Fri"],
    xlabel="Time of the week",
    ylabel="Normalized energy transfer",
)
_ = ax.legend()

# %%
# Keep in mind that those predicted percentiles are just estimations from a
# model. One can still improve the quality of such estimations by:
#
# - collecting more data-points;
# - better tuning of the model hyperparameters (see
#   :ref:`sphx_glr_auto_examples_ensemble_plot_gradient_boosting_quantile.py`);
# - engineering more predictive features from the same data (see
#   :ref:`sphx_glr_auto_examples_applications_plot_cyclical_feature_engineering.py`).
#
# Monotonic constraints
# =====================
#
# Given specific domain knowledge that requires the relationship between a
# feature and the target to be monotonically increasing or decreasing, one can
# enforce such behaviour in the predictions of a HGBT model using monotonic
# constraints. This makes the model more interpretable and can reduce its
# variance (and potentially mitigate overfitting) at the risk of increasing
# bias. Monotonic constraints can also be used to enforce specific regulatory
# requirements, ensure compliance and align with ethical considerations.
#
# In the present example, the policy of transfering energy from Victoria to New
# South Wales is meant to alleviate price fluctuations, meaning that the model
# predictions have to enforce such goal, i.e. transfer should increase with
# price and demand in New South Wales, but also decrease with price and demand
# in Victoria, in order to benefit both populations.
#
# If the training data has feature names, it’s possible to specify the monotonic
# constraints by passing a dictionary with the convention:
#
# - 1: monotonic increase
# - 0: no constraint
# - -1: monotonic decrease
#
# Else, one can pass an array-like encoding the above convention by position.

from sklearn.inspection import PartialDependenceDisplay

monotonic_cst = {
    "date": 0,
    "day": 0,
    "period": 0,
    "nswdemand": 1,
    "nswprice": 1,
    "vicdemand": -1,
    "vicprice": -1,
}
hgbt_no_cst = HistGradientBoostingRegressor().fit(X, y)
hgbt_cst = HistGradientBoostingRegressor(monotonic_cst=monotonic_cst).fit(X, y)

fig, ax = plt.subplots(nrows=2, figsize=(15, 10))
disp = PartialDependenceDisplay.from_estimator(
    hgbt_no_cst,
    X,
    features=["nswdemand", "nswprice"],
    line_kw={"linewidth": 2, "label": "unconstrained", "color": "tab:blue"},
    ax=ax[0],
)
PartialDependenceDisplay.from_estimator(
    hgbt_cst,
    X,
    features=["nswdemand", "nswprice"],
    line_kw={"linewidth": 2, "label": "constrained", "color": "tab:orange"},
    ax=disp.axes_,
)
disp = PartialDependenceDisplay.from_estimator(
    hgbt_no_cst,
    X,
    features=["vicdemand", "vicprice"],
    line_kw={"linewidth": 2, "label": "unconstrained", "color": "tab:blue"},
    ax=ax[1],
)
PartialDependenceDisplay.from_estimator(
    hgbt_cst,
    X,
    features=["vicdemand", "vicprice"],
    line_kw={"linewidth": 2, "label": "constrained", "color": "tab:orange"},
    ax=disp.axes_,
)
_ = plt.legend()

# %%
# Indeed, we can verify that the predictive quality of the model is not
# significantly degraded by introducing the monotonic constraints:

cv_results = cross_validate(
    hgbt_no_cst,
    X,
    y,
    cv=ts_cv,
    scoring="neg_root_mean_squared_error",
)
rmse = -cv_results["test_score"]
print(f"RMSE without constraints = {rmse.mean():.2f} +/- {rmse.std():.2f}")

cv_results = cross_validate(
    hgbt_cst,
    X,
    y,
    cv=ts_cv,
    scoring="neg_root_mean_squared_error",
)
rmse = -cv_results["test_score"]
print(f"RMSE with constraints    = {rmse.mean():.2f} +/- {rmse.std():.2f}")

# %%
# That being said, notice the comparison is between to different models that may
# be optimized by a different combination of hyperparameters. That is the reason
# why we do no use the `common_params` in this section as done before.
