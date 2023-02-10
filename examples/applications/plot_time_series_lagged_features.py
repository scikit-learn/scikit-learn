"""
===========================================
Lagged Features for Time Series Forecasting
===========================================

This example demonstrates how pandas-engineered lagged features can be used
for time series forecasting with
:class:`~sklearn.ensemble.HistGradientBoostingRegressor` on the Bike Sharing
Demand dataset.

See the example on
:ref:`sphx_glr_auto_examples_applications_plot_cyclical_feature_engineering.py`
for some data exploration on this dataset and a demo on periodic feature
engineering.

"""

# %%
# Analyzing the Bike Sharing Demand Dataset
# -----------------------------------------
#
# We start by loading the data from the OpenML repository.
from sklearn.datasets import fetch_openml
import numpy as np
import pandas as pd

bike_sharing = fetch_openml(
    "Bike_Sharing_Demand", version=2, as_frame=True, parser="pandas"
)
df = bike_sharing.frame
df

# %%
# Next, we take a look at the statistical summary of the dataset
# so that we can better understand the data that we are working with.
summary = pd.DataFrame(df.describe())
summary = (
    summary.style.background_gradient()
    .set_table_attributes("style = 'display: inline'")
    .set_caption("Statistics of the Dataset")
    .set_table_styles([{"selector": "caption", "props": [("font-size", "16px")]}])
)
summary

# %%
# Taking a look at the total count of rented bikes for the first
# week in our dataset.
import matplotlib.pyplot as plt
import seaborn as sns

bike_count = df["count"]
bike_count[: 7 * 24].plot(figsize=(15, 6))
_ = plt.title("First week of bike rental count data")

# %%
# Next, let us look at the fraction of the rented fleet vs the number
# number of hours.
fig, ax = plt.subplots(figsize=(12, 4))
bike_count.hist(bins=30, ax=ax)
_ = ax.set(
    xlabel="Fraction of rented fleet demand",
    ylabel="Number of hours",
)

# %%
# Let us look at the count of the seasons `fall`, `spring`, `summer`
# and `winter` present in the dataset.
plt.figure(figsize=(15, 10))
ax = sns.countplot(data=df, x="season")

bbox_args = dict(boxstyle="round", fc="0.9")
for p in ax.patches:
    ax.annotate(
        "{:.0f}".format(p.get_height()),
        (p.get_x() + 0.3, p.get_height() + 75),
        color="black",
        bbox=bbox_args,
        fontsize=15,
    )

plt.title("Count for Season", fontsize=20)
plt.xlabel("Season", fontsize=15)
plt.ylabel("Count", fontsize=15)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.show()

# %%
# To check for outliers in the features of the dataset we create box-plots
# for `temp`, `feel_temp`, `humidity` and `windspeed`.
plt.figure(figsize=(15, 15))


def create_boxplot(feature, color):
    sns.boxplot(data=df[[feature]], color=color)
    plt.title("Box Plot for " + feature.title(), fontsize=20)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)


plt.subplot(221)
create_boxplot("temp", "blue")

plt.subplot(222)
create_boxplot("feel_temp", "orange")

plt.subplot(223)
create_boxplot("humidity", "green")

plt.subplot(224)
create_boxplot("windspeed", "red")

# %%
# Exploring the hourly count of bikes over weekdays:
plt.figure(figsize=(15, 10))
sns.pointplot(x=df["hour"].astype(int), y=df["count"], hue=df["weekday"].astype(int))
plt.title("Hourly Count for Weekdays", fontsize=20)
plt.xlabel("Hour", fontsize=15)
plt.ylabel("Count", fontsize=15)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.show()

# %%
# Exploring the hourly count of bikes over different seasons:
plt.figure(figsize=(15, 10))
sns.pointplot(x=df["hour"].astype(int), y=df["count"], hue=df["season"])
plt.title("Hourly Count for Seasons", fontsize=20)
plt.xlabel("Hour", fontsize=15)
plt.ylabel("Count", fontsize=15)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.show()

# %%
# Generating Pandas-Engineered Lagged Features
# --------------------------------------------
# Let's consider the problem of predicting the demand at the
# next hour given past demands. Since the demand is a continuous
# variable, one could intuitively use any regression model. However, we do
# not have the usual `(X_train, y_train)` dataset. Instead, we just have
# the `y_train` demand data sequentially organized by time.
count = bike_count
lagged_df = pd.concat(
    [
        count,
        count.shift(1).rename("lagged_count_1h"),
        count.shift(2).rename("lagged_count_2h"),
        count.shift(3).rename("lagged_count_3h"),
        count.shift(24).rename("lagged_count_1d"),
        count.shift(24 + 1).rename("lagged_count_1d_1h"),
        count.shift(7 * 24).rename("lagged_count_7d"),
        count.shift(7 * 24 + 1).rename("lagged_count_7d_1h"),
        count.shift(1).rolling(24).mean().rename("lagged_mean_24h"),
        count.shift(1).rolling(24).max().rename("lagged_max_24h"),
        count.shift(1).rolling(24).min().rename("lagged_min_24h"),
        count.shift(1).rolling(7 * 24).mean().rename("lagged_mean_7d"),
        count.shift(1).rolling(7 * 24).max().rename("lagged_max_7d"),
        count.shift(1).rolling(7 * 24).min().rename("lagged_min_7d"),
    ],
    axis="columns",
)
lagged_df.tail(10)

# %%
# Watch out however, the first lines have undefined values because their own
# past is unknown. This depends on how much lag we used:
lagged_df.head(10)

# %%
# We can now separate the lagged features in a matrix `X` and the target variable
# (the counts to predict) in an array of the same first dimension `y`.
lagged_df = lagged_df.dropna()
X = lagged_df.drop("count", axis="columns")
y = lagged_df["count"]
print("X shape: {}\ny shape: {}".format(X.shape, y.shape))

# %%
# Naive evaluation of the next hour bike demand regression
# --------------------------------------------------------
# Let's randomly split our tabularized dataset to train a gradient
# boosting regression tree (GBRT) model and evaluate it using Mean
# Absolute Percentage Error (MAPE).
from sklearn.model_selection import train_test_split
from sklearn.ensemble import HistGradientBoostingRegressor

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

model = HistGradientBoostingRegressor().fit(X_train, y_train)

# %%
# Taking a look at the performance of the model.
from sklearn.metrics import mean_absolute_percentage_error

y_pred = model.predict(X_test)
mean_absolute_percentage_error(y_test, y_pred)

# %%
# Proper Next Hour Forecasting Evaluation
# ---------------------------------------
# Let's use a proper evaluation splitting strategies that takes into account
# the temporal structure of the dataset to evaluate our model's ability to
# predict data points in the future (to avoid cheating by reading values from
# the lagged features in the training set).
from sklearn.model_selection import TimeSeriesSplit

ts_cv = TimeSeriesSplit(
    n_splits=3,  # to keep the notebook fast enough on common laptops
    gap=48,  # 2 days data gap between train and test
    max_train_size=10000,  # keep train sets of comparable sizes
    test_size=3000,  # for 2 or 3 digits of precision in scores
)
all_splits = list(ts_cv.split(X, y))

# %%
# Training the model and evaluating its performance based on MAPE.
train_idx, test_idx = all_splits[0]
X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

model = HistGradientBoostingRegressor().fit(X_train, y_train)
y_pred = model.predict(X_test)
mean_absolute_percentage_error(y_test, y_pred)

# %%
# The error rate of this model is better than our naive shuffling
# train-test split. This is quite expected but maybe the first split
# is easier to predict (more regular) than the others. Let's assess this
# variability of our error evaluation with proper cross-validation:
from sklearn.model_selection import cross_val_score

cv_mape_scores = -cross_val_score(
    model, X, y, cv=ts_cv, scoring="neg_mean_absolute_percentage_error"
)
cv_mape_scores

# %%
# The variability across splits is quite large! In a real life setting
# it would be advised to use more splits to better assess the variability.
# Let's report the mean CV scores and their standard deviation from now on.
print(f"CV MAPE: {cv_mape_scores.mean():.3f} ± {cv_mape_scores.std():.3f}")

# %%
# To get a finer evaluation of our models we can compute and report several
# cross-validation metrics at once using a dedicated helper function:
from time import perf_counter
from sklearn.model_selection import cross_validate
from sklearn.metrics import mean_pinball_loss
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import mean_squared_error


def evaluate(model, X, y, cv):
    def score_func(estimator, X, y):
        y_pred = estimator.predict(X)
        return {
            "mean_absolute_percentage_error": mean_absolute_percentage_error(y, y_pred),
            "root_mean_squared_error": np.sqrt(mean_squared_error(y, y_pred)),
            "mean_absolute_error": mean_absolute_error(y, y_pred),
            "mean_pinball_05_loss": mean_pinball_loss(y, y_pred, alpha=0.05),
            "mean_pinball_50_loss": mean_pinball_loss(y, y_pred, alpha=0.50),
            "mean_pinball_95_loss": mean_pinball_loss(y, y_pred, alpha=0.95),
        }

    tic = perf_counter()
    cv_results = cross_validate(
        model,
        X,
        y,
        cv=cv,
        scoring=score_func,
    )
    toc = perf_counter()
    for key, value in cv_results.items():
        if key.startswith("test_"):
            print(f"{key[5:]}: {value.mean():.3f} ± {value.std():.3f}")
    print(f"\ndone in {toc - tic:.3f} s")


gbrt_mse = HistGradientBoostingRegressor(loss="squared_error")
evaluate(gbrt_mse, X, y, cv=ts_cv)

# %%
# Model evaluation using `loss="poisson"`
gbrt_poisson = HistGradientBoostingRegressor(loss="poisson")
evaluate(gbrt_poisson, X, y, cv=ts_cv)

# %%
# Modeling Predictive Uncertainty via Quantile Regression
# -------------------------------------------------------
# Instead of modeling the expected value of the distribution of
# `Y|X` like the least squares and Poisson losses do, one could try to
# estimate quantiles of the conditional distribution.
#
# `Y|X=xi` is expected to be a random variable for a given data point `xi`
# because we expect that the number of rentals cannot be 100% accurately
# predicted from the features. It can be influenced by other variables
# not properly captured by the existing lagged features. For instance
# whether or not it will rain in the next hour cannot be fully be anticipated
# from the past hours bike rental data. This is what we call aleatoric
# uncertainty.
#
# Quantile regression makes it possible to give a finer description of that
# distribution without making strong assumptions on its shape.
#
# The conditional 5th percentile (a.k.a. 0.05-quantile) can be estimated with:
gbrt_percentile_05 = HistGradientBoostingRegressor(loss="quantile", quantile=0.05)
evaluate(gbrt_percentile_05, X, y, cv=ts_cv)

# %%
# The conditional median (0.50-quantile) can be estimated with:
gbrt_median = HistGradientBoostingRegressor(loss="quantile", quantile=0.5)
evaluate(gbrt_median, X, y, cv=ts_cv)

# %%
# And finally the 0.95 quantile:
gbrt_percentile_95 = HistGradientBoostingRegressor(loss="quantile", quantile=0.95)
evaluate(gbrt_percentile_95, X, y, cv=ts_cv)

# %%
# A Qualitative Look at the Predictions
# -------------------------------------
# We can now visualize the performance of the model with regards
# to the 5th percentile, median and the 95th percentile:
all_splits = list(ts_cv.split(X, y))
train_idx, test_idx = all_splits[0]

X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

max_iter = 50
gbrt_mean_poisson = HistGradientBoostingRegressor(loss="poisson", max_iter=max_iter)
gbrt_mean_poisson.fit(X_train, y_train)
mean_predictions = gbrt_mean_poisson.predict(X_test)

gbrt_median = HistGradientBoostingRegressor(
    loss="quantile", quantile=0.5, max_iter=max_iter
)
gbrt_median.fit(X_train, y_train)
median_predictions = gbrt_median.predict(X_test)

gbrt_percentile_5 = HistGradientBoostingRegressor(
    loss="quantile", quantile=0.05, max_iter=max_iter
)
gbrt_percentile_5.fit(X_train, y_train)
percentile_5_predictions = gbrt_percentile_5.predict(X_test)

gbrt_percentile_95 = HistGradientBoostingRegressor(
    loss="quantile", quantile=0.95, max_iter=max_iter
)
gbrt_percentile_95.fit(X_train, y_train)
percentile_95_predictions = gbrt_percentile_95.predict(X_test)

# %%
last_hours = slice(-96, None)
fig, ax = plt.subplots(figsize=(15, 7))
plt.title("Predictions by regression models")
ax.plot(
    y_test.values[last_hours],
    "x-",
    alpha=0.2,
    label="Actual demand",
    color="black",
)
ax.plot(
    median_predictions[last_hours],
    "^-",
    label="GBRT median",
)
ax.plot(
    mean_predictions[last_hours],
    "x-",
    label="GBRT mean (Poisson)",
)
ax.fill_between(
    np.arange(96),
    percentile_5_predictions[last_hours],
    percentile_95_predictions[last_hours],
    alpha=0.3,
    label="GBRT 90% interval",
)
_ = ax.legend()

# %%
# Looking at the performance of non-linear regression models vs
# perfect models:
fig, axes = plt.subplots(ncols=3, figsize=(15, 6), sharey=True)
fig.suptitle("Non-linear regression models")
predictions = [
    median_predictions,
    percentile_5_predictions,
    percentile_95_predictions,
]
labels = [
    "Median",
    "5th percentile",
    "95th percentile",
]
for ax, pred, label in zip(axes, predictions, labels):
    ax.scatter(pred, y_test.values, alpha=0.3, label=label)
    ax.plot([0, y.max()], [0, y.max()], "--", label="Perfect model")
    ax.set(
        xlim=(0, y.max()),
        ylim=(0, y.max()),
        xlabel="Predicted demand",
        ylabel="True demand",
    )
    ax.legend()

plt.show()

# %%
# Concluding Remarks
# ------------------
# Through this example we explored time series forecasting using lagged
# features. We compared a naive regression (using the standardized
# `train_test_split`) with a proper time series evaluation strategy using
# `TimeSeriesSplit`. We observed that the model trained using `TimeSerieSplit`
# gives a better Mean Average Percentage Error (MAPE) rate than the naive
# method. We also analyzed the predictive uncertainty of our model via
# Quantile Regression. Predictions based on the 5th, 50th and 95th
# percentile using `loss="quantile"` provide us with a better understanding
# of the forecasts made by our time series regression model.
