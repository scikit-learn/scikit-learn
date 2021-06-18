
"""
==========================
Cyclic feature engineering
==========================

This notebook introduces different strategies to leverage time-related features
for a bike sharing demand regression task that is highly dependent on business
cycles (days, weeks, months) and yearly season cycles.

In the process we introduce how to perform periodic feature engineering using
the :class:`sklearn.preprocessing.SplineTransformer` class and its
`extrapolation="periodic"` option.

"""
# %%
# Bike Sharing Demand Data Exploration
# ------------------------------------
#
# We start by loading the data from the OpenML repository.
from sklearn.datasets import fetch_openml

bike_sharing = fetch_openml("Bike_Sharing_Demand", version=2, as_frame=True)
df = bike_sharing.frame

# %%
# To get a quick understanding of the periodic patterns of the data, let us
# have at the average weekly demand.
#
# Note that the week starts on a Sunday, during the week-end. We can clearly
# distinguish the commute patterns in the morning and evenings of the work days
# and the leisure use of the bikes on the weekends with a more spread peak
# demand around the middle of the days:
import matplotlib.pyplot as plt


average_week_demand = df.groupby(["weekday", "hour"]).mean()["count"]
average_week_demand.plot(figsize=(12, 4))
_ = plt.title("Average number of bike rentals during the week")


# %%
#
# The target of the prediction problem is the absolute count of bike rentals on
# a hourly basis:
df["count"].max()

# %% [markdown]
#
# Let's rescale the target variable (absolute value of the demand) to predict a
# relative demand so that the mean absolute error is more easily interpreted as
# a fraction of the maximum demand and avoid numerical issues with least
# squares regression.

# %%
y = df["count"] / 1000
_ = y.hist(bins=30)

# %%
# The input feature data frame is a time annotated hourly log of variables
# describing the weather conditions. It includes both numerical and categorical
# variables. Note that the time information has already been expanded into
# several complementary columns.
#
X = df.drop("count", axis="columns")
X
# %%
# We first introspect the distribution of the categorical variables, starting
# with `"weather"`:
#
X["weather"].value_counts()

# %%
# Since there are only a few `"heavy_rain"` events, we cannot use this category
# to train machine learning models. Instead it we simplify the representation
# by collapsing those into the `"rain"` category.
#
X["weather"][X["weather"] == "heavy_rain"] = "rain"
# %%
X["weather"].value_counts()

# %%
# As expected season variable is well balanced:
#
X["season"].value_counts()

# %%
# Time-based Cross-Validation Splits
# ----------------------------------
#
# Since the dataset is a time-ordered event log (hourly demand), we will use a
# time-sensitive cross-validation splitter to evaluate our demand forecasting
# model as realistically as possible. We use a gap of 2 days between the train
# and test side of the splits. We also limit the training set size to make the
# performance of the CV folds more stable.
#
# 1000 test datapoints should be enough to quantify the performance of the
# model. This reprent a bit less than a month and a half of contiguous test
# data:

from sklearn.model_selection import TimeSeriesSplit

ts_cv = TimeSeriesSplit(
    n_splits=5, gap=48, max_train_size=10000, test_size=1000,
)

# %%
# Let us manually inspect the various splits to check that the
# `TimeSeriesSplit` works as we expect, starting with the first split:
all_splits = list(ts_cv.split(X, y))
train_0, test_0 = all_splits[0]

# %%
X.iloc[test_0]

# %%
X.iloc[train_0]

# %%
# We now inspect the last split:
train_4, test_4 = all_splits[4]

# %%
X.iloc[test_4]

# %%
X.iloc[train_4]

# %%
# All is well. We are now ready to do some predictive modeling!
#
# Gradient Boosting
# -----------------
#
# Gradient Boosting Regression with decision trees is often flexible enough
# to efficiently handle heteorogenous tabular data with a mix of categorical
# and numerical features as long as the number of samples is large enough.
#
# Here we do minimal ordinal encoding for the categorical variables and then
# let the model know that it should treat those a categorical variables using
# a dedicated split rule.
#
# The numerical variable need no preprocessing and for the sake of simplicity
# we only try the default hyper-parameters for this model:
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OrdinalEncoder
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.model_selection import cross_validate


categorical_columns = [
    "weather",
    "season",
    "holiday",
    "workingday",
]
categories = [
    ["clear", "misty", "rain", "heavy_rain"],
    ["spring", "summer", "fall", "winter"],
    ["False", "True"],
    ["False", "True"],
]
ordinal_encoder = OrdinalEncoder(categories=categories)


gbrt_pipeline = make_pipeline(
    ColumnTransformer(
        transformers=[
            ("categorical", ordinal_encoder, categorical_columns),
        ],
        remainder="passthrough",
    ),
    HistGradientBoostingRegressor(
        categorical_features=range(4),
    ),
)

# %%
#
# Lets evaluate our gradient boosting model with the mean absolute error of the
# relative demand averaged accross our 5 time-based cross-validation splits:


def evaluate(model, X, y, cv):
    cv_results = cross_validate(
        model,
        X,
        y,
        cv=ts_cv,
        scoring="neg_mean_absolute_error"
    )
    scores = -cv_results["test_score"]
    print(f"Mean Absolute Error: {scores.mean():.3f} +/- {scores.std():.3f}")


evaluate(gbrt_pipeline, X, y, cv=ts_cv)

# %%
# This models has an average error around 4 to 5% of the maximum demand. This
# is quite good for a first trial without any hyper-parameters tuning! We just
# had to make the categorical variables explicit. Note that the time related
# features are passed as such but this not much of a problem for tree-based
# models as they can learn a non-monotonic relationship between ordinal input
# features and the target.
#
# This is not the case for linear regression model as we will see in the
# following.
#
# Naive Linear Regression
# -----------------------
#
# As usual for linear models, categorical variables need to be one-hot encoded.
# As the numerical features are all approximately on similar scales, we let
# them go through unscaled:
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import RidgeCV
import numpy as np


one_hot_encoder = OneHotEncoder(handle_unknown="ignore")
alphas = np.logspace(-6, 6, 25)
naive_linear_pipeline = make_pipeline(
    ColumnTransformer(
        transformers=[
            ("categorical", one_hot_encoder, categorical_columns),
        ],
        remainder="passthrough",
    ),
    RidgeCV(alphas=alphas),
)


evaluate(naive_linear_pipeline, X, y, cv=ts_cv)

# %%
#
# The performance is a not good: around 14% of the max demand on average. This
# is more than three times higher than the average error of the gradient
# boosting model. We can suspect that the naive original encoding of the
# periodic time-related features might prevent the linear regression model to
# properly leverage the time information: linear regression cannot model a
# non-monotonic relationship between the input features and the target.
#
# Trigonometric Features
# ----------------------
#
# As a first attempt, we can try to encode each of those periodic features
# using a sine and cosine transform with the matching period.
#
# Each ordinal time feature is transformed into 2 features that together encode
# equivalent information in a non-monotonic way, and more importantly without
# jump between the first and the last value of the periodic range.
from sklearn.preprocessing import FunctionTransformer


def sin_transformer(period):
    return FunctionTransformer(lambda x: np.sin(x / period * 2 * np.pi))


def cos_transformer(period):
    return FunctionTransformer(lambda x: np.cos(x / period * 2 * np.pi))


# %%
#
# Let us visualize the effect of this feature expansion on some synthetic hour
# data with a bit of extrapolation beyond hour=23:
import pandas as pd

hour_df = pd.DataFrame(
    np.linspace(0, 26, 1000).reshape(-1, 1),
    columns=["hour"],
)
hour_df["hour_sin"] = sin_transformer(24).fit_transform(hour_df)["hour"]
hour_df["hour_cos"] = cos_transformer(24).fit_transform(hour_df)["hour"]
hour_df.plot(x="hour")
_ = plt.title("Trigonometric encoding for the 'hour' feature")

# %%
#
# We can now build a feature extraction pipeline using this strategy:
cyclic_cossin_transformer = ColumnTransformer(
    transformers=[
        ("categorical", one_hot_encoder, categorical_columns),
        ("month_sin", sin_transformer(12), ["month"]),
        ("month_cos", cos_transformer(12), ["month"]),
        ("weekday_sin", sin_transformer(7), ["weekday"]),
        ("weekday_cos", cos_transformer(7), ["weekday"]),
        ("hour_sin", sin_transformer(24), ["hour"]),
        ("hour_cos", cos_transformer(24), ["hour"]),
    ],
    remainder="passthrough",
)
cyclic_cossin_linear_pipeline = make_pipeline(
    cyclic_cossin_transformer,
    RidgeCV(alphas=alphas),
)
evaluate(cyclic_cossin_linear_pipeline, X, y, cv=ts_cv)


# %%
#
# Unfortunately this simple feature engineering does not seem to drastically
# improve the performance of our linear regression model. We will further
# analyze possible reasons for this disappointing outcome at the end of this
# notebook.
#
# Periodic Spline Features
# ------------------------
#
# Instead we can try an alternative encoding of the periodic time-related
# features using spline transformations with a large-enough number of knots,
# and as a result a larger number of expanded features:
from sklearn.preprocessing import SplineTransformer


def periodic_spline_transformer(period, n_knots=None, degree=3):
    if n_knots is None:
        n_knots = period
    return SplineTransformer(
        degree=degree,
        n_knots=n_knots,
        knots=np.linspace(0, period, n_knots).reshape(n_knots, 1),
        extrapolation="periodic",
    )


# %%
#
# Again, let us visualize the effect of this feature expansion on some synthetic hour
# data with a bit of extrapolation beyond hour=23:
hour_df = pd.DataFrame(
    np.linspace(0, 26, 1000).reshape(-1, 1),
    columns=["hour"],
)
splines = periodic_spline_transformer(24, n_knots=12).fit_transform(hour_df)
splines_df = pd.DataFrame(
    splines,
    columns=[f"spline_{i}" for i in range(splines.shape[1])],
)
pd.concat([hour_df, splines_df], axis="columns").plot(x="hour")
_ = plt.title("Periodic spline-based encoding for the 'hour' feature")


# %%
# We can now build a predictive pipeline using this alternative periodic
# feature engineering strategy.
#
# For the "hours" columns we use only 12 knots for a period of 24 hours to
# avoid over-representing this feature compared to months (12 natural knots)
# and weekday (7 natural knots).
cyclic_spline_transformer = ColumnTransformer(
    transformers=[
        ("categorical", one_hot_encoder, categorical_columns),
        ("cyclic_month", periodic_spline_transformer(12), ["month"]),
        ("cyclic_weekday", periodic_spline_transformer(7), ["weekday"]),
        ("cyclic_hour", periodic_spline_transformer(24, n_knots=12), ["hour"]),
    ],
    remainder="passthrough",
)
cyclic_spline_linear_pipeline = make_pipeline(
    cyclic_spline_transformer,
    RidgeCV(alphas=alphas),
)
evaluate(cyclic_spline_linear_pipeline, X, y, cv=ts_cv)


# %%
# Spline features make it possible for the linear model to successfully
# leverage the periodic time-related features and reduce the error from ~14% to
# ~10% of the maximum demand.
#
# However linear models cannot capture non-linear interactions between features
# even if the features them-selves are marginally non-linear as is the case
# with features constructed by `SplineTransformer`.
#
# An example of a such a non-linear interaction that we would like to model
# could be the impact of the rain that might not be the same during the working
# days and the week-ends and holidays for instance.
#
# To capture such interactions we could either use a partial polynomial
# expansion `PolynomicalFeatures(degree=2, interaction_only=True)` or use the
# Nystroem method to compute an approximate polynomial kernel expansion.
#
# Let us try the latter for computational reasons:
from sklearn.preprocessing import MinMaxScaler
from sklearn.kernel_approximation import Nystroem


cyclic_spline_poly_pipeline = make_pipeline(
    cyclic_spline_transformer,
    MinMaxScaler(),
    Nystroem(kernel="poly", degree=2, n_components=300, random_state=0),
    RidgeCV(alphas=alphas),
)
evaluate(cyclic_spline_poly_pipeline, X, y, cv=ts_cv)

# %%
# We observe that this model performance can almost rival the performance of
# the gradient boosted trees with an average error around 6% of the maximum
# demand.
#
# Note that while the final step of this pipeline is a linear regression model,
# the intermediate steps such as the spline feature extraction and the Nyström
# kernel approximation are highly non-linear. As a result the compound
# pipeline is much more expressive than a simple linear model.
#
# Qualitative Analysis of the Models Predictions
# ----------------------------------------------
#
# Let us have a qualitative look at the predictions of this model and the gradient
# boosted trees:

# %%
gbrt_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
gbrt_predictions = gbrt_pipeline.predict(X.iloc[test_0])
# %%
cyclic_spline_poly_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_spline_poly_predictions = cyclic_spline_poly_pipeline.predict(X.iloc[test_0])

# %%
plt.figure(figsize=(6, 6))
plt.scatter(y.iloc[test_0].values, gbrt_predictions,
            alpha=0.2, label="Gradient Boosted Trees")
plt.scatter(y.iloc[test_0].values, cyclic_spline_poly_predictions,
            alpha=0.2, label="Splines + polynomial kernel regression")
plt.plot([0, 1], [0, 1], "--", label="Perfect model")
plt.xlim(0, 1)
plt.ylim(0, 1)
plt.xlabel("True demand")
plt.ylabel("Predicted demand")
_ = plt.legend()


# %%
# We can observe that both models under-estimate the high demand events, but
# gradient boosting a bit less so. The low demand events are well predicted on
# average by gradient boosting while the spline polynomial regression pipeline
# seems to systematic over-estimate demand in that regime.
#
# We can qualitatively confirm this global analysis by zooming on the last 96
# hours (4 days) of the test set:

last_hours = slice(-96, None)
plt.figure(figsize=(12, 4))
plt.plot(
    gbrt_predictions[last_hours],
    "x-",
    label="Gradient Boosted Trees",
)
plt.plot(
    cyclic_spline_poly_predictions[last_hours],
    "x-",
    label="Splines + polynomial kernel regression",
)
plt.title("Predictions by non-linear models")
plt.plot(
    y.iloc[test_0].values[last_hours],
    "x-",
    alpha=0.2,
    label="Actual demand",
    color="black",
)
_ = plt.legend()

# %%
# Let us know compare the predictions of the 3 linear models (without kernel
# approximation).
#
# Here we want to visualize the impact of the feature engineering choices on
# the time related-shape of the predictions.
naive_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
naive_linear_predictions = naive_linear_pipeline.predict(X.iloc[test_0])

cyclic_cossin_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_cossin_linear_predictions = cyclic_cossin_linear_pipeline.predict(
    X.iloc[test_0]
)

cyclic_spline_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_spline_linear_preictions = cyclic_spline_linear_pipeline.predict(
    X.iloc[test_0]
)

# %%
plt.figure(figsize=(12, 4))
last_hours = slice(-96, None)
plt.plot(naive_linear_predictions[last_hours], "x-", label="Raw time features")
plt.plot(
    cyclic_cossin_linear_predictions[last_hours],
    "x-",
    label="Trigonometric time features"
)
plt.plot(
    cyclic_spline_linear_preictions[last_hours],
    "x-",
    label="Spline-based time features"
)
plt.plot(
    y.iloc[test_0].values[last_hours],
    "x-",
    alpha=0.2,
    label="Actual demand",
    color="black",
)
plt.title("Predictions by linear models")
_ = plt.legend()

# %%
# We can draw the following conclusions from the above plot:
#
# - the raw ordinal time-related features are problematic because they do not
#   capture the natural periodicity: we observe a big jump in the predictions
#   at the end of each day when the hour features goes from 23 back to 0. We
#   can expect similar artifacts at the end of each week or each year.
#
# - as expected, the trigonometric features (sine and cosine) do not have this
#   discontinuties at midnight but the linear regression model still fails to
#   leverage those features to properly model intra-day variations. Using
#   trigonmetric features for higher harmonics or additional trigonometric
#   features for the natural period with different phases could potentially fix
#   this problem.
#
# - the periodic spline-based features fix those two problems at a time: they
#   give more expressivity to the linear model by making it possible to focus
#   on a specific hours thanks to the use of 12 knots. Furthermore the
#   `extrapolation="periodic"` option enforces a smooth representation between
#   `hour=23` and `hour=0`.
#
# Finally the lack of modeling of interactions between features seems to be a
# problem even for the model with spline-based features as its performance is
# qualitatively much worse than the same model with a polynomial kernel
# approximation on top as shown previously.
