
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
# Data exploration on the Bike Sharing Demand dataset
# ---------------------------------------------------
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
X['weather'].replace(to_replace='heavy_rain', value="rain", inplace=True)
# %%
X["weather"].value_counts()

# %%
# As expected season variable is well balanced:
#
X["season"].value_counts()

# %%
# Time-based cross-validation
# ---------------------------
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
# Naive linear regression
# -----------------------
#
# As usual for linear models, categorical variables need to be one-hot encoded.
# As the numerical features are all approximately on similar scales, we let
# them go through unscaled, including the time-related features:
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import RidgeCV
import numpy as np


one_hot_encoder = OneHotEncoder(handle_unknown="ignore", sparse=False)
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
# The performance is a not good: the everage error is around 14% of the maximum
# demand. This is more than three times higher than the average error of the
# gradient boosting model. We can suspect that the naive original encoding of
# the periodic time-related features might prevent the linear regression model
# to properly leverage the time information: linear regression cannot model a
# non-monotonic relationship between the input features and the target.
#
# For example, the raw numerical encoding of the "hour" feature prevents the
# linear model to model that an increase of hour in the morning from 6 to 8
# should have a strong positive impact on the number of bike rentals while a
# increase of similar magnitude in the evening from 18 to 20 should have a
# strong negative impact on the predicted number of bike rentals.
#
# Time-steps as categories
# ------------------------
#
# Since the time features are encoded in a discrete manner using integers (24
# unique values in the "hours" feature), we could decide to treat those as
# categorical variables and ignore any assumption implied by the ordering of
# the hour values using a one-hot encoding.
#
# Using one-hot encoding for the time features gives the linear model a lot
# more flexibility as we introduce one additional feature per discrete time
# level.

one_hot_linear_pipeline = make_pipeline(
    ColumnTransformer(
        transformers=[
            ("categorical", one_hot_encoder, categorical_columns),
            ("one_hot_time", one_hot_encoder, ["hour", "weekday", "month"]),
        ],
        remainder="passthrough",
    ),
    RidgeCV(alphas=alphas),
)

evaluate(one_hot_linear_pipeline, X, y, cv=ts_cv)

# The average error rate of this model is 10% which is much better than using
# the original ordinal encoding of the time feature, confirming our intuition
# that the linear regression model benefit from the added flexibility to not
# treat time progression in a monotonic manner.
#
# However, this introduces a very large number of new features. If the time of
# the day was represented in minutes since the start of the day instead of
# hours, one-hot encoding would have introduced 1440 features instead of 24.
# This could cause some significant overfitting. To avoid this we could use
# :func:`sklearn.preprocessing.KBinsDiscretizer` instead to re-bin the number
# of levels of fine-grained ordinal or numerical variables while still
# benefitting from the non-monotonic expressivity advantages of one-hot
# encoding.
#
# Finally, we also observe than one-hot encoding completely ignores the
# ordering of the hour levels while this could be an interesting inductive bias
# to preserve to some level. In the following we try to explore smooth, non
# monotonic encoding that locally preserve the relative ordering of time
# features.
#
# Trigonometric features
# ----------------------
#
# As a first attempt, we can try to encode each of those periodic features
# using a sine and cosine transform with the matching period.
#
# Each ordinal time feature is transformed into 2 features that together encode
# equivalent information in a non-monotonic way, and more importantly without
# any jump between the first and the last value of the periodic range.
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
# The performance of our linear regression model with this simple feature
# engineering is a bit better than using the original ordinal time features but
# worse than using the one-hot encoded time features. We will further analyze
# possible reasons for this disappointing outcome at the end of this notebook.
#
# Periodic spline features
# ------------------------
#
# We can try an alternative encoding of the periodic time-related features
# using spline transformations with a large-enough number of knots, and as a
# result a larger number of expanded features:
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
# Again, let us visualize the effect of this feature expansion on some
# synthetic hour data with a bit of extrapolation beyond hour=23:
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
# Thanks for the use of the `extrapolation="periodic"` parameter, we observe
# that the feature encoding stays smooth when extrapolating beyond midnight.
#
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
# ~10% of the maximum demand, which is similar to what we observed with the
# one-hot encoded features.
#
# Qualitative analysis of the impact of features on linear models predictions
# ---------------------------------------------------------------------------
#
# Here we want to visualize the impact of the feature engineering choices on
# the time related-shape of the predictions.
#
# To do so we consider an arbitrary time-based split to compare the predictions
# on a range of held out data points.
naive_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
naive_linear_predictions = naive_linear_pipeline.predict(X.iloc[test_0])

one_hot_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
one_hot_linear_predictions = one_hot_linear_pipeline.predict(X.iloc[test_0])

cyclic_cossin_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_cossin_linear_predictions = cyclic_cossin_linear_pipeline.predict(
    X.iloc[test_0]
)

cyclic_spline_linear_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_spline_linear_preictions = cyclic_spline_linear_pipeline.predict(
    X.iloc[test_0]
)

# %%
# We visualize those predictions by zooming on the last 96 hours (4 days) of
# the test set to get some qualitative insights:
plt.figure(figsize=(12, 4))
last_hours = slice(-96, None)
plt.plot(naive_linear_predictions[last_hours], "x-", label="Ordinal time features")
plt.plot(
    one_hot_linear_predictions[last_hours],
    "x-",
    label="One-hot time features"
)
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
#   discontinuities at midnight but the linear regression model still fails to
#   leverage those features to properly model intra-day variations. Using
#   trigonmetric features for higher harmonics or additional trigonometric
#   features for the natural period with different phases could potentially fix
#   this problem.
#
# - the periodic spline-based features fix those two problems at once: they
#   give more expressivity to the linear model by making it possible to focus
#   on a specific hours thanks to the use of 12 knots. Furthermore the
#   `extrapolation="periodic"` option enforces a smooth representation between
#   `hour=23` and `hour=0`.
#
# - the one-hot encoded features behave similarly to the periodic spline-based
#   features but are more spiky: for instance they can better model the morning
#   peak during the weak days since this peak last fewer than an hour. But
#   remember that what can be an advantage for linear models is not necessarily
#   one for kernel models.
#
# Finally we observe that none of the linear models can approximate the true
# bike rentals demand, especially for the peaks that can be very sharp at rush
# hours during the working days but much flatter during the week-ends: the most
# accurate linear models based on splines or one-hot encoding tend to forecast
# peaks of commuting-related bike rentals even on the week-ends and
# under-estimate the commuting-related events during the working days.
#
# These systematic prediction errors reveal a form of under-fitting and can be
# explained by the lack of non-additive modeling of the interactions between
# features (in this case "workingday" and features derived from "hours"). This
# issue will be addressed in the following section.
#
# Modeling non-linear feature interactions with kernels
# -----------------------------------------------------
#
# The previous analysis highlighted the need to model the interactions between
# "workingday" and "hours". Another example of a such a non-linear interactions
# that we would like to model could be the impact of the rain that might not be
# the same during the working days and the week-ends and holidays for instance.
#
# Linear models cannot capture non-linear interactions between features even if
# the features them-selves are marginally non-linear as is the case with
# features constructed by `SplineTransformer` (or one-hot encoding or binning).
#
# To capture such interactions we could either use a partial polynomial
# expansion `PolynomicalFeatures(degree=2, interaction_only=True)` or use the
# Nyström method to compute an approximate polynomial kernel expansion.
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
#
# We observe that this model can almost rival the performance of the gradient
# boosted trees with an average error around 6% of the maximum demand.
#
# Note that while the final step of this pipeline is a linear regression model,
# the intermediate steps such as the spline feature extraction and the Nyström
# kernel approximation are highly non-linear. As a result the compound pipeline
# is much more expressive than a simple linear regression model.
#
# For the sake of completeness, we also evaluate the combination of one-hot
# encoding and kernel approximation:

one_hot_poly_pipeline = make_pipeline(
    ColumnTransformer(
        transformers=[
            ("categorical", one_hot_encoder, categorical_columns),
            ("one_hot_time", one_hot_encoder, ["hour", "weekday", "month"]),
        ],
        remainder="passthrough",
    ),
    Nystroem(kernel="poly", degree=2, n_components=300, random_state=0),
    RidgeCV(alphas=alphas),
)
evaluate(one_hot_poly_pipeline, X, y, cv=ts_cv)


# %%
# While one-hot features were competitive with spline-based features when using
# linear models, this is no longer the case when using a low-rank approximation
# of a non-linear kernel: this can be explained by the fact that spline
# features are smoother and allow the kernel approximation to find a more
# expressive decision function.
#
# Let us now have a qualitative look at the predictions of the kernel models
# and of the gradient boosted trees that should be able to better model
# non-linear interactions between features:
gbrt_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
gbrt_predictions = gbrt_pipeline.predict(X.iloc[test_0])

one_hot_poly_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
one_hot_poly_predictions = one_hot_poly_pipeline.predict(X.iloc[test_0])

cyclic_spline_poly_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_spline_poly_predictions = cyclic_spline_poly_pipeline.predict(X.iloc[test_0])

# %%
# Again we zoom on the last 4 days of the test set:

last_hours = slice(-96, None)
plt.figure(figsize=(12, 4))
plt.plot(
    gbrt_predictions[last_hours],
    "x-",
    label="Gradient Boosted Trees",
)
plt.plot(
    one_hot_poly_predictions[last_hours],
    "x-",
    label="One-hot + polynomial kernel",
)
plt.plot(
    cyclic_spline_poly_predictions[last_hours],
    "x-",
    label="Splines + polynomial kernel",
)
plt.title("Predictions by non-linear regression models")
plt.plot(
    y.iloc[test_0].values[last_hours],
    "x-",
    alpha=0.2,
    label="Actual demand",
    color="black",
)
_ = plt.legend()


# %%
# First note that trees can naturally model non-linear feature interactions
# since by default decision trees are allowed to grow beyond a depth of 2
# levels.
#
# Here we can observe that the combinations of spline features and non-linear
# kernels works quite well and can almost rival the accuracy of the gradient
# boosting regression trees.
#
# On the contrary, one-hot time features do not perform that well with the low
# rank kernel model. In particular they significantly over-estimate the low
# demand hours more than the competing models.
#
# We also observe than none of the models can successully predict some of the
# peak rentals at the rush hours during the working days. It is possible that
# access to additional features would be required to further improve the
# accuracy of the predictions. For instance, it could be useful to have access
# to the geographical repartition of the fleet at any point in time or the
# fraction of bikes that are immobilized because they need servicing.
#
# Let us finally get a more quantative look at the prediction errors of those
# three models using the true vs predicted demand scatter plots:
fig, axes = plt.subplots(ncols=3, figsize=(12, 4), sharey=True)
plt.suptitle("Non-linear regression models")
predictions = [
    one_hot_poly_predictions,
    cyclic_spline_poly_predictions,
    gbrt_predictions,
]
labels = [
    "One hot + polynomial kernel",
    "Splines + polynomial kernel",
    "Gradient Boosted Trees"
]
for ax, pred, label in zip(axes, predictions, labels):
    ax.scatter(y.iloc[test_0].values, pred, alpha=0.3, label=label)
    ax.plot([0, 1], [0, 1], "--", label="Perfect model")
    ax.set(
        xlim=(0, 1),
        ylim=(0, 1),
        xlabel="True demand",
        ylabel="Predicted demand",
    )
    ax.legend()


# %%
# This visualization confirms the conclusions we draw on the previous plot.
#
# All models under-estimate the high demand events (working days rush hours),
# but gradient boosting a bit less so. The low demand events are well predicted
# on average by gradient boosting while the one-hot polynomial regression
# pipeline seems to systematically over-estimate demand in that regime. Overall
# the predictions of the gradient boosted trees are closer to the diagonal than
# for the kernel models.
#
# Finally, we note that we could have obtained slightly better results for
# kernel models by using more components (higher rank kernel approximation) at
# the cost of longer fit and prediction times. For `n_components`, the
# performance of the one-hot features would even match the spline features.
#
# The `Nystroem` + `RidgeCV` classifier could also have been replaced by
# `sklearn.neural_network.MLPRegressor` with one or two hidden layers and we
# would have obtained quite similar results.
