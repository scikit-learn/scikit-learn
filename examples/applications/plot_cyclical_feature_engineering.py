
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
from sklearn.datasets import fetch_openml

X, y = fetch_openml(
    "Bike_Sharing_Demand", version=2, as_frame=True, return_X_y=True
)
# %%
y.max()

# %% [markdown]
#
# Let's rescale the target variable (absolute value of the demand) to predict a
# relative demand so that the mean absolute error is more easily interpreted as
# a fraction of the maximum demand and avoid numerical issues.

# %%
y /= 1000
y.hist(bins=30)

# %%
# The input feature data frame is a time annotated hourly log of variables
# describing the weather conditions. It includes both numerical and categorical
# variables. Note that the time information has already been expanded into
# several complementary columns.
#
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
# The season variable is well behaved:
#
X["season"].value_counts()

# %%
# Since the dataset is a time-ordered event log (hourly demand), we will use a
# time-sensitive cross-validation splitter to evaluate our demand forecasting
# model as realistically as possible:

from sklearn.model_selection import TimeSeriesSplit

ts_cv = TimeSeriesSplit(
    n_splits=5, gap=48, max_train_size=10000, test_size=1000,
)

all_splits = list(ts_cv.split(X, y))
train_0, test_0 = all_splits[0]

# %%
X.iloc[test_0]

# %%
X.iloc[train_0]

# %%
train_4, test_4 = all_splits[4]

# %%
X.iloc[test_4]

# %%
X.iloc[train_4]

# %%
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
# we only try the default hyper-parameeters for this model:
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
# is quite good for a model with default hyper-parameter that just required to
# make the categorical variables explicit. Note that the time related features
# are passed as such but this not much of a problem for tree-based models as
# they can learn a non-monotonic relationship between ordinal input features
# and the target.
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
# The performance is a not good: around 15% of the max demand on average. This
# is three times higher than the gradient boosting model. We can suspect that
# the non-monotonic original encoding of the periodic time-related features
# might prevent the linear regression model to preperly leverage the time
# information.
#
# Trigonometric Features
# ----------------------
#
# As a first attempt, we can try to encode each of those periodic features
# using a sine and cosine transform with the matching period.
from sklearn.preprocessing import FunctionTransformer


def sin_transformer(period):
    return FunctionTransformer(lambda x: np.sin(x / period * 2 * np.pi))


def cos_transformer(period):
    return FunctionTransformer(lambda x: np.sin(x / period * 2 * np.pi))


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
# Unfortunately this simple feature engineering does not seem to significantly
# improve the performance of our linear regression model.
#
# Periodic Spline Features
# ------------------------
#
# Instead we can try an alternative encoding of the periodic time-related
# features using spline transformations with a large-enough number of knots,
# and as a result a larger number of expanded features:
from sklearn.preprocessing import SplineTransformer

# one knot for each month in year
cyclic_month_splines = SplineTransformer(
    degree=3,
    n_knots=12,
    knots=np.linspace(0, 12, 12).reshape(12, 1),
    extrapolation="periodic",
)
# one knot for each day in week
cyclic_weekday_splines = SplineTransformer(
    degree=3,
    n_knots=7,
    knots=np.linspace(0, 7, 7).reshape(7, 1),
    extrapolation="periodic",
)
# one knot every 2 hours in the day
cyclic_hour_splines = SplineTransformer(
    degree=3,
    n_knots=24,
    knots=np.linspace(0, 24, 24).reshape(24, 1),
    extrapolation="periodic",
)

cyclic_spline_transformer = ColumnTransformer(
    transformers=[
        ("categorical", one_hot_encoder, categorical_columns),
        ("cyclic_month", cyclic_month_splines, ["month"]),
        ("cyclic_weekday", cyclic_weekday_splines, ["weekday"]),
        ("cyclic_hour", cyclic_hour_splines, ["hour"]),
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
# leverage the periodic time-related features. However such linear models
# cannot capture interactions between features (for instance the impact of rain
# might not be the same during the work days and the week-ends).
#
# To capture such interaction we could either use a partial polynomial
# expansion `PolynomicalFeatures(degree=2, interaction_only=True)` or use the
# Nystroem method to compute an approximate polynomial kernel expansion. Let
# try the latter for computational reasons:
from sklearn.preprocessing import MinMaxScaler
from sklearn.kernel_approximation import Nystroem


cyclic_spline_poly_pipeline = make_pipeline(
    cyclic_spline_transformer,
    MinMaxScaler(),
    Nystroem(kernel="poly", degree=2, n_components=300),
    RidgeCV(alphas=alphas),
)
evaluate(cyclic_spline_poly_pipeline, X, y, cv=ts_cv)

# %%
# We observe that this model performance can almost rival the performance of
# the gradient boosted trees.
#
# Note that while the final step of this pipeline is a linear regression model,
# the intermediate steps such as the spline feature extraction and the Nyström
# kernel approximation are highly non-linear. As the result the compound model
# is much more expressive than a simple linear model.
#
# Let have a qualitative look at the predictions of this model and the gradient
# boosted trees:

# %%
gbrt_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
gbrt_predictions = gbrt_pipeline.predict(X.iloc[test_0])
# %%
cyclic_spline_poly_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_spline_poly_predictions = cyclic_spline_poly_pipeline.predict(X.iloc[test_0])

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(6, 6))
plt.scatter(y.iloc[test_0].values, gbrt_predictions,
            alpha=0.2, label="GBDT")
plt.scatter(y.iloc[test_0].values, cyclic_spline_poly_predictions,
            alpha=0.2, label="Spline poly reg.")
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
# We can confirm this global analysis by zooming on the last 96 hours (4 days)
# of the test set:

plt.figure(figsize=(12, 4))
last_hours = slice(-96, None)
plt.plot(y.iloc[test_0].values[last_hours], "x-", label="True")
plt.plot(gbrt_predictions[last_hours], "x-", label="GBDT predictions")
plt.plot(cyclic_spline_poly_predictions[last_hours], "x-",
         label="Spline poly predictions")
_ = plt.legend()

# %%
