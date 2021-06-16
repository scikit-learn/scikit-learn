
# %%
from sklearn.datasets import fetch_openml

X, y = fetch_openml(
    "Bike_Sharing_Demand", version=2, as_frame=True, return_X_y=True
)
# %%
y.hist(bins=30)
y.max()

# %%
y /= 1000
y.hist(bins=30)

# %%
X
# %%
X["weather"].value_counts()

# %%
X["weather"][X["weather"] == "heavy_rain"] = "rain"
# %%
X["weather"].value_counts()

# %%
X["season"].value_counts()

# %%
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
%%time
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
%%time
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import RidgeCV
from sklearn.linear_model import SGDRegressor
import numpy as np


one_hot_encoder = OneHotEncoder(handle_unknown="ignore")
naive_linear_pipeline = make_pipeline(
    ColumnTransformer(
        transformers=[
            ("categorical", one_hot_encoder, categorical_columns),
        ],
        remainder="passthrough",
    ),
    RidgeCV(alphas=np.logspace(-6, 6, 24)),
)


evaluate(naive_linear_pipeline, X, y, cv=ts_cv)


# %%
%%time
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
    RidgeCV(alphas=np.logspace(-6, 6, 24)),
)
evaluate(cyclic_spline_linear_pipeline, X, y, cv=ts_cv)


# %%
%%time
from sklearn.preprocessing import MinMaxScaler
from sklearn.kernel_approximation import Nystroem


cyclic_spline_poly_pipeline = make_pipeline(
    cyclic_spline_transformer,
    MinMaxScaler(),
    Nystroem(kernel="poly", degree=2, n_components=300),
    RidgeCV(alphas=np.logspace(-6, 6, 24)),
)
evaluate(cyclic_spline_poly_pipeline, X, y, cv=ts_cv)

# %%
gbrt_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
gbrt_predictions = gbrt_pipeline.predict(X.iloc[test_0])
# %%
cyclic_spline_poly_pipeline.fit(X.iloc[train_0], y.iloc[train_0])
cyclic_spline_poly_predictions = cyclic_spline_poly_pipeline.predict(X.iloc[test_0])

# %%
plt.scatter(y.iloc[test_0].values, gbrt_predictions, alpha=0.1)
plt.scatter(y.iloc[test_0].values, cyclic_spline_poly_predictions, alpha=0.1)
plt.plot([0, 1], [0, 1], "--")
plt.xlim(0, 1)
plt.ylim(0, 1)

# %%
import matplotlib.pyplot as plt

plt.figure(figsize=(12, 4))
last_hours = slice(-96, None)
plt.plot(y.iloc[test_0].values[last_hours], label="True")
plt.plot(gbrt_predictions[last_hours], label="GBDT predictions")
plt.plot(cyclic_spline_poly_predictions[last_hours], label="Spline poly predictions")
plt.legend();
# %%
