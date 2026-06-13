import math
import numbers
import time

import numpy as np
from sklearn.tree import DecisionTreeRegressor


def is_scalar_nan(x):
    return (
        not isinstance(x, numbers.Integral)
        and isinstance(x, numbers.Real)
        and math.isnan(x)
    )


def old_check(X):
    X = np.asarray(X, dtype=object)
    missing_mask = np.fromiter(
        (is_scalar_nan(x) for x in X.ravel()), dtype=bool, count=X.size
    )
    return np.any(missing_mask)


def new_encoded_check(X_encoded):
    return np.isnan(X_encoded).any()


def bench(name, func, repeats=3):
    best = float("inf")
    last = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        last = func()
        best = min(best, time.perf_counter() - t0)
    print(f"{name}: {best:.4f}s -> {last}")


def raises_missing(est, X):
    try:
        est._transform_categorical_features(X)
    except ValueError as exc:
        return str(exc)
    return "no error"


rng = np.random.RandomState(0)

X = rng.randn(100_000, 100)
X32 = X.astype(np.float32, copy=False)

print("Reported-size numeric check, shape=", X.shape)
bench("old scalar scan", lambda: old_check(X), repeats=1)
bench("new numeric isnan check", lambda: new_encoded_check(X32), repeats=5)

X_nan = X32.copy()
X_nan[-1, -1] = np.nan
bench("new numeric isnan check with one nan", lambda: new_encoded_check(X_nan), repeats=5)

print("\nTree categorical transform path")
X_cat = rng.randint(0, 8, size=(100_000, 20)).astype(np.float64)
est = DecisionTreeRegressor(
    categorical_features=np.arange(X_cat.shape[1]), random_state=0
)
est.is_categorical_ = np.ones(X_cat.shape[1], dtype=bool)
est._fit_categorical_features(X_cat)

bench(
    "new _transform_categorical_features clean",
    lambda: est._transform_categorical_features(X_cat).shape,
    repeats=5,
)

X_cat_nan = X_cat.copy()
X_cat_nan[-1, -1] = np.nan
bench(
    "new _transform_categorical_features with nan",
    lambda: raises_missing(est, X_cat_nan),
    repeats=5,
)