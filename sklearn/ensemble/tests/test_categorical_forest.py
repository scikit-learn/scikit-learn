"""
Testing for the forest module (sklearn.ensemble.forest).
"""

# Authors: Gilles Louppe,
#          Brian Holt,
#          Andreas Mueller,
#          Arnaud Joly
# License: BSD 3 clause

from typing import Any, Dict

import joblib
import numpy as np
import pytest

from sklearn import datasets
from sklearn.ensemble import (
    ExtraTreesClassifier,
    ExtraTreesRegressor,
    RandomForestClassifier,
    RandomForestRegressor,
    RandomTreesEmbedding,
)
from sklearn.utils.validation import check_random_state

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# Larger classification sample used for testing feature importances
X_large, y_large = datasets.make_classification(
    n_samples=500,
    n_features=10,
    n_informative=3,
    n_redundant=0,
    n_repeated=0,
    shuffle=False,
    random_state=0,
)

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = check_random_state(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# Make regression dataset
X_reg, y_reg = datasets.make_regression(n_samples=500, n_features=10, random_state=1)

# also make a hastie_10_2 dataset
hastie_X, hastie_y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
hastie_X = hastie_X.astype(np.float32)

# Get the default backend in joblib to test parallelism and interaction with
# different backends
DEFAULT_JOBLIB_BACKEND = joblib.parallel.get_active_backend()[0].__class__

FOREST_CLASSIFIERS = {
    "ExtraTreesClassifier": ExtraTreesClassifier,
    "RandomForestClassifier": RandomForestClassifier,
}

FOREST_REGRESSORS = {
    "ExtraTreesRegressor": ExtraTreesRegressor,
    "RandomForestRegressor": RandomForestRegressor,
}

FOREST_TRANSFORMERS = {
    "RandomTreesEmbedding": RandomTreesEmbedding,
}

FOREST_ESTIMATORS: Dict[str, Any] = dict()
FOREST_ESTIMATORS.update(FOREST_CLASSIFIERS)
FOREST_ESTIMATORS.update(FOREST_REGRESSORS)
FOREST_ESTIMATORS.update(FOREST_TRANSFORMERS)

FOREST_CLASSIFIERS_REGRESSORS: Dict[str, Any] = FOREST_CLASSIFIERS.copy()
FOREST_CLASSIFIERS_REGRESSORS.update(FOREST_REGRESSORS)


pytest.mark.parametrize("model", FOREST_CLASSIFIERS_REGRESSORS)


def _make_categorical(
    n_rows: int,
    n_numerical: int,
    n_categorical: int,
    cat_size: int,
    n_num_meaningful: int,
    n_cat_meaningful: int,
    regression: bool,
    return_tuple: bool,
    random_state: int,
):
    from sklearn.preprocessing import OneHotEncoder

    rng = np.random.RandomState(random_state)

    numeric = rng.standard_normal((n_rows, n_numerical))
    categorical = rng.randint(0, cat_size, (n_rows, n_categorical))
    categorical_ohe = OneHotEncoder(categories="auto").fit_transform(
        categorical[:, :n_cat_meaningful]
    )

    data_meaningful = np.hstack(
        (numeric[:, :n_num_meaningful], categorical_ohe.todense())
    )
    _, cols = data_meaningful.shape
    coefs = rng.standard_normal(cols)
    y = np.dot(data_meaningful, coefs)
    y = np.asarray(y).reshape(-1)
    X = np.hstack((numeric, categorical))

    if not regression:
        y = (y < y.mean()).astype(int)

    meaningful_features = np.r_[
        np.arange(n_num_meaningful), np.arange(n_cat_meaningful) + n_numerical
    ]

    if return_tuple:
        return X, y, meaningful_features
    else:
        return {"X": X, "y": y, "meaningful_features": meaningful_features}


@pytest.mark.parametrize("model", FOREST_CLASSIFIERS_REGRESSORS)
@pytest.mark.parametrize(
    "data_params",
    [
        {
            "n_rows": 10000,
            "n_numerical": 10,
            "n_categorical": 5,
            "cat_size": 3,
            "n_num_meaningful": 1,
            "n_cat_meaningful": 2,
        },
        {
            "n_rows": 1000,
            "n_numerical": 0,
            "n_categorical": 5,
            "cat_size": 3,
            "n_num_meaningful": 0,
            "n_cat_meaningful": 3,
        },
        {
            "n_rows": 1000,
            "n_numerical": 5,
            "n_categorical": 5,
            "cat_size": 64,
            "n_num_meaningful": 0,
            "n_cat_meaningful": 2,
        },
        {
            "n_rows": 1000,
            "n_numerical": 5,
            "n_categorical": 5,
            "cat_size": 3,
            "n_num_meaningful": 0,
            "n_cat_meaningful": 3,
        },
    ],
)
def test_categorical_data(model, data_params):
    # DecisionTrees are too slow for large category sizes.
    if data_params["cat_size"] > 8 and "RandomForest" in model:
        pass

    X, y, meaningful_features = _make_categorical(
        **data_params,
        regression=model in FOREST_REGRESSORS,
        return_tuple=True,
        random_state=42,
    )
    rows, cols = X.shape
    categorical_features = (
        np.arange(data_params["n_categorical"]) + data_params["n_numerical"]
    )

    model = FOREST_CLASSIFIERS_REGRESSORS[model](
        random_state=42, categorical=categorical_features, n_estimators=100
    ).fit(X, y)
    fi = model.feature_importances_
    bad_features = np.array([True] * cols)
    bad_features[meaningful_features] = False

    good_ones = fi[meaningful_features]
    print(good_ones)
    bad_ones = fi[bad_features]
    print(bad_ones)

    # all good features should be more important than all bad features.
    assert np.all([np.all(x > bad_ones) for x in good_ones])
