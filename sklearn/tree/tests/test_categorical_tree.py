import numpy as np
import pytest
from scipy.sparse import csc_matrix

from sklearn import datasets
from sklearn.random_projection import _sparse_random_matrix
from sklearn.tree import (
    DecisionTreeClassifier,
    DecisionTreeRegressor,
    ExtraTreeClassifier,
    ExtraTreeRegressor,
)
from sklearn.utils.validation import check_random_state

CLF_CRITERIONS = ("gini", "log_loss")
REG_CRITERIONS = ("squared_error", "absolute_error", "friedman_mse", "poisson")

CLF_TREES = {
    "DecisionTreeClassifier": DecisionTreeClassifier,
    "ExtraTreeClassifier": ExtraTreeClassifier,
}

REG_TREES = {
    "DecisionTreeRegressor": DecisionTreeRegressor,
    "ExtraTreeRegressor": ExtraTreeRegressor,
}

ALL_TREES: dict = dict()
ALL_TREES.update(CLF_TREES)
ALL_TREES.update(REG_TREES)


X_small = np.array(
    [
        [0, 0, 4, 0, 0, 0, 1, -14, 0, -4, 0, 0, 0, 0],
        [0, 0, 5, 3, 0, -4, 0, 0, 1, -5, 0.2, 0, 4, 1],
        [-1, -1, 0, 0, -4.5, 0, 0, 2.1, 1, 0, 0, -4.5, 0, 1],
        [-1, -1, 0, -1.2, 0, 0, 0, 0, 0, 0, 0.2, 0, 0, 1],
        [-1, -1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1],
        [-1, -2, 0, 4, -3, 10, 4, 0, -3.2, 0, 4, 3, -4, 1],
        [2.11, 0, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0.5, 0, -3, 1],
        [2.11, 0, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0, 0, -2, 1],
        [2.11, 8, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0, 0, -2, 1],
        [2.11, 8, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0.5, 0, -1, 0],
        [2, 8, 5, 1, 0.5, -4, 10, 0, 1, -5, 3, 0, 2, 0],
        [2, 0, 1, 1, 1, -1, 1, 0, 0, -2, 3, 0, 1, 0],
        [2, 0, 1, 2, 3, -1, 10, 2, 0, -1, 1, 2, 2, 0],
        [1, 1, 0, 2, 2, -1, 1, 2, 0, -5, 1, 2, 3, 0],
        [3, 1, 0, 3, 0, -4, 10, 0, 1, -5, 3, 0, 3, 1],
        [2.11, 8, -6, -0.5, 0, 1, 0, 0, -3.2, 6, 0.5, 0, -3, 1],
        [2.11, 8, -6, -0.5, 0, 1, 0, 0, -3.2, 6, 1.5, 1, -1, -1],
        [2.11, 8, -6, -0.5, 0, 10, 0, 0, -3.2, 6, 0.5, 0, -1, -1],
        [2, 0, 5, 1, 0.5, -2, 10, 0, 1, -5, 3, 1, 0, -1],
        [2, 0, 1, 1, 1, -2, 1, 0, 0, -2, 0, 0, 0, 1],
        [2, 1, 1, 1, 2, -1, 10, 2, 0, -1, 0, 2, 1, 1],
        [1, 1, 0, 0, 1, -3, 1, 2, 0, -5, 1, 2, 1, 1],
        [3, 1, 0, 1, 0, -4, 1, 0, 1, -2, 0, 0, 1, 0],
    ]
)

y_small = [1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0]
y_small_reg = [
    1.0,
    2.1,
    1.2,
    0.05,
    10,
    2.4,
    3.1,
    1.01,
    0.01,
    2.98,
    3.1,
    1.1,
    0.0,
    1.2,
    2,
    11,
    0,
    0,
    4.5,
    0.201,
    1.06,
    0.9,
    0,
]

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(1)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the diabetes dataset
# and randomly permute it
diabetes = datasets.load_diabetes()
perm = rng.permutation(diabetes.target.size)
diabetes.data = diabetes.data[perm]
diabetes.target = diabetes.target[perm]

digits = datasets.load_digits()
perm = rng.permutation(digits.target.size)
digits.data = digits.data[perm]
digits.target = digits.target[perm]

random_state = check_random_state(0)
X_multilabel, y_multilabel = datasets.make_multilabel_classification(
    random_state=0, n_samples=30, n_features=10
)

# NB: despite their names X_sparse_* are numpy arrays (and not sparse matrices)
X_sparse_pos = random_state.uniform(size=(20, 5))
X_sparse_pos[X_sparse_pos <= 0.8] = 0.0
y_random = random_state.randint(0, 4, size=(20,))
X_sparse_mix = _sparse_random_matrix(20, 10, density=0.25, random_state=0).toarray()


DATASETS = {
    "iris": {"X": iris.data, "y": iris.target},
    "diabetes": {"X": diabetes.data, "y": diabetes.target},
    "digits": {"X": digits.data, "y": digits.target},
    "toy": {"X": X, "y": y},
    "clf_small": {"X": X_small, "y": y_small},
    "reg_small": {"X": X_small, "y": y_small_reg},
    "multilabel": {"X": X_multilabel, "y": y_multilabel},
    "sparse-pos": {"X": X_sparse_pos, "y": y_random},
    "sparse-neg": {"X": -X_sparse_pos, "y": y_random},
    "sparse-mix": {"X": X_sparse_mix, "y": y_random},
    "zeros": {"X": np.zeros((20, 3)), "y": y_random},
}
for name in DATASETS:
    DATASETS[name]["X_sparse"] = csc_matrix(DATASETS[name]["X"])


@pytest.mark.parametrize("name", ALL_TREES)
@pytest.mark.parametrize(
    "categorical",
    ["invalid string", [[0]], [False, False, False], [1, 2], [-3], [0, 0, 1]],
)
def test_invalid_categorical(name, categorical):
    Tree = ALL_TREES[name]
    if categorical == "invalid string":
        with pytest.raises(ValueError, match="The 'categorical' parameter"):
            Tree(categorical=categorical).fit(X, y)
    else:
        with pytest.raises(ValueError, match="Invalid value for categorical"):
            Tree(categorical=categorical).fit(X, y)


@pytest.mark.parametrize("name", ALL_TREES)
def test_no_sparse_with_categorical(name):
    # Currently we do not support sparse categorical features
    X, y, X_sparse = [DATASETS["clf_small"][z] for z in ["X", "y", "X_sparse"]]
    Tree = ALL_TREES[name]
    with pytest.raises(
        NotImplementedError, match="Categorical features not supported with sparse"
    ):
        Tree(categorical=[6, 10]).fit(X_sparse, y)

    with pytest.raises(
        NotImplementedError, match="Categorical features not supported with sparse"
    ):
        Tree(categorical=[6, 10]).fit(X, y).predict(X_sparse)


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


@pytest.mark.parametrize("model", ALL_TREES)
@pytest.mark.parametrize(
    "data_params",
    [
        {
            "n_rows": 1000,
            "n_numerical": 5,
            "n_categorical": 5,
            "cat_size": 3,
            "n_num_meaningful": 2,
            "n_cat_meaningful": 3,
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
    if data_params["cat_size"] > 8 and "DecisionTree" in model:
        pass

    X, y, meaningful_features = _make_categorical(
        **data_params, regression=model in REG_TREES, return_tuple=True, random_state=42
    )
    rows, cols = X.shape
    categorical_features = (
        np.arange(data_params["n_categorical"]) + data_params["n_numerical"]
    )

    model = ALL_TREES[model](random_state=43, categorical=categorical_features).fit(
        X, y
    )
    fi = model.feature_importances_
    bad_features = np.array([True] * cols)
    bad_features[meaningful_features] = False

    good_ones = fi[meaningful_features]
    bad_ones = fi[bad_features]

    # all good features should be more important than all bad features.
    # XXX: or at least a large fraction of them
    # assert np.all([np.all(x > bad_ones) for x in good_ones])
    assert np.mean([np.all(x > bad_ones) for x in good_ones]) > 0.7

    leaves = model.tree_.children_left < 0
    assert np.all(model.tree_.impurity[leaves] < 1e-6)
