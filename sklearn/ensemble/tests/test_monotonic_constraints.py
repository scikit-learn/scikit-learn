import re

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from sklearn.datasets import make_classification, make_regression
from sklearn.ensemble import GradientBoostingClassifier, GradientBoostingRegressor
from sklearn.utils._testing import _convert_container
from sklearn.utils.fixes import CSC_CONTAINERS


@pytest.mark.parametrize("depth_first_builder", (True, False))
@pytest.mark.parametrize("criterion", ("friedman_mse", "squared_error"))
@pytest.mark.parametrize("csc_container", [None] + CSC_CONTAINERS)
def test_monotonic_constraints_regressor(
    depth_first_builder, criterion, global_random_seed, csc_container
):
    n_samples = 1000
    n_samples_train = 900
    # Build a regression task using 5 informative features
    X, y = make_regression(
        n_samples=n_samples,
        n_features=5,
        n_informative=5,
        random_state=global_random_seed,
    )
    train = np.arange(n_samples_train)
    test = np.arange(n_samples_train, n_samples)
    X_train = X[train]
    y_train = y[train]
    X_test = np.copy(X[test])
    X_test_incr = np.copy(X_test)
    X_test_decr = np.copy(X_test)
    X_test_incr[:, 0] += 10
    X_test_decr[:, 1] += 10
    monotonic_cst = np.zeros(X.shape[1])
    monotonic_cst[0] = 1
    monotonic_cst[1] = -1

    if csc_container is not None:
        X_train = csc_container(X_train)

    params = {
        "monotonic_cst": monotonic_cst,
        "criterion": criterion,
        "random_state": global_random_seed,
        "n_estimators": 5,
    }
    if depth_first_builder:
        params["max_depth"] = None
    else:
        params["max_depth"] = 8
        params["max_leaf_nodes"] = n_samples_train

    est = GradientBoostingRegressor(**params)
    est.fit(X_train, y_train)
    y = est.predict(X_test)

    # Monotonic increase constraint
    y_incr = est.predict(X_test_incr)
    # y_incr should always be greater than y
    assert np.all(y_incr >= y)

    # Monotonic decrease constraint
    y_decr = est.predict(X_test_decr)
    # y_decr should always be lower than y
    assert np.all(y_decr <= y)


@pytest.mark.parametrize("depth_first_builder", (True, False))
@pytest.mark.parametrize("csc_container", [None] + CSC_CONTAINERS)
def test_monotonic_constraints_classifier(
    depth_first_builder, global_random_seed, csc_container
):
    n_samples = 1000
    n_samples_train = 900
    X, y = make_classification(
        n_samples=n_samples,
        n_classes=2,
        n_features=5,
        n_informative=5,
        n_redundant=0,
        random_state=global_random_seed,
    )
    X_train, y_train = X[:n_samples_train], y[:n_samples_train]
    X_test, _ = X[n_samples_train:], y[n_samples_train:]

    X_test_0incr, X_test_0decr = np.copy(X_test), np.copy(X_test)
    X_test_1incr, X_test_1decr = np.copy(X_test), np.copy(X_test)
    X_test_0incr[:, 0] += 10
    X_test_0decr[:, 0] -= 10
    X_test_1incr[:, 1] += 10
    X_test_1decr[:, 1] -= 10
    monotonic_cst = np.zeros(X.shape[1])
    monotonic_cst[0] = 1
    monotonic_cst[1] = -1

    if csc_container is not None:
        X_train = csc_container(X_train)

    params = {
        "monotonic_cst": monotonic_cst,
        "random_state": global_random_seed,
        "n_estimators": 5,
    }
    if depth_first_builder:
        params["max_depth"] = None
    else:
        params["max_depth"] = 8
        params["max_leaf_nodes"] = n_samples_train

    est = GradientBoostingClassifier(**params)
    est.fit(X_train, y_train)
    proba_test = est.predict_proba(X_test)

    assert np.logical_and(
        proba_test >= 0.0, proba_test <= 1.0
    ).all(), "Probability should always be in [0, 1] range."
    assert_array_almost_equal(proba_test.sum(axis=1), 1.0)

    # Monotonic increase constraint, it applies to the positive class
    assert np.all(est.predict_proba(X_test_0incr)[:, 1] >= proba_test[:, 1])
    assert np.all(est.predict_proba(X_test_0decr)[:, 1] <= proba_test[:, 1])

    # Monotonic decrease constraint, it applies to the positive class
    assert np.all(est.predict_proba(X_test_1incr)[:, 1] <= proba_test[:, 1])
    assert np.all(est.predict_proba(X_test_1decr)[:, 1] >= proba_test[:, 1])


@pytest.mark.parametrize("use_feature_names", (True, False))
def test_predictions(global_random_seed, use_feature_names):
    # This is the same test as the one for HistGradientBoostingRegressor
    rng = np.random.RandomState(global_random_seed)

    n_samples = 1000
    f_0 = rng.rand(n_samples)  # positive correlation with y
    f_1 = rng.rand(n_samples)  # negative correslation with y
    X = np.c_[f_0, f_1]
    columns_name = ["f_0", "f_1"]
    constructor_name = "dataframe" if use_feature_names else "array"
    X = _convert_container(X, constructor_name, columns_name=columns_name)

    noise = rng.normal(loc=0.0, scale=0.01, size=n_samples)
    y = 5 * f_0 + np.sin(10 * np.pi * f_0) - 5 * f_1 - np.cos(10 * np.pi * f_1) + noise

    if use_feature_names:
        monotonic_cst = {"f_0": +1, "f_1": -1}
    else:
        monotonic_cst = [+1, -1]

    gbdt = GradientBoostingRegressor(monotonic_cst=monotonic_cst)
    gbdt.fit(X, y)

    linspace = np.linspace(0, 1, 100)
    sin = np.sin(linspace)
    constant = np.full_like(linspace, fill_value=0.5)

    # First feature (POS)
    # assert pred is all increasing when f_0 is all increasing
    X = np.c_[linspace, constant]
    X = _convert_container(X, constructor_name, columns_name=columns_name)
    pred = gbdt.predict(X)
    assert (np.diff(pred) >= 0.0).all()
    # assert pred actually follows the variations of f_0
    X = np.c_[sin, constant]
    X = _convert_container(X, constructor_name, columns_name=columns_name)
    pred = gbdt.predict(X)
    assert np.all((np.diff(pred) >= 0) == (np.diff(sin) >= 0))

    # Second feature (NEG)
    # assert pred is all decreasing when f_1 is all increasing
    X = np.c_[constant, linspace]
    X = _convert_container(X, constructor_name, columns_name=columns_name)
    pred = gbdt.predict(X)
    assert (np.diff(pred) <= 0.0).all()
    # assert pred actually follows the inverse variations of f_1
    X = np.c_[constant, sin]
    X = _convert_container(X, constructor_name, columns_name=columns_name)
    pred = gbdt.predict(X)
    assert ((np.diff(pred) <= 0) == (np.diff(sin) >= 0)).all()


def test_multiclass_raises():
    X, y = make_classification(
        n_samples=100, n_features=5, n_classes=3, n_informative=3, random_state=0
    )
    y[0] = 0
    monotonic_cst = np.zeros(X.shape[1])
    monotonic_cst[0] = -1
    monotonic_cst[1] = 1
    est = GradientBoostingClassifier(
        max_depth=None, monotonic_cst=monotonic_cst, random_state=0
    )

    msg = "Monotonicity constraints are not supported with multiclass classification"
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)


def test_multiple_output_raises():
    X = [[1, 2, 3, 4, 5], [6, 7, 8, 9, 10]]
    y = [[1, 0, 1, 0, 1], [1, 0, 1, 0, 1]]

    est = GradientBoostingClassifier(
        max_depth=None, monotonic_cst=np.array([-1, 1]), random_state=0
    )
    msg = "Monotonicity constraints are not supported with multiple output"
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)


def test_missing_values_raises():
    X, y = make_classification(
        n_samples=100, n_features=5, n_classes=2, n_informative=3, random_state=0
    )
    X[0, 0] = np.nan
    monotonic_cst = np.zeros(X.shape[1])
    monotonic_cst[0] = 1
    est = GradientBoostingClassifier(
        max_depth=None, monotonic_cst=monotonic_cst, random_state=0
    )

    msg = "Input X contains NaN"
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)


def test_bad_monotonic_cst_raises():
    X = [[1, 2], [3, 4], [5, 6], [7, 8], [9, 10]]
    y = [1, 0, 1, 0, 1]

    msg = "monotonic_cst has shape (3,) but the input data X has 2 features."
    est = GradientBoostingClassifier(
        max_depth=None, monotonic_cst=np.array([-1, 1, 0]), random_state=0
    )
    with pytest.raises(ValueError, match=re.escape(msg)):
        est.fit(X, y)

    msg = "monotonic_cst must be an array-like of -1, 0 or 1."
    est = GradientBoostingClassifier(
        max_depth=None, monotonic_cst=np.array([-2, 2]), random_state=0
    )
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)

    est = GradientBoostingClassifier(
        max_depth=None, monotonic_cst=np.array([-1, 0.8]), random_state=0
    )
    with pytest.raises(ValueError, match=msg + "(.*)0.8]"):
        est.fit(X, y)


def test_bad_monotonic_cst_related_to_feature_names():
    pd = pytest.importorskip("pandas")
    X = pd.DataFrame({"a": [0, 1, 2], "b": [0, 1, 2]})
    y = np.array([0, 1, 0])

    monotonic_cst = {"d": 1, "a": 1, "c": -1}
    gbdt = GradientBoostingRegressor(monotonic_cst=monotonic_cst)
    expected_msg = re.escape(
        "monotonic_cst contains 2 unexpected feature names: ['c', 'd']."
    )
    with pytest.raises(ValueError, match=expected_msg):
        gbdt.fit(X, y)

    monotonic_cst = {k: 1 for k in "abcdefghijklmnopqrstuvwxyz"}
    gbdt = GradientBoostingRegressor(monotonic_cst=monotonic_cst)
    expected_msg = re.escape(
        "monotonic_cst contains 24 unexpected feature names: "
        "['c', 'd', 'e', 'f', 'g', '...']."
    )
    with pytest.raises(ValueError, match=expected_msg):
        gbdt.fit(X, y)

    monotonic_cst = {"a": 1}
    gbdt = GradientBoostingRegressor(monotonic_cst=monotonic_cst)
    expected_msg = re.escape(
        "GradientBoostingRegressor was not fitted on data with feature "
        "names. Pass monotonic_cst as an integer array instead."
    )
    with pytest.raises(ValueError, match=expected_msg):
        gbdt.fit(X.values, y)

    monotonic_cst = {"b": -1, "a": "+"}
    gbdt = GradientBoostingRegressor(monotonic_cst=monotonic_cst)
    expected_msg = re.escape("monotonic_cst['a'] must be either -1, 0 or 1. Got '+'.")
    with pytest.raises(ValueError, match=expected_msg):
        gbdt.fit(X, y)
