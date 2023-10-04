import numpy as np
import pytest

from sklearn.cluster import KMeans
from sklearn.datasets import make_classification, make_regression
from sklearn.inspection._pd_utils import (
    _check_feature_names,
    _get_feature_index,
    _robust_predict_for_scatter,
)
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.utils._testing import _convert_container


@pytest.mark.parametrize(
    "feature_names, array_type, expected_feature_names",
    [
        (None, "array", ["x0", "x1", "x2"]),
        (None, "dataframe", ["a", "b", "c"]),
        (np.array(["a", "b", "c"]), "array", ["a", "b", "c"]),
    ],
)
def test_check_feature_names(feature_names, array_type, expected_feature_names):
    X = np.random.randn(10, 3)
    column_names = ["a", "b", "c"]
    X = _convert_container(X, constructor_name=array_type, columns_name=column_names)
    feature_names_validated = _check_feature_names(X, feature_names)
    assert feature_names_validated == expected_feature_names


def test_check_feature_names_error():
    X = np.random.randn(10, 3)
    feature_names = ["a", "b", "c", "a"]
    msg = "feature_names should not contain duplicates."
    with pytest.raises(ValueError, match=msg):
        _check_feature_names(X, feature_names)


@pytest.mark.parametrize("fx, idx", [(0, 0), (1, 1), ("a", 0), ("b", 1), ("c", 2)])
def test_get_feature_index(fx, idx):
    feature_names = ["a", "b", "c"]
    assert _get_feature_index(fx, feature_names) == idx


@pytest.mark.parametrize(
    "fx, feature_names, err_msg",
    [
        ("a", None, "Cannot plot partial dependence for feature 'a'"),
        ("d", ["a", "b", "c"], "Feature 'd' not in feature_names"),
    ],
)
def test_get_feature_names_error(fx, feature_names, err_msg):
    with pytest.raises(ValueError, match=err_msg):
        _get_feature_index(fx, feature_names)


class TestRobustPredictForRegression:
    @pytest.fixture
    def est(self):
        X, y = make_regression(n_samples=100, n_features=5, random_state=42)
        return LinearRegression().fit(X, y)

    def test_predict(self, est):
        X = np.random.rand(100, 5)
        y_pred = _robust_predict_for_scatter(X, est, "auto")
        assert y_pred.shape == (100, 1)

    def test_not_fitted(self):
        X = np.random.rand(100, 5)
        est = LinearRegression()
        with pytest.raises(ValueError):
            _robust_predict_for_scatter(X, est, "auto")


class TestRobustPredictForClassification:
    @pytest.fixture
    def lr(self):
        X, y = make_classification(n_samples=100, n_features=5, random_state=42)
        return LogisticRegression().fit(X, y)

    @pytest.fixture
    def km(self):
        X, y = make_classification(n_samples=100, n_features=5, random_state=42)
        return KMeans(n_init=10).fit(X, y)

    def test_predict_proba(self, lr):
        X = np.random.rand(100, 5)
        y_pred = _robust_predict_for_scatter(X, lr, "predict_proba")
        assert y_pred.shape == (100, 1)

    def test_decision_function(self, lr):
        X = np.random.rand(100, 5)
        y_pred = _robust_predict_for_scatter(X, lr, "decision_function")
        assert y_pred.shape == (100, 1)

    def test_none_auto(self, km):
        with pytest.raises(ValueError):
            X = np.random.rand(100, 5)
            _ = _robust_predict_for_scatter(X, km, "auto")

    def test_none_predict_proba(self, km):
        with pytest.raises(ValueError):
            X = np.random.rand(100, 5)
            _ = _robust_predict_for_scatter(X, km, "predict_proba")

    def test_none_decision_function(self, km):
        with pytest.raises(ValueError):
            X = np.random.rand(100, 5)
            _ = _robust_predict_for_scatter(X, km, "decision_function")
