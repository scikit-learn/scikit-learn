import pytest
import scipy
import numpy as np
from numpy.testing import assert_array_equal

from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor


@pytest.mark.parametrize('n_features_to_select', (0, 5))
def test_bad_n_features_to_select(n_features_to_select):
    X, y = make_regression(n_features=5)
    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select)
    with pytest.raises(ValueError, match="must be in"):
        sfs.fit(X, y)


@pytest.mark.parametrize('forward', (True, False))
@pytest.mark.parametrize('n_features_to_select', (1, 5, 9, None))
def test_n_features_to_select(forward, n_features_to_select):
    # Make sure n_features_to_select is respected

    X, y = make_regression(n_features=10)
    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select,
                                    forward=forward, cv=2)
    sfs.fit(X, y)
    if n_features_to_select is None:
        n_features_to_select = 5  # n_features // 2
    assert sfs.get_support(indices=True).shape[0] == n_features_to_select
    assert sfs.n_features_to_select_ == n_features_to_select
    assert sfs.transform(X).shape[1] == n_features_to_select


@pytest.mark.parametrize('seed', range(10))
@pytest.mark.parametrize('forward', (True, False))
@pytest.mark.parametrize('n_features_to_select, expected_selected_features', [
    (2, [0, 2]),  # f1 is dropped since it has no predictive power
    (1, [2]),  # f2 is more predictive than f0 so it's kept
])
def test_sanity(seed, forward, n_features_to_select,
                expected_selected_features):
    # Basic sanity check: 3 features, only f0 and f2 are correlated with the
    # target, f2 having a stronger correlation than f0. We expect f1 to be
    # dropped, and f2 to always be selected.

    rng = np.random.RandomState(seed)
    n_samples = 100
    X = rng.randn(n_samples, 3)
    y = 3 * X[:, 0] - 10 * X[:, 2]

    sfs = SequentialFeatureSelector(LinearRegression(),
                                    n_features_to_select=n_features_to_select,
                                    forward=forward, cv=2)
    sfs.fit(X, y)
    assert_array_equal(sfs.get_support(indices=True),
                       expected_selected_features)


def test_sparse_support():
    # Make sure sparse data is supported

    X, y = make_regression(n_features=10)
    X = scipy.sparse.csr_matrix(X)
    sfs = SequentialFeatureSelector(LinearRegression(), cv=2)
    sfs.fit(X, y)
    sfs.transform(X)


def test_nan_support():
    # Make sure nans are OK if the underlying estimator supports nans

    rng = np.random.RandomState(0)
    n_samples, n_features = 100, 10
    X, y = make_regression(n_samples, n_features, random_state=0)
    nan_mask = rng.randint(0, 2, size=(n_samples, n_features), dtype=bool)
    X[nan_mask] = np.nan
    sfs = SequentialFeatureSelector(HistGradientBoostingRegressor(), cv=2)
    sfs.fit(X, y)
    sfs.transform(X)

    with pytest.raises(ValueError, match='Input contains NaN'):
        # LinearRegression does not support nans
        SequentialFeatureSelector(LinearRegression(), cv=2).fit(X, y)
