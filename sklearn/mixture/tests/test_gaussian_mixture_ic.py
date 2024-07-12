"""Testing for GaussianMixtureIC"""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_equal

from sklearn.exceptions import NotFittedError
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixtureIC
from sklearn.utils._param_validation import InvalidParameterError


def _test_wrong_inputs(X, error_type, **kws):
    with pytest.raises(error_type):
        gmIC = GaussianMixtureIC(**kws)
        gmIC.fit(X)


def _test_right_inputs(X, **kws):
    gmIC = GaussianMixtureIC(**kws)
    gmIC.fit(X)


def test_n_components():
    X = np.random.normal(0, 1, size=(100, 3))

    # min_components must be less than 1
    _test_wrong_inputs(X, ValueError, min_components=0)

    # min_components must be an integer
    _test_wrong_inputs(X, TypeError, min_components="1")

    # max_components must be at least min_components
    _test_wrong_inputs(X, ValueError, max_components=0)

    # max_components must be an integer
    _test_wrong_inputs(X, TypeError, max_components="1")

    # max_components must be at most n_samples
    _test_wrong_inputs(X, ValueError, max_components=101)

    # min_components must be at most n_samples
    _test_wrong_inputs(X, ValueError, **{"min_components": 101, "max_components": 102})


def test_input_param():
    X = np.random.normal(0, 1, size=(100, 3))

    # covariance type is not an array, string or list
    _test_wrong_inputs(X, InvalidParameterError, covariance_type=1)

    # covariance type is not in ['spherical', 'diag', 'tied', 'full', 'all']
    _test_wrong_inputs(X, InvalidParameterError, covariance_type="1")

    # several but not all covariance types in ['spherical', 'diag', 'tied', 'full']
    _test_right_inputs(X, covariance_type=["spherical", "diag"])

    # covariance type is 'all'
    _test_right_inputs(X, covariance_type="all")

    # criterion is not "aic" or "bic"
    _test_wrong_inputs(X, ValueError, criterion="cic")

    # n_init is not an integer
    _test_wrong_inputs(X, TypeError, n_init="1")

    # n_init must be at least 1
    _test_wrong_inputs(X, ValueError, n_init=0)


def test_predict_without_fit():
    X = np.random.normal(0, 1, size=(100, 3))

    with pytest.raises(NotFittedError):
        gmIC = GaussianMixtureIC(min_components=2)
        gmIC.predict(X)


def _test_two_class(**kws):
    """
    Easily separable two gaussian problem.
    """
    np.random.seed(1)

    n = 100
    d = 3

    X1 = np.random.normal(2, 0.5, size=(n, d))
    X2 = np.random.normal(-2, 0.5, size=(n, d))
    X = np.vstack((X1, X2))
    y = np.repeat([0, 1], n)

    # test BIC
    gmIC = GaussianMixtureIC(max_components=5, criterion="bic", **kws)
    gmIC.fit(X, y)
    n_components = gmIC.n_components_

    # Assert that the two cluster model is the best
    assert_equal(n_components, 2)

    # Assert that we get perfect clustering
    ari = adjusted_rand_score(y, gmIC.fit_predict(X))
    assert_allclose(ari, 1)

    # test AIC
    gmIC = GaussianMixtureIC(max_components=5, criterion="aic", **kws)
    gmIC.fit(X, y)
    n_components = gmIC.n_components_

    # AIC gets the number of components wrong
    assert_equal(n_components >= 1, True)
    assert_equal(n_components <= 5, True)


def test_two_class():
    _test_two_class()


def test_two_class_sequential_v_parallel():
    """
    Testing independence of results from the execution mode
    (sequential vs. parallel using ``joblib.Parallel``).
    """
    np.random.seed(1)

    n = 100
    d = 3

    X1 = np.random.normal(2, 0.75, size=(n, d))
    X2 = np.random.normal(-2, 0.5, size=(n, d))
    X = np.vstack((X1, X2))

    gmIC_parallel = GaussianMixtureIC(max_components=5, criterion="bic", n_jobs=-1)
    preds_parallel = gmIC_parallel.fit_predict(X)

    gmIC_sequential = GaussianMixtureIC(max_components=5, criterion="bic", n_jobs=1)
    preds_sequential = gmIC_sequential.fit_predict(X)

    # Results obtained with sequential and parallel executions
    # must be identical
    assert_equal(preds_parallel, preds_sequential)
