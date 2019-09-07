# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD 3 clause

import pytest

import numpy as np
from scipy import sparse
from scipy import linalg

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_allclose

from sklearn.linear_model.base import LinearRegression
from sklearn.linear_model.base import _preprocess_data
from sklearn.linear_model.base import _rescale_data
from sklearn.linear_model.base import make_dataset
from sklearn.utils import check_random_state
from sklearn.datasets.samples_generator import make_sparse_uncorrelated
from sklearn.datasets.samples_generator import make_regression
from sklearn.datasets import load_iris

rng = np.random.RandomState(0)
rtol = 1e-6


def test_linear_regression():
    # Test LinearRegression on a simple dataset.
    # a simple dataset
    X = [[1], [2]]
    Y = [1, 2]

    reg = LinearRegression()
    reg.fit(X, Y)

    assert_array_almost_equal(reg.coef_, [1])
    assert_array_almost_equal(reg.intercept_, [0])
    assert_array_almost_equal(reg.predict(X), [1, 2])

    # test it also for degenerate input
    X = [[1]]
    Y = [0]

    reg = LinearRegression()
    reg.fit(X, Y)
    assert_array_almost_equal(reg.coef_, [0])
    assert_array_almost_equal(reg.intercept_, [0])
    assert_array_almost_equal(reg.predict(X), [0])


def test_linear_regression_sample_weights():
    # TODO: loop over sparse data as well

    rng = np.random.RandomState(0)

    # It would not work with under-determined systems
    for n_samples, n_features in ((6, 5), ):

        y = rng.randn(n_samples)
        X = rng.randn(n_samples, n_features)
        sample_weight = 1.0 + rng.rand(n_samples)

        for intercept in (True, False):

            # LinearRegression with explicit sample_weight
            reg = LinearRegression(fit_intercept=intercept)
            reg.fit(X, y, sample_weight=sample_weight)
            coefs1 = reg.coef_
            inter1 = reg.intercept_

            assert reg.coef_.shape == (X.shape[1], )  # sanity checks
            assert reg.score(X, y) > 0.5

            # Closed form of the weighted least square
            # theta = (X^T W X)^(-1) * X^T W y
            W = np.diag(sample_weight)
            if intercept is False:
                X_aug = X
            else:
                dummy_column = np.ones(shape=(n_samples, 1))
                X_aug = np.concatenate((dummy_column, X), axis=1)

            coefs2 = linalg.solve(X_aug.T.dot(W).dot(X_aug),
                                  X_aug.T.dot(W).dot(y))

            if intercept is False:
                assert_array_almost_equal(coefs1, coefs2)
            else:
                assert_array_almost_equal(coefs1, coefs2[1:])
                assert_almost_equal(inter1, coefs2[0])


def test_raises_value_error_if_sample_weights_greater_than_1d():
    # Sample weights must be either scalar or 1D

    n_sampless = [2, 3]
    n_featuress = [3, 2]

    for n_samples, n_features in zip(n_sampless, n_featuress):
        X = rng.randn(n_samples, n_features)
        y = rng.randn(n_samples)
        sample_weights_OK = rng.randn(n_samples) ** 2 + 1
        sample_weights_OK_1 = 1.
        sample_weights_OK_2 = 2.

        reg = LinearRegression()

        # make sure the "OK" sample weights actually work
        reg.fit(X, y, sample_weights_OK)
        reg.fit(X, y, sample_weights_OK_1)
        reg.fit(X, y, sample_weights_OK_2)


def test_fit_intercept():
    # Test assertions on betas shape.
    X2 = np.array([[0.38349978, 0.61650022],
                   [0.58853682, 0.41146318]])
    X3 = np.array([[0.27677969, 0.70693172, 0.01628859],
                   [0.08385139, 0.20692515, 0.70922346]])
    y = np.array([1, 1])

    lr2_without_intercept = LinearRegression(fit_intercept=False).fit(X2, y)
    lr2_with_intercept = LinearRegression().fit(X2, y)

    lr3_without_intercept = LinearRegression(fit_intercept=False).fit(X3, y)
    lr3_with_intercept = LinearRegression().fit(X3, y)

    assert (lr2_with_intercept.coef_.shape ==
                 lr2_without_intercept.coef_.shape)
    assert (lr3_with_intercept.coef_.shape ==
                 lr3_without_intercept.coef_.shape)
    assert (lr2_without_intercept.coef_.ndim ==
                 lr3_without_intercept.coef_.ndim)


def test_linear_regression_sparse(random_state=0):
    # Test that linear regression also works with sparse data
    random_state = check_random_state(random_state)
    for i in range(10):
        n = 100
        X = sparse.eye(n, n)
        beta = random_state.rand(n)
        y = X * beta[:, np.newaxis]

        ols = LinearRegression()
        ols.fit(X, y.ravel())
        assert_array_almost_equal(beta, ols.coef_ + ols.intercept_)

        assert_array_almost_equal(ols.predict(X) - y.ravel(), 0)


@pytest.mark.parametrize('normalize', [True, False])
@pytest.mark.parametrize('fit_intercept', [True, False])
def test_linear_regression_sparse_equal_dense(normalize, fit_intercept):
    # Test that linear regression agrees between sparse and dense
    rng = check_random_state(0)
    n_samples = 200
    n_features = 2
    X = rng.randn(n_samples, n_features)
    X[X < 0.1] = 0.
    Xcsr = sparse.csr_matrix(X)
    y = rng.rand(n_samples)
    params = dict(normalize=normalize, fit_intercept=fit_intercept)
    clf_dense = LinearRegression(**params)
    clf_sparse = LinearRegression(**params)
    clf_dense.fit(X, y)
    clf_sparse.fit(Xcsr, y)
    assert clf_dense.intercept_ == pytest.approx(clf_sparse.intercept_)
    assert_allclose(clf_dense.coef_, clf_sparse.coef_)


def test_linear_regression_multiple_outcome(random_state=0):
    # Test multiple-outcome linear regressions
    X, y = make_regression(random_state=random_state)

    Y = np.vstack((y, y)).T
    n_features = X.shape[1]

    reg = LinearRegression()
    reg.fit((X), Y)
    assert reg.coef_.shape == (2, n_features)
    Y_pred = reg.predict(X)
    reg.fit(X, y)
    y_pred = reg.predict(X)
    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T, Y_pred, decimal=3)


def test_linear_regression_sparse_multiple_outcome(random_state=0):
    # Test multiple-outcome linear regressions with sparse data
    random_state = check_random_state(random_state)
    X, y = make_sparse_uncorrelated(random_state=random_state)
    X = sparse.coo_matrix(X)
    Y = np.vstack((y, y)).T
    n_features = X.shape[1]

    ols = LinearRegression()
    ols.fit(X, Y)
    assert ols.coef_.shape == (2, n_features)
    Y_pred = ols.predict(X)
    ols.fit(X, y.ravel())
    y_pred = ols.predict(X)
    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T, Y_pred, decimal=3)


def test_preprocess_data():
    n_samples = 200
    n_features = 2
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    expected_X_mean = np.mean(X, axis=0)
    expected_X_norm = np.std(X, axis=0) * np.sqrt(X.shape[0])
    expected_y_mean = np.mean(y, axis=0)

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=False, normalize=False)
    assert_array_almost_equal(X_mean, np.zeros(n_features))
    assert_array_almost_equal(y_mean, 0)
    assert_array_almost_equal(X_norm, np.ones(n_features))
    assert_array_almost_equal(Xt, X)
    assert_array_almost_equal(yt, y)

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=True, normalize=False)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_norm, np.ones(n_features))
    assert_array_almost_equal(Xt, X - expected_X_mean)
    assert_array_almost_equal(yt, y - expected_y_mean)

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=True, normalize=True)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_norm, expected_X_norm)
    assert_array_almost_equal(Xt, (X - expected_X_mean) / expected_X_norm)
    assert_array_almost_equal(yt, y - expected_y_mean)


def test_preprocess_data_multioutput():
    n_samples = 200
    n_features = 3
    n_outputs = 2
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples, n_outputs)
    expected_y_mean = np.mean(y, axis=0)

    args = [X, sparse.csc_matrix(X)]
    for X in args:
        _, yt, _, y_mean, _ = _preprocess_data(X, y, fit_intercept=False,
                                               normalize=False)
        assert_array_almost_equal(y_mean, np.zeros(n_outputs))
        assert_array_almost_equal(yt, y)

        _, yt, _, y_mean, _ = _preprocess_data(X, y, fit_intercept=True,
                                               normalize=False)
        assert_array_almost_equal(y_mean, expected_y_mean)
        assert_array_almost_equal(yt, y - y_mean)

        _, yt, _, y_mean, _ = _preprocess_data(X, y, fit_intercept=True,
                                               normalize=True)
        assert_array_almost_equal(y_mean, expected_y_mean)
        assert_array_almost_equal(yt, y - y_mean)


def test_preprocess_data_weighted():
    n_samples = 200
    n_features = 2
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    sample_weight = rng.rand(n_samples)
    expected_X_mean = np.average(X, axis=0, weights=sample_weight)
    expected_y_mean = np.average(y, axis=0, weights=sample_weight)

    # XXX: if normalize=True, should we expect a weighted standard deviation?
    #      Currently not weighted, but calculated with respect to weighted mean
    expected_X_norm = (np.sqrt(X.shape[0]) *
                       np.mean((X - expected_X_mean) ** 2, axis=0) ** .5)

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=True, normalize=False,
                         sample_weight=sample_weight)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_norm, np.ones(n_features))
    assert_array_almost_equal(Xt, X - expected_X_mean)
    assert_array_almost_equal(yt, y - expected_y_mean)

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=True, normalize=True,
                         sample_weight=sample_weight)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_norm, expected_X_norm)
    assert_array_almost_equal(Xt, (X - expected_X_mean) / expected_X_norm)
    assert_array_almost_equal(yt, y - expected_y_mean)


def test_sparse_preprocess_data_with_return_mean():
    n_samples = 200
    n_features = 2
    # random_state not supported yet in sparse.rand
    X = sparse.rand(n_samples, n_features, density=.5)  # , random_state=rng
    X = X.tolil()
    y = rng.rand(n_samples)
    XA = X.toarray()
    expected_X_norm = np.std(XA, axis=0) * np.sqrt(X.shape[0])

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=False, normalize=False,
                         return_mean=True)
    assert_array_almost_equal(X_mean, np.zeros(n_features))
    assert_array_almost_equal(y_mean, 0)
    assert_array_almost_equal(X_norm, np.ones(n_features))
    assert_array_almost_equal(Xt.A, XA)
    assert_array_almost_equal(yt, y)

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=True, normalize=False,
                         return_mean=True)
    assert_array_almost_equal(X_mean, np.mean(XA, axis=0))
    assert_array_almost_equal(y_mean, np.mean(y, axis=0))
    assert_array_almost_equal(X_norm, np.ones(n_features))
    assert_array_almost_equal(Xt.A, XA)
    assert_array_almost_equal(yt, y - np.mean(y, axis=0))

    Xt, yt, X_mean, y_mean, X_norm = \
        _preprocess_data(X, y, fit_intercept=True, normalize=True,
                         return_mean=True)
    assert_array_almost_equal(X_mean, np.mean(XA, axis=0))
    assert_array_almost_equal(y_mean, np.mean(y, axis=0))
    assert_array_almost_equal(X_norm, expected_X_norm)
    assert_array_almost_equal(Xt.A, XA / expected_X_norm)
    assert_array_almost_equal(yt, y - np.mean(y, axis=0))


def test_csr_preprocess_data():
    # Test output format of _preprocess_data, when input is csr
    X, y = make_regression()
    X[X < 2.5] = 0.0
    csr = sparse.csr_matrix(X)
    csr_, y, _, _, _ = _preprocess_data(csr, y, True)
    assert csr_.getformat() == 'csr'


@pytest.mark.parametrize('is_sparse', (True, False))
@pytest.mark.parametrize('to_copy', (True, False))
def test_preprocess_copy_data_no_checks(is_sparse, to_copy):
    X, y = make_regression()
    X[X < 2.5] = 0.0

    if is_sparse:
        X = sparse.csr_matrix(X)

    X_, y_, _, _, _ = _preprocess_data(X, y, True,
                                       copy=to_copy, check_input=False)

    if to_copy and is_sparse:
        assert not np.may_share_memory(X_.data, X.data)
    elif to_copy:
        assert not np.may_share_memory(X_, X)
    elif is_sparse:
        assert np.may_share_memory(X_.data, X.data)
    else:
        assert np.may_share_memory(X_, X)


def test_dtype_preprocess_data():
    n_samples = 200
    n_features = 2
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)

    X_32 = np.asarray(X, dtype=np.float32)
    y_32 = np.asarray(y, dtype=np.float32)
    X_64 = np.asarray(X, dtype=np.float64)
    y_64 = np.asarray(y, dtype=np.float64)

    for fit_intercept in [True, False]:
        for normalize in [True, False]:

            Xt_32, yt_32, X_mean_32, y_mean_32, X_norm_32 = _preprocess_data(
                X_32, y_32, fit_intercept=fit_intercept, normalize=normalize,
                return_mean=True)

            Xt_64, yt_64, X_mean_64, y_mean_64, X_norm_64 = _preprocess_data(
                X_64, y_64, fit_intercept=fit_intercept, normalize=normalize,
                return_mean=True)

            Xt_3264, yt_3264, X_mean_3264, y_mean_3264, X_norm_3264 = (
                _preprocess_data(X_32, y_64, fit_intercept=fit_intercept,
                                 normalize=normalize, return_mean=True))

            Xt_6432, yt_6432, X_mean_6432, y_mean_6432, X_norm_6432 = (
                _preprocess_data(X_64, y_32, fit_intercept=fit_intercept,
                                 normalize=normalize, return_mean=True))

            assert Xt_32.dtype == np.float32
            assert yt_32.dtype == np.float32
            assert X_mean_32.dtype == np.float32
            assert y_mean_32.dtype == np.float32
            assert X_norm_32.dtype == np.float32

            assert Xt_64.dtype == np.float64
            assert yt_64.dtype == np.float64
            assert X_mean_64.dtype == np.float64
            assert y_mean_64.dtype == np.float64
            assert X_norm_64.dtype == np.float64

            assert Xt_3264.dtype == np.float32
            assert yt_3264.dtype == np.float32
            assert X_mean_3264.dtype == np.float32
            assert y_mean_3264.dtype == np.float32
            assert X_norm_3264.dtype == np.float32

            assert Xt_6432.dtype == np.float64
            assert yt_6432.dtype == np.float64
            assert X_mean_6432.dtype == np.float64
            assert y_mean_6432.dtype == np.float64
            assert X_norm_6432.dtype == np.float64

            assert X_32.dtype == np.float32
            assert y_32.dtype == np.float32
            assert X_64.dtype == np.float64
            assert y_64.dtype == np.float64

            assert_array_almost_equal(Xt_32, Xt_64)
            assert_array_almost_equal(yt_32, yt_64)
            assert_array_almost_equal(X_mean_32, X_mean_64)
            assert_array_almost_equal(y_mean_32, y_mean_64)
            assert_array_almost_equal(X_norm_32, X_norm_64)


def test_rescale_data():
    n_samples = 200
    n_features = 2

    sample_weight = 1.0 + rng.rand(n_samples)
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    rescaled_X, rescaled_y = _rescale_data(X, y, sample_weight)
    rescaled_X2 = X * np.sqrt(sample_weight)[:, np.newaxis]
    rescaled_y2 = y * np.sqrt(sample_weight)
    assert_array_almost_equal(rescaled_X, rescaled_X2)
    assert_array_almost_equal(rescaled_y, rescaled_y2)


def test_fused_types_make_dataset():
    iris = load_iris()

    X_32 = iris.data.astype(np.float32)
    y_32 = iris.target.astype(np.float32)
    X_csr_32 = sparse.csr_matrix(X_32)
    sample_weight_32 = np.arange(y_32.size, dtype=np.float32)

    X_64 = iris.data.astype(np.float64)
    y_64 = iris.target.astype(np.float64)
    X_csr_64 = sparse.csr_matrix(X_64)
    sample_weight_64 = np.arange(y_64.size, dtype=np.float64)

    # array
    dataset_32, _ = make_dataset(X_32, y_32, sample_weight_32)
    dataset_64, _ = make_dataset(X_64, y_64, sample_weight_64)
    xi_32, yi_32, _, _ = dataset_32._next_py()
    xi_64, yi_64, _, _ = dataset_64._next_py()
    xi_data_32, _, _ = xi_32
    xi_data_64, _, _ = xi_64

    assert xi_data_32.dtype == np.float32
    assert xi_data_64.dtype == np.float64
    assert_allclose(yi_64, yi_32, rtol=rtol)

    # csr
    datasetcsr_32, _ = make_dataset(X_csr_32, y_32, sample_weight_32)
    datasetcsr_64, _ = make_dataset(X_csr_64, y_64, sample_weight_64)
    xicsr_32, yicsr_32, _, _ = datasetcsr_32._next_py()
    xicsr_64, yicsr_64, _, _ = datasetcsr_64._next_py()
    xicsr_data_32, _, _ = xicsr_32
    xicsr_data_64, _, _ = xicsr_64

    assert xicsr_data_32.dtype == np.float32
    assert xicsr_data_64.dtype == np.float64

    assert_allclose(xicsr_data_64, xicsr_data_32, rtol=rtol)
    assert_allclose(yicsr_64, yicsr_32, rtol=rtol)

    assert_array_equal(xi_data_32, xicsr_data_32)
    assert_array_equal(xi_data_64, xicsr_data_64)
    assert_array_equal(yi_32, yicsr_32)
    assert_array_equal(yi_64, yicsr_64)
