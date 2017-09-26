# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD 3 clause

import numpy as np
from scipy import sparse
from scipy import linalg
from itertools import product


from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import ignore_warnings

from sklearn.linear_model.base import LinearRegression
from sklearn.linear_model.base import _preprocess_data
from sklearn.linear_model.base import sparse_center_data, center_data
from sklearn.linear_model.base import _rescale_data
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_greater
from sklearn.datasets.samples_generator import make_sparse_uncorrelated
from sklearn.datasets.samples_generator import make_regression

rng = np.random.RandomState(0)


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

            assert_equal(reg.coef_.shape, (X.shape[1], ))  # sanity checks
            assert_greater(reg.score(X, y), 0.5)

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
    lr2_with_intercept = LinearRegression(fit_intercept=True).fit(X2, y)

    lr3_without_intercept = LinearRegression(fit_intercept=False).fit(X3, y)
    lr3_with_intercept = LinearRegression(fit_intercept=True).fit(X3, y)

    assert_equal(lr2_with_intercept.coef_.shape,
                 lr2_without_intercept.coef_.shape)
    assert_equal(lr3_with_intercept.coef_.shape,
                 lr3_without_intercept.coef_.shape)
    assert_equal(lr2_without_intercept.coef_.ndim,
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


def test_linear_regression_multiple_outcome(random_state=0):
    # Test multiple-outcome linear regressions
    X, y = make_regression(random_state=random_state)

    Y = np.vstack((y, y)).T
    n_features = X.shape[1]

    reg = LinearRegression(fit_intercept=True)
    reg.fit((X), Y)
    assert_equal(reg.coef_.shape, (2, n_features))
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
    assert_equal(ols.coef_.shape, (2, n_features))
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
    assert_equal(csr_.getformat(), 'csr')


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

            assert_equal(Xt_32.dtype, np.float32)
            assert_equal(yt_32.dtype, np.float32)
            assert_equal(X_mean_32.dtype, np.float32)
            assert_equal(y_mean_32.dtype, np.float32)
            assert_equal(X_norm_32.dtype, np.float32)

            assert_equal(Xt_64.dtype, np.float64)
            assert_equal(yt_64.dtype, np.float64)
            assert_equal(X_mean_64.dtype, np.float64)
            assert_equal(y_mean_64.dtype, np.float64)
            assert_equal(X_norm_64.dtype, np.float64)

            assert_equal(Xt_3264.dtype, np.float32)
            assert_equal(yt_3264.dtype, np.float32)
            assert_equal(X_mean_3264.dtype, np.float32)
            assert_equal(y_mean_3264.dtype, np.float32)
            assert_equal(X_norm_3264.dtype, np.float32)

            assert_equal(Xt_6432.dtype, np.float64)
            assert_equal(yt_6432.dtype, np.float64)
            assert_equal(X_mean_6432.dtype, np.float64)
            assert_equal(y_mean_6432.dtype, np.float64)
            assert_equal(X_norm_6432.dtype, np.float64)

            assert_equal(X_32.dtype, np.float32)
            assert_equal(y_32.dtype, np.float32)
            assert_equal(X_64.dtype, np.float64)
            assert_equal(y_64.dtype, np.float64)

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


@ignore_warnings  # all deprecation warnings
def test_deprecation_center_data():
    n_samples = 200
    n_features = 2

    w = 1.0 + rng.rand(n_samples)
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)

    param_grid = product([True, False], [True, False], [True, False],
                         [None, w])

    for (fit_intercept, normalize, copy, sample_weight) in param_grid:

        XX = X.copy()  # such that we can try copy=False as well

        X1, y1, X1_mean, X1_var, y1_mean = \
            center_data(XX, y, fit_intercept=fit_intercept,
                        normalize=normalize, copy=copy,
                        sample_weight=sample_weight)

        XX = X.copy()

        X2, y2, X2_mean, X2_var, y2_mean = \
            _preprocess_data(XX, y, fit_intercept=fit_intercept,
                             normalize=normalize, copy=copy,
                             sample_weight=sample_weight)

        assert_array_almost_equal(X1, X2)
        assert_array_almost_equal(y1, y2)
        assert_array_almost_equal(X1_mean, X2_mean)
        assert_array_almost_equal(X1_var, X2_var)
        assert_array_almost_equal(y1_mean, y2_mean)

    # Sparse cases
    X = sparse.csr_matrix(X)

    for (fit_intercept, normalize, copy, sample_weight) in param_grid:

        X1, y1, X1_mean, X1_var, y1_mean = \
            center_data(X, y, fit_intercept=fit_intercept, normalize=normalize,
                        copy=copy, sample_weight=sample_weight)

        X2, y2, X2_mean, X2_var, y2_mean = \
            _preprocess_data(X, y, fit_intercept=fit_intercept,
                             normalize=normalize, copy=copy,
                             sample_weight=sample_weight, return_mean=False)

        assert_array_almost_equal(X1.toarray(), X2.toarray())
        assert_array_almost_equal(y1, y2)
        assert_array_almost_equal(X1_mean, X2_mean)
        assert_array_almost_equal(X1_var, X2_var)
        assert_array_almost_equal(y1_mean, y2_mean)

    for (fit_intercept, normalize) in product([True, False], [True, False]):

        X1, y1, X1_mean, X1_var, y1_mean = \
            sparse_center_data(X, y, fit_intercept=fit_intercept,
                               normalize=normalize)

        X2, y2, X2_mean, X2_var, y2_mean = \
            _preprocess_data(X, y, fit_intercept=fit_intercept,
                             normalize=normalize, return_mean=True)

        assert_array_almost_equal(X1.toarray(), X2.toarray())
        assert_array_almost_equal(y1, y2)
        assert_array_almost_equal(X1_mean, X2_mean)
        assert_array_almost_equal(X1_var, X2_var)
        assert_array_almost_equal(y1_mean, y2_mean)
