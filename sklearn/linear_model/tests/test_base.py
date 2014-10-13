# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#
# License: BSD 3 clause

import numpy as np
from scipy import sparse

from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal

from sklearn.linear_model.base import LinearRegression
from sklearn.linear_model.base import center_data, sparse_center_data
from sklearn.utils import check_random_state
from sklearn.datasets.samples_generator import make_sparse_uncorrelated
from sklearn.datasets.samples_generator import make_regression


def test_linear_regression():
    """
    Test LinearRegression on a simple dataset.
    """
    # a simple dataset
    X = [[1], [2]]
    Y = [1, 2]

    clf = LinearRegression()
    clf.fit(X, Y)

    assert_array_almost_equal(clf.coef_, [1])
    assert_array_almost_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.predict(X), [1, 2])

    # test it also for degenerate input
    X = [[1]]
    Y = [0]

    clf = LinearRegression()
    clf.fit(X, Y)
    assert_array_almost_equal(clf.coef_, [0])
    assert_array_almost_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.predict(X), [0])


def test_linear_regression_n_jobs():
    """
    Test for the n_jobs parameter on the fit method and the constructor
    """
    X = [[1], [2]]
    Y = [1, 2]
    clf = LinearRegression()
    clf_fit = clf.fit(X, Y, 4)
    assert_equal(clf_fit.n_jobs, clf.n_jobs)
    assert_equal(clf.n_jobs, 1)


def test_fit_intercept():
    """
    Test assertions on betas shape.
    """
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
    "Test that linear regression also works with sparse data"
    random_state = check_random_state(random_state)
    n = 100
    X = sparse.eye(n, n)
    beta = random_state.rand(n)
    y = X * beta[:, np.newaxis]

    ols = LinearRegression()
    ols.fit(X, y.ravel())
    assert_array_almost_equal(beta, ols.coef_ + ols.intercept_)
    assert_array_almost_equal(ols.residues_, 0)


def test_linear_regression_multiple_outcome(random_state=0):
    "Test multiple-outcome linear regressions"
    X, y = make_regression(random_state=random_state)

    Y = np.vstack((y, y)).T
    n_features = X.shape[1]

    clf = LinearRegression(fit_intercept=True)
    clf.fit((X), Y)
    assert_equal(clf.coef_.shape, (2, n_features))
    Y_pred = clf.predict(X)
    clf.fit(X, y)
    y_pred = clf.predict(X)
    assert_array_almost_equal(np.vstack((y_pred, y_pred)).T, Y_pred, decimal=3)


def test_linear_regression_sparse_multiple_outcome(random_state=0):
    "Test multiple-outcome linear regressions with sparse data"
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


def test_center_data():
    n_samples = 200
    n_features = 2
    rng = check_random_state(0)
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    expected_X_mean = np.mean(X, axis=0)
    # XXX: currently scaled to variance=n_samples
    expected_X_std = np.std(X, axis=0) * np.sqrt(X.shape[0])
    expected_y_mean = np.mean(y, axis=0)

    Xt, yt, X_mean, y_mean, X_std = center_data(X, y, fit_intercept=False,
                                                normalize=False)
    assert_array_almost_equal(X_mean, np.zeros(n_features))
    assert_array_almost_equal(y_mean, 0)
    assert_array_almost_equal(X_std, np.ones(n_features))
    assert_array_almost_equal(Xt, X)
    assert_array_almost_equal(yt, y)

    Xt, yt, X_mean, y_mean, X_std = center_data(X, y, fit_intercept=True,
                                                normalize=False)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_std, np.ones(n_features))
    assert_array_almost_equal(Xt, X - expected_X_mean)
    assert_array_almost_equal(yt, y - expected_y_mean)

    Xt, yt, X_mean, y_mean, X_std = center_data(X, y, fit_intercept=True,
                                                normalize=True)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_std, expected_X_std)
    assert_array_almost_equal(Xt, (X - expected_X_mean) / expected_X_std)
    assert_array_almost_equal(yt, y - expected_y_mean)


def test_center_data_multioutput():
    n_samples = 200
    n_features = 3
    n_outputs = 2
    rng = check_random_state(0)
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples, n_outputs)
    expected_y_mean = np.mean(y, axis=0)

    args = [(center_data, X), (sparse_center_data, sparse.csc_matrix(X))]
    for center, X in args:
        _, yt, _, y_mean, _ = center(X, y, fit_intercept=False,
                                     normalize=False)
        assert_array_almost_equal(y_mean, np.zeros(n_outputs))
        assert_array_almost_equal(yt, y)

        _, yt, _, y_mean, _ = center(X, y, fit_intercept=True,
                                     normalize=False)
        assert_array_almost_equal(y_mean, expected_y_mean)
        assert_array_almost_equal(yt, y - y_mean)

        _, yt, _, y_mean, _ = center(X, y, fit_intercept=True,
                                     normalize=True)
        assert_array_almost_equal(y_mean, expected_y_mean)
        assert_array_almost_equal(yt, y - y_mean)


def test_center_data_weighted():
    n_samples = 200
    n_features = 2
    rng = check_random_state(0)
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)
    sample_weight = rng.rand(n_samples)
    expected_X_mean = np.average(X, axis=0, weights=sample_weight)
    expected_y_mean = np.average(y, axis=0, weights=sample_weight)

    # XXX: if normalize=True, should we expect a weighted standard deviation?
    #      Currently not weighted, but calculated with respect to weighted mean
    # XXX: currently scaled to variance=n_samples
    expected_X_std = (np.sqrt(X.shape[0]) *
                      np.mean((X - expected_X_mean) ** 2, axis=0) ** .5)

    Xt, yt, X_mean, y_mean, X_std = center_data(X, y, fit_intercept=True,
                                                normalize=False,
                                                sample_weight=sample_weight)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_std, np.ones(n_features))
    assert_array_almost_equal(Xt, X - expected_X_mean)
    assert_array_almost_equal(yt, y - expected_y_mean)

    Xt, yt, X_mean, y_mean, X_std = center_data(X, y, fit_intercept=True,
                                                normalize=True,
                                                sample_weight=sample_weight)
    assert_array_almost_equal(X_mean, expected_X_mean)
    assert_array_almost_equal(y_mean, expected_y_mean)
    assert_array_almost_equal(X_std, expected_X_std)
    assert_array_almost_equal(Xt, (X - expected_X_mean) / expected_X_std)
    assert_array_almost_equal(yt, y - expected_y_mean)


def test_sparse_center_data():
    n_samples = 200
    n_features = 2
    rng = check_random_state(0)
    # random_state not supported yet in sparse.rand
    X = sparse.rand(n_samples, n_features, density=.5)  # , random_state=rng
    X = X.tolil()
    y = rng.rand(n_samples)
    XA = X.toarray()
    # XXX: currently scaled to variance=n_samples
    expected_X_std = np.std(XA, axis=0) * np.sqrt(X.shape[0])

    Xt, yt, X_mean, y_mean, X_std = sparse_center_data(X, y,
                                                       fit_intercept=False,
                                                       normalize=False)
    assert_array_almost_equal(X_mean, np.zeros(n_features))
    assert_array_almost_equal(y_mean, 0)
    assert_array_almost_equal(X_std, np.ones(n_features))
    assert_array_almost_equal(Xt.A, XA)
    assert_array_almost_equal(yt, y)

    Xt, yt, X_mean, y_mean, X_std = sparse_center_data(X, y,
                                                       fit_intercept=True,
                                                       normalize=False)
    assert_array_almost_equal(X_mean, np.mean(XA, axis=0))
    assert_array_almost_equal(y_mean, np.mean(y, axis=0))
    assert_array_almost_equal(X_std, np.ones(n_features))
    assert_array_almost_equal(Xt.A, XA)
    assert_array_almost_equal(yt, y - np.mean(y, axis=0))

    Xt, yt, X_mean, y_mean, X_std = sparse_center_data(X, y,
                                                       fit_intercept=True,
                                                       normalize=True)
    assert_array_almost_equal(X_mean, np.mean(XA, axis=0))
    assert_array_almost_equal(y_mean, np.mean(y, axis=0))
    assert_array_almost_equal(X_std, expected_X_std)
    assert_array_almost_equal(Xt.A, XA / expected_X_std)
    assert_array_almost_equal(yt, y - np.mean(y, axis=0))


def test_csr_sparse_center_data():
    """Test output format of sparse_center_data, when input is csr"""
    X, y = make_regression()
    X[X < 2.5] = 0.0
    csr = sparse.csr_matrix(X)
    csr_, y, _, _, _ = sparse_center_data(csr, y, True)
    assert_equal(csr_.getformat(), 'csr')
