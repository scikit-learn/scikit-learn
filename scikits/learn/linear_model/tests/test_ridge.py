import numpy as np

from numpy.testing import assert_almost_equal, assert_array_almost_equal, \
                          assert_equal, assert_array_equal

from scikits.learn import datasets
from scikits.learn.metrics import mean_square_error

from scikits.learn.linear_model.base import LinearRegression

from scikits.learn.linear_model.ridge import Ridge
from scikits.learn.linear_model.sparse.ridge import Ridge as SpRidge

from scikits.learn.linear_model.ridge import RidgeLOO
from scikits.learn.linear_model.sparse.ridge import RidgeLOO as SpRidgeLOO

diabetes = datasets.load_diabetes()

np.random.seed(0)

def test_ridge():
    """Ridge regression convergence test using score

    TODO: for this test to be robust, we should use a dataset instead
    of np.random.
    """
    alpha = 1.0

    # With more samples than features
    n_samples, n_features = 6, 5
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
    assert ridge.score(X, y) > 0.5

    # With more features than samples
    n_samples, n_features = 5, 10
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
    assert ridge.score(X, y) > .9

def test_toy_ridge_object():
    """Test BayesianRegression ridge classifier

    TODO: test also n_samples > n_features
    """
    X = np.array([[1], [2]])
    Y = np.array([1, 2])
    clf = Ridge(alpha=0.0)
    clf.fit(X, Y)
    X_test = [[1], [2], [3], [4]]
    assert_almost_equal(clf.predict(X_test), [1., 2, 3, 4])

    assert_equal(len(clf.coef_.shape), 1)
    assert_equal(type(clf.intercept_), np.float64)

    Y = np.vstack((Y,Y)).T

    clf.fit(X, Y)
    X_test = [[1], [2], [3], [4]]

    assert_equal(len(clf.coef_.shape), 2)
    assert_equal(type(clf.intercept_), np.ndarray)


def test_ridge_vs_lstsq():
    """On alpha=0., Ridge and OLS yield the same solution."""

    # we need more samples than features
    n_samples, n_features = 5, 4
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    ridge = Ridge(alpha=0.)
    ols = LinearRegression()

    ridge.fit(X, y)
    ols.fit (X, y)
    assert_almost_equal(ridge.coef_, ols.coef_)

    ridge.fit(X, y, fit_intercept=False)
    ols.fit (X, y, fit_intercept=False)
    assert_almost_equal(ridge.coef_, ols.coef_)

def _test_ridge_loo(ridge_class, ridge_loo_class):
    X, y = diabetes.data, diabetes.target
    n_samples = X.shape[0]
    ind = np.arange(n_samples)
    np.random.shuffle(ind)
    ind = ind[:200]
    X, y = X[ind], y[ind]
    n_samples = X.shape[0]

    ret = []

    ridge_loo = ridge_loo_class(fit_intercept=False)
    ridge = ridge_class(fit_intercept=False)

    # efficient LOO
    K, v, Q = ridge_loo._pre_compute(X, y)
    errors, c = ridge_loo._errors(v, Q, y, 1.0)
    values, c = ridge_loo._values(K, v, Q, y, 1.0)

    # brute-force LOO: remove one example at a time
    errors2 = []
    values2 = []
    for i in range(n_samples):
        sel = np.arange(n_samples) != i
        X_new = X[sel]
        y_new = y[sel]
        ridge.fit(X_new, y_new)
        value = ridge.predict([X[i]])[0]
        error = (y[i] - value) ** 2
        errors2.append(error)
        values2.append(value)

    # check that efficient and brute-force LOO give same results
    assert_almost_equal(errors, errors2)
    assert_almost_equal(values, values2)

    # check best alpha
    ridge_loo.fit(X, y)
    best_alpha = ridge_loo.best_alpha
    ret.append(best_alpha)

    # check that we get same best alpha with custom loss_func
    ridge_loo2 = RidgeLOO(fit_intercept=False, loss_func=mean_square_error)
    ridge_loo2.fit(X, y)
    assert_equal(ridge_loo2.best_alpha, best_alpha)

    # check that we get same best alpha with sample weights
    ridge_loo.fit(X, y, sample_weight=np.ones(n_samples))
    assert_equal(ridge_loo.best_alpha, best_alpha)

    # simulate several responses
    Y = np.vstack((y,y)).T

    ridge_loo.fit(X, Y)
    Y_pred = ridge_loo.predict(X)
    ridge_loo.fit(X, y)
    y_pred = ridge_loo.predict(X)

    assert_array_almost_equal(np.vstack((y_pred,y_pred)).T,
                              Y_pred)

def test_ridge_loo():
    ret_dense = _test_ridge_loo(Ridge, RidgeLOO)
    ret_sp = _test_ridge_loo(SpRidge, SpRidgeLOO)
    assert_array_equal(ret_dense, ret_sp)

