import numpy as np
import scipy.sparse as sp

from numpy.testing import assert_almost_equal, assert_array_almost_equal, \
                          assert_equal, assert_array_equal

from scikits.learn import datasets
from scikits.learn.metrics import mean_square_error

from scikits.learn.linear_model.base import LinearRegression

from scikits.learn.linear_model.ridge import Ridge
from scikits.learn.linear_model.ridge import _RidgeGCV
from scikits.learn.linear_model.ridge import RidgeCV
from scikits.learn.linear_model.ridge import RidgeClassifier
from scikits.learn.linear_model.ridge import RidgeClassifierCV


from scikits.learn.cross_val import KFold

diabetes = datasets.load_diabetes()

X_diabetes, y_diabetes = diabetes.data, diabetes.target
ind = np.arange(X_diabetes.shape[0])
np.random.shuffle(ind)
ind = ind[:200]
X_diabetes, y_diabetes = X_diabetes[ind], y_diabetes[ind]

iris = datasets.load_iris()

X_iris = sp.csr_matrix(iris.data)
y_iris = iris.target

np.random.seed(0)

DENSE_FILTER = lambda X: X
SPARSE_FILTER = lambda X: sp.csr_matrix(X)

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

    ridge.fit(X, y, sample_weight=np.ones(n_samples))
    assert ridge.score(X, y) > 0.5

    # With more features than samples
    n_samples, n_features = 5, 10
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)
    ridge = Ridge(alpha=alpha)
    ridge.fit(X, y)
    assert ridge.score(X, y) > .9

    ridge.fit(X, y, sample_weight=np.ones(n_samples))
    assert ridge.score(X, y) > 0.9

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

def _test_ridge_loo(filter_):
    # test that can work with both dense or sparse matrices
    n_samples = X_diabetes.shape[0]

    ret = []

    ridge_gcv = _RidgeGCV(fit_intercept=False)
    ridge = Ridge(fit_intercept=False)

    # generalized cross-validation (efficient leave-one-out)
    K, v, Q = ridge_gcv._pre_compute(X_diabetes, y_diabetes)
    errors, c = ridge_gcv._errors(v, Q, y_diabetes, 1.0)
    values, c = ridge_gcv._values(K, v, Q, y_diabetes, 1.0)

    # brute-force leave-one-out: remove one example at a time
    errors2 = []
    values2 = []
    for i in range(n_samples):
        sel = np.arange(n_samples) != i
        X_new = X_diabetes[sel]
        y_new = y_diabetes[sel]
        ridge.fit(X_new, y_new)
        value = ridge.predict([X_diabetes[i]])[0]
        error = (y_diabetes[i] - value) ** 2
        errors2.append(error)
        values2.append(value)

    # check that efficient and brute-force LOO give same results
    assert_almost_equal(errors, errors2)
    assert_almost_equal(values, values2)

    # check best alpha
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes)
    best_alpha = ridge_gcv.best_alpha
    ret.append(best_alpha)

    # check that we get same best alpha with custom loss_func
    ridge_gcv2 = _RidgeGCV(fit_intercept=False, loss_func=mean_square_error)
    ridge_gcv2.fit(filter_(X_diabetes), y_diabetes)
    assert_equal(ridge_gcv2.best_alpha, best_alpha)

    # check that we get same best alpha with sample weights
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes,
                  sample_weight=np.ones(n_samples))
    assert_equal(ridge_gcv.best_alpha, best_alpha)

    # simulate several responses
    Y = np.vstack((y_diabetes,y_diabetes)).T

    ridge_gcv.fit(filter_(X_diabetes), Y)
    Y_pred = ridge_gcv.predict(filter_(X_diabetes))
    ridge_gcv.fit(filter_(X_diabetes), y_diabetes)
    y_pred = ridge_gcv.predict(filter_(X_diabetes))

    assert_array_almost_equal(np.vstack((y_pred,y_pred)).T,
                              Y_pred)

    return ret

def _test_ridge_cv(filter_):
    n_samples = X_diabetes.shape[0]

    ridge_cv = RidgeCV()
    ridge_cv.fit(filter_(X_diabetes), y_diabetes)
    ridge_cv.predict(filter_(X_diabetes))

    assert_equal(len(ridge_cv.coef_.shape), 1)
    assert_equal(type(ridge_cv.intercept_), np.float64)

    cv = KFold(n_samples, 5)
    ridge_cv.fit(filter_(X_diabetes), y_diabetes, cv=cv)
    ridge_cv.predict(filter_(X_diabetes))

    assert_equal(len(ridge_cv.coef_.shape), 1)
    assert_equal(type(ridge_cv.intercept_), np.float64)

def _test_ridge_diabetes(filter_):
    ridge = Ridge(fit_intercept=False)
    ridge.fit(filter_(X_diabetes), y_diabetes)
    return np.round(ridge.score(filter_(X_diabetes), y_diabetes), 5)

def _test_multi_ridge_diabetes(filter_):
    # simulate several responses
    Y = np.vstack((y_diabetes,y_diabetes)).T

    ridge = Ridge(fit_intercept=False)
    ridge.fit(filter_(X_diabetes), Y)
    Y_pred = ridge.predict(filter_(X_diabetes))
    ridge.fit(filter_(X_diabetes), y_diabetes)
    y_pred = ridge.predict(filter_(X_diabetes))
    assert_array_almost_equal(np.vstack((y_pred,y_pred)).T,
                              Y_pred)

def _test_ridge_classifiers(filter_):
    for clf in (RidgeClassifier(), RidgeClassifierCV()):
        clf.fit(filter_(X_iris), y_iris)
        y_pred = clf.predict(filter_(X_iris))
        assert np.mean(y_iris == y_pred) >= 0.8

    clf = RidgeClassifierCV()
    n_samples = X_iris.shape[0]
    cv = KFold(n_samples, 5)
    clf.fit(filter_(X_iris), y_iris, cv=cv)
    y_pred = clf.predict(filter_(X_iris))
    assert np.mean(y_iris == y_pred) >= 0.8

def test_dense_sparse():
    for test_func in (_test_ridge_loo,
                      _test_ridge_cv,
                      _test_ridge_diabetes,
                     _test_multi_ridge_diabetes,
                     _test_ridge_classifiers):
        # test dense matrix
        ret_dense = test_func(DENSE_FILTER)
        # test sparse matrix
        ret_sp = test_func(SPARSE_FILTER)
        # test that the outputs are the same
        assert_array_equal(ret_dense, ret_sp)



