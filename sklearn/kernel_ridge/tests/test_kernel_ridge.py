"""
Testing for Kernel Ridge
(sklearn.kernel_ridge)

Authors: Mathieu Blondel <mathieu@mblondel.org>
         Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
         Carlos Perales <sir.perales@gmail.com>
License: BSD 3 clause

"""
import numpy as np
import scipy.sparse as sp
from numpy.testing import assert_array_equal, assert_array_almost_equal
from sklearn import datasets, base
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises, ignore_warnings
from sklearn.externals import six
from sklearn.kernel_ridge import KernelRidgeClassifier
from sklearn.linear_model import Ridge
from sklearn.kernel_ridge import KernelRidge

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

# also load the iris dataset
iris = datasets.load_iris()
rng = check_random_state(42)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# regression data
X_reg, y_reg = datasets.make_regression(n_features=10, random_state=0)
X_reg_csr = sp.csr_matrix(X_reg)
X_reg_csc = sp.csc_matrix(X_reg)
Y_reg = np.array([y_reg, y_reg]).T


# Tests for classification
def test_unicode_kernel():
    # Test that a unicode kernel name does not cause a TypeError
    if six.PY2:
        # Test unicode (same as str on python3)
        clf = KernelRidgeClassifier(kernel=u'linear')
        clf.fit(X, Y)
        clf.predict(T)

    # Test default behavior on both versions
    clf = KernelRidgeClassifier(kernel='linear')
    clf.fit(X, Y)
    clf.predict(T)


def test_krc_clone_with_callable_kernel():
    # create KernelRidgeClassifier with callable linear kernel,
    # check that results are the same as with built-in linear kernel
    krc_callable = KernelRidgeClassifier(kernel=lambda x, y: np.dot(x, y.T))
    # clone for checking clonability with lambda functions..
    krc_cloned = base.clone(krc_callable)
    krc_cloned.fit(iris.data, iris.target)

    krc_builtin = KernelRidgeClassifier(kernel='linear')
    krc_builtin.fit(iris.data, iris.target)

    assert_array_almost_equal(krc_cloned.dual_coef_,
                              krc_builtin.dual_coef_)
    assert_array_equal(krc_cloned.predict(iris.data),
                       krc_builtin.predict(iris.data))
    assert_array_almost_equal(krc_cloned.predict(iris.data),
                              krc_builtin.predict(iris.data))


def test_krc_bad_kernel():
    clf = KernelRidgeClassifier(kernel=lambda x, y: x)
    assert_raises(ValueError, clf.fit, X, Y)


def test_bad_input():
    # Test that it gives proper exception on negative alpha
    # alpha cannot be negative because Kernel Ridge come from
    # a minimization problem where all the terms must be positive
    assert_raises(ValueError, KernelRidgeClassifier(alpha=-0.5).fit, X, Y)


# Test for regression
def test_kernel_ridge():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X_reg, y_reg).predict(X_reg)
    pred2 = KernelRidge(kernel="linear",
                        alpha=1).fit(X_reg, y_reg).predict(X_reg)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_csr():
    pred = Ridge(alpha=1, fit_intercept=False,
                 solver="cholesky").fit(X_reg_csr, y_reg).predict(X_reg_csr)
    pred2 = KernelRidge(kernel="linear",
                        alpha=1).fit(X_reg_csr, y_reg).predict(X_reg_csr)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_csc():
    pred = Ridge(alpha=1, fit_intercept=False,
                 solver="cholesky").fit(X_reg_csr, y_reg).predict(X_reg_csr)
    pred2 = KernelRidge(kernel="linear",
                        alpha=1).fit(X_reg_csr, y_reg).predict(X_reg_csr)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_singular_kernel():
    # alpha=0 causes a LinAlgError in computing the dual coefficients,
    # which causes a fallback to a lstsq solver. This is tested here.
    pred = Ridge(alpha=0, fit_intercept=False).fit(X_reg, y_reg).predict(X_reg)
    kr = KernelRidge(kernel="linear", alpha=0)
    ignore_warnings(kr.fit)(X_reg, y_reg)
    pred2 = kr.predict(X_reg)
    assert_array_almost_equal(pred, pred2)


def test_kernel_ridge_multi_output():
    pred = Ridge(alpha=1, fit_intercept=False).fit(X_reg, Y_reg).predict(X_reg)
    pred2 = KernelRidge(kernel="linear",
                        alpha=1).fit(X_reg, Y_reg).predict(X_reg)
    assert_array_almost_equal(pred, pred2)

    pred3 = KernelRidge(kernel="linear",
                        alpha=1).fit(X_reg, y_reg).predict(X_reg)
    pred3 = np.array([pred3, pred3]).T
    assert_array_almost_equal(pred2, pred3)
