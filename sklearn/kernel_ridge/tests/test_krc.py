"""
Testing for Kernel Ridge Classifier
(sklearn.kernel_ridge)

Author: Carlos Perales <sir.perales@gmail.com>
License: BSD 3 clause

"""
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.testing import assert_almost_equal
from sklearn import datasets, base
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises
from sklearn.externals import six
from sklearn.kernel_ridge import KernelRidgeClassifier

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


def test_tweak_params():
    # Make sure some tweaking of parameters works.
    clf = KernelRidgeClassifier(kernel='linear', gamma=0.0, alpha=0.5)
    clf.fit(X, Y)
    dual_coef = np.array([0.13333333, -0.57777778, 0.13333333,
                          0.57777778, -0.13333333, -0.13333333])
    assert_almost_equal(clf.dual_coef_, dual_coef, decimal=8)
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    dual_coef = np.array([-1.06666667, -0.71111111, -1.06666667,
                          -1.28888889, -3.93333333, -0.93333333])
    clf.dual_coef_ = dual_coef
    assert_almost_equal(clf.predict([[-.1, -.1]]), [2])


def test_bad_input():
    # Test that it gives proper exception on deficient input
    # impossible value of C

    assert_raises(ValueError, KernelRidgeClassifier(alpha=-0.5).fit, X, Y)

    clf = KernelRidgeClassifier()
    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)

    # error for precomputed kernels
    clf = KernelRidgeClassifier(kernel='precomputed')
    assert_raises(ValueError, clf.fit, X, Y)

    Xt = np.array(X).T
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = KernelRidgeClassifier()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)


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
    KernelRidgeClassifier_callable = \
        KernelRidgeClassifier(kernel=lambda x, y: np.dot(x, y.T))
    # clone for checking clonability with lambda functions..
    KernelRidgeClassifier_cloned = base.clone(KernelRidgeClassifier_callable)
    KernelRidgeClassifier_cloned.fit(iris.data, iris.target)

    KernelRidgeClassifier_builtin = KernelRidgeClassifier(kernel='linear')
    KernelRidgeClassifier_builtin.fit(iris.data, iris.target)

    assert_array_almost_equal(KernelRidgeClassifier_cloned.dual_coef_,
                              KernelRidgeClassifier_builtin.dual_coef_)
    assert_array_equal(KernelRidgeClassifier_cloned.predict(iris.data),
                       KernelRidgeClassifier_builtin.predict(iris.data))
    assert_array_almost_equal(KernelRidgeClassifier_cloned.predict(iris.data),
                              KernelRidgeClassifier_builtin.predict(iris.data))


def test_krc_bad_kernel():
    clf = KernelRidgeClassifier(kernel=lambda x, y: x)
    assert_raises(ValueError, clf.fit, X, Y)
