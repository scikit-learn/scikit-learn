"""
Testing for Kernel Ridge Classifier
(sklearn.kernel_ridge)

Authors: Carlos Perales <sir.perales@gmail.com>
License: BSD 3 clause

"""
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.testing import assert_almost_equal
from sklearn import datasets, base
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal, assert_true
from sklearn.utils.testing import assert_greater
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


def test_kernel_ridge_classifier_parameters():
    # Test parameters.
    clf = KernelRidgeClassifier(kernel='rbf').fit(X, Y)
    assert_array_equal(clf.predict(X), Y)


def test_precomputed():
    # KernelRidgeClassifier with a precomputed kernel.
    # # We test it with check_estimator first.
    # Then we test it with a toy dataset and with iris.
    clf = KernelRidgeClassifier(kernel='precomputed')
    # Gram matrix for train data (square matrix)
    # (we use just a linear kernel)
    K = np.dot(X, np.array(X).T)
    clf.fit(K, Y)
    # Gram matrix for test data (rectangular matrix)
    KT = np.dot(T, np.array(X).T)
    pred = clf.predict(KT)
    assert_raises(ValueError, clf.predict, KT.T)

    pred = clf.predict(KT)
    assert_array_equal(pred, true_result)

    # same as before, but using a callable function instead of the kernel
    # matrix. kernel is just a linear kernel
    def kfunc(x, y):
        return np.dot(x, y.T)
    clf = KernelRidgeClassifier(kernel=kfunc)
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_equal(pred, true_result)

    # test a precomputed kernel with the iris dataset
    # and check parameters against a linear KernelRidgeClassifier
    clf = KernelRidgeClassifier(kernel='precomputed')
    clf2 = KernelRidgeClassifier(kernel='linear')
    K = np.dot(iris.data, iris.data.T)
    clf.fit(K, iris.target)
    clf2.fit(iris.data, iris.target)
    pred = clf.predict(K)
    assert_array_almost_equal(clf.dual_coef_, clf2.dual_coef_)
    assert clf.gamma == clf2.gamma
    assert_almost_equal(np.mean(pred == iris.target), .82, decimal=2)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    K = np.zeros_like(K)
    for i in range(len(iris.data)):
        for j in range(clf.dual_coef_.shape[0]):
            K[i, j] = np.dot(iris.data[i], iris.data[j])

    pred = clf.predict(K)
    assert_almost_equal(np.mean(pred == iris.target), .82, decimal=2)

    clf = KernelRidgeClassifier(kernel=kfunc)
    clf.fit(iris.data, iris.target)
    assert_almost_equal(np.mean(pred == iris.target), .82, decimal=2)


def test_tweak_params():
    # Make sure some tweaking of parameters works.
    clf = KernelRidgeClassifier(kernel='linear', gamma=0.0, alpha=0.5)
    clf.fit(X, Y)
    dual_coef = np.array([1.06666667, 0.71111111, 1.06666667,
                          1.28888889, 0.93333333, 0.93333333])
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


def test_linear_elm():
    clf = KernelRidgeClassifier(kernel='linear').fit(X, Y)

    # by default should have dual_coef
    assert_true(clf.dual_coef_.all())

    assert_array_equal(clf.predict(T), true_result)


def test_linear_elm_iris():
    # Test that Kernel ELM with linear kernel
    # gives plausible predictions on the iris dataset
    # Also, test symbolic class names (classes_).
    target = iris.target_names[iris.target]
    clf = KernelRidgeClassifier(kernel='linear').fit(iris.data, target)
    assert_equal(set(clf.classes_), set(iris.target_names))
    assert_greater(np.mean(clf.predict(iris.data) == target), 0.79)


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
