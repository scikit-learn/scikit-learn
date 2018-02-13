"""
Testing for Kernel Ridge Classification
(sklearn.kernel_ridge)

Authors: Carlos Perales <sir.perales@gmail.com>
License: BSD 3 clause

"""
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.testing import assert_almost_equal
from scipy import sparse
from sklearn import datasets, base
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal, assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_raises
from sklearn.externals import six
from sklearn.utils.estimator_checks import check_estimator
from sklearn.kernel_ridge import KRC

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


def test_check_estimator():
    check_estimator(KRC)


def test_krc_parameters():
    # Test parameters.
    clf = KRC(kernel='linear').fit(X, Y)
    assert_array_equal(clf.predict(X), Y)


def test_iris():
    # Check consistency on dataset iris.

    # shuffle the dataset so that labels are not ordered
    for k in ('linear', 'rbf'):
        clf = KRC(kernel=k).fit(iris.data, iris.target)
        assert_greater(np.mean(clf.predict(iris.data) == iris.target), 0.75)

    assert_array_equal(clf.classes_, np.sort(clf.classes_))


def test_precomputed():
    # KRC with a precomputed kernel.
    # We test it with a toy dataset and with iris.
    clf = KRC(kernel='precomputed')
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
    clf = KRC(kernel=kfunc)
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_equal(pred, true_result)

    # test a precomputed kernel with the iris dataset
    # and check parameters against a linear KRC
    clf = KRC(kernel='precomputed')
    clf2 = KRC(kernel='linear')
    K = np.dot(iris.data, iris.data.T)
    clf.fit(K, iris.target)
    clf2.fit(iris.data, iris.target)
    pred = clf.predict(K)
    assert_array_almost_equal(clf.dual_coef_, clf2.dual_coef_)
    assert clf.gamma == clf2.gamma
    assert_almost_equal(np.mean(pred == iris.target), .8, decimal=2)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    K = np.zeros_like(K)
    for i in range(len(iris.data)):
        for j in range(clf.h_):
            K[i, j] = np.dot(iris.data[i], iris.data[j])

    pred = clf.predict(K)
    assert_almost_equal(np.mean(pred == iris.target), .8, decimal=2)

    clf = KRC(kernel=kfunc)
    clf.fit(iris.data, iris.target)
    assert_almost_equal(np.mean(pred == iris.target), .8, decimal=2)


def test_tweak_params():
    # Make sure some tweaking of parameters works.
    clf = KRC(kernel='linear', gamma=0.0, alpha=0.5)
    clf.fit(X, Y)
    # dual_coef = np.array([[0.47826087, 0.52173913],
    #                       [0.65217391, 0.34782609],
    #                       [0.47826087, 0.52173913],
    #                       [0.34782609, 0.65217391],
    #                       [0.52173913, 0.47826087],
    #                       [0.52173913, 0.47826087]])
    dual_coef = np.array([[0.93333333, 1.06666667],
                          [1.28888889, 0.71111111],
                          [0.93333333, 1.06666667],
                          [0.71111111, 1.28888889],
                          [1.06666667, 0.93333333],
                          [1.06666667, 0.93333333]])
    assert_almost_equal(clf.dual_coef_, dual_coef, decimal=8)
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    dual_coef = np.array([[0.47826087, 0.52173913],
                          [0.52173913, 0.47826087],
                          [0.34782609, 0.65217391],
                          [0.47826087, 0.52173913],
                          [0.65217391, 0.34782609],
                          [0.52173913, 0.47826087]])
    clf.dual_coef_ = dual_coef
    assert_almost_equal(clf.predict([[-.1, -.1]]), [2])


def test_bad_input():
    # Test that it gives proper exception on deficient input
    # impossible value of C

    assert_raises(ValueError, KRC(alpha=-0.5).fit, X, Y)

    clf = KRC()
    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)

    # error for precomputed kernels
    clf = KRC(kernel='precomputed')
    assert_raises(ValueError, clf.fit, X, Y)

    # predict with sparse input when trained with dense
    clf = KRC().fit(X, Y)
    assert_raises(ValueError, clf.predict, sparse.lil_matrix(X))

    Xt = np.array(X).T
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = KRC()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)


def test_unicode_kernel():
    # Test that a unicode kernel name does not cause a TypeError
    if six.PY2:
        # Test unicode (same as str on python3)
        clf = KRC(kernel=u'linear')
        clf.fit(X, Y)
        clf.predict(T)

    # Test default behavior on both versions
    clf = KRC(kernel='linear')
    clf.fit(X, Y)
    clf.predict(T)


def test_linear_elm():
    clf = KRC(kernel='linear').fit(X, Y)

    # by default should have dual_coef
    assert_true(clf.dual_coef_.all())

    assert_array_equal(clf.predict(T), true_result)


def test_linear_elm_iris():
    # Test that Kernel ELM with linear kernel
    # gives plausible predictions on the iris dataset
    # Also, test symbolic class names (classes_).
    target = iris.target_names[iris.target]
    clf = KRC(kernel='linear').fit(iris.data, target)
    assert_equal(set(clf.classes_), set(iris.target_names))
    assert_greater(np.mean(clf.predict(iris.data) == target), 0.79)


def test_krc_clone_with_callable_kernel():
    # create KRC with callable linear kernel, check that results are the same
    # as with built-in linear kernel
    krc_callable = KRC(kernel=lambda x, y: np.dot(x, y.T))
    # clone for checking clonability with lambda functions..
    krc_cloned = base.clone(krc_callable)
    krc_cloned.fit(iris.data, iris.target)

    krc_builtin = KRC(kernel='linear')
    krc_builtin.fit(iris.data, iris.target)

    assert_array_almost_equal(krc_cloned.dual_coef_,
                              krc_builtin.dual_coef_)
    assert_array_equal(krc_cloned.predict(iris.data),
                       krc_builtin.predict(iris.data))
    assert_array_almost_equal(krc_cloned.predict(iris.data),
                              krc_builtin.predict(iris.data))


def test_krc_bad_kernel():
    krc = KRC(kernel=lambda x, y: x)
    assert_raises(ValueError, krc.fit, X, Y)


def test_unfitted():
    X = "foo!"  # input validation not required when KRC is not fitted

    clf = KRC()
    assert_raises_regexp(Exception, r"not\b.*\bfitted\b",
                         clf.predict, X)
