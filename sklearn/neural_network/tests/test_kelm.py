"""
Testing for Kernel Extreme Learning Machine module (sklearn.neural_network.kelm)

"""
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from numpy.testing import assert_almost_equal
from scipy import sparse
from sklearn import datasets, base
from sklearn.neural_network import KernelELM
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_equal, assert_true
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_raises_regexp
from sklearn.utils.testing import assert_raises
from sklearn.externals import six
from sklearn.utils.estimator_checks import check_estimator

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
    check_estimator(KernelELM)


def test_kelm_parameters():
    # Test parameters.
    clf = KernelELM(kernel='linear').fit(X, Y)
    assert_array_equal(clf.predict(X), Y)


def test_iris():
    # Check consistency on dataset iris.

    # shuffle the dataset so that labels are not ordered
    for k in ('linear', 'rbf'):
        clf = KernelELM(kernel=k).fit(iris.data, iris.target)
        assert_greater(np.mean(clf.predict(iris.data) == iris.target), 0.75)

    assert_array_equal(clf.classes_, np.sort(clf.classes_))


def test_precomputed():
    # KELM with a precomputed kernel.
    # We test it with a toy dataset and with iris.
    clf = KernelELM(kernel='precomputed')
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

    kfunc = lambda x, y: np.dot(x, y.T)
    clf = KernelELM(kernel=kfunc)
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_equal(pred, true_result)

    # test a precomputed kernel with the iris dataset
    # and check parameters against a linear KELM
    clf = KernelELM(kernel='precomputed')
    clf2 = KernelELM(kernel='linear', gamma=0.0)
    K = np.dot(iris.data, iris.data.T)
    clf.fit(K, iris.target)
    clf2.fit(iris.data, iris.target)
    pred = clf.predict(K)
    assert_array_almost_equal(clf.beta_, clf2.beta_)
    assert_array_almost_equal(clf.gamma_, clf2.gamma_)
    assert_almost_equal(np.mean(pred == iris.target), .8, decimal=2)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    K = np.zeros_like(K)
    for i in range(len(iris.data)):
        for j in range(clf.h_):
            K[i, j] = np.dot(iris.data[i], iris.data[j])

    pred = clf.predict(K)
    assert_almost_equal(np.mean(pred == iris.target), .8, decimal=2)

    clf = KernelELM(kernel=kfunc)
    clf.fit(iris.data, iris.target)
    assert_almost_equal(np.mean(pred == iris.target), .8, decimal=2)


def test_tweak_params():
    # Make sure some tweaking of parameters works.
    clf = KernelELM(kernel='linear', C=1.0)
    clf.fit(X, Y)
    beta = np.array([[0.47826087, 0.52173913],
                     [0.65217391, 0.34782609],
                     [0.47826087, 0.52173913],
                     [0.34782609, 0.65217391],
                     [0.52173913, 0.47826087],
                     [0.52173913, 0.47826087]])
    assert_almost_equal(clf.beta_, beta, decimal=8)
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    beta = np.array([[0.47826087, 0.52173913],
                     [0.52173913, 0.47826087],
                     [0.34782609, 0.65217391],
                     [0.47826087, 0.52173913],
                     [0.65217391, 0.34782609],
                     [0.52173913, 0.47826087]])
    clf.beta_ = beta
    assert_almost_equal(clf.predict([[-.1, -.1]]), [2])


def test_bad_input():
    # Test that it gives proper exception on deficient input
    # impossible value of C

    assert_raises(ValueError, KernelELM(C=-1).fit, X, Y)

    clf = KernelELM()
    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)

    # error for precomputed kernels
    clf = KernelELM(kernel='precomputed')
    assert_raises(ValueError, clf.fit, X, Y)

    # predict with sparse input when trained with dense
    clf = KernelELM().fit(X, Y)
    assert_raises(ValueError, clf.predict, sparse.lil_matrix(X))

    Xt = np.array(X).T
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = KernelELM()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)


def test_unicode_kernel():
    # Test that a unicode kernel name does not cause a TypeError
    if six.PY2:
        # Test unicode (same as str on python3)
        clf = KernelELM(kernel=u'linear')
        clf.fit(X, Y)
        clf.predict(T)

    # Test default behavior on both versions
    clf = KernelELM(kernel='linear')
    clf.fit(X, Y)
    clf.predict(T)


def test_linear_elm():
    clf = KernelELM(kernel='linear').fit(X, Y)

    # by default should have beta
    assert_true(clf.beta_.all())

    assert_array_equal(clf.predict(T), true_result)


def test_linear_elm_iris():
    # Test that Kernel ELM with linear kernel
    # gives plausible predictions on the iris dataset
    # Also, test symbolic class names (classes_).
    target = iris.target_names[iris.target]
    clf = KernelELM(kernel='linear').fit(iris.data, target)
    assert_equal(set(clf.classes_), set(iris.target_names))
    assert_greater(np.mean(clf.predict(iris.data) == target), 0.79)


def test_kelm_clone_with_callable_kernel():
    # create KELM with callable linear kernel, check that results are the same
    # as with built-in linear kernel
    kelm_callable = KernelELM(kernel=lambda x, y: np.dot(x, y.T))
    # clone for checking clonability with lambda functions..
    kelm_cloned = base.clone(kelm_callable)
    kelm_cloned.fit(iris.data, iris.target)

    kelm_builtin = KernelELM(kernel='linear')
    kelm_builtin.fit(iris.data, iris.target)

    assert_array_almost_equal(kelm_cloned.beta_,
                              kelm_builtin.beta_)
    assert_array_equal(kelm_cloned.predict(iris.data),
                       kelm_builtin.predict(iris.data))
    assert_array_almost_equal(kelm_cloned.predict(iris.data),
                              kelm_builtin.predict(iris.data))


def test_kelm_bad_kernel():
    kelm = KernelELM(kernel=lambda x, y: x)
    assert_raises(ValueError, kelm.fit, X, Y)


def test_unfitted():
    X = "foo!"  # input validation not required when KELM is not fitted

    clf = KernelELM()
    assert_raises_regexp(Exception, r"not\b.*\bfitted\b",
                         clf.predict, X)


if __name__ == '__main__':
    test_linear_elm_iris()
