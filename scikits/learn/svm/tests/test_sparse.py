import numpy as np
import scipy.sparse
from scikits.learn import datasets, svm
from numpy.testing import assert_array_almost_equal, \
     assert_array_equal, assert_equal, assert_raises

# test sample 1
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

# test sample 2
X2 = [[0, 0, 0], [1, 1, 1], [2, 0, 0, ],
      [0, 0, 2], [3, 3, 3]]
Y2 = [1, 2, 2, 2, 3]
T2 = [[-1, -1, -1], [1, 1, 1], [2, 2, 2]]
true_result2 = [1, 2, 3]


iris = datasets.load_iris()


def test_SVC():
    """Check that sparse SVC gives the same result as SVC"""

    clf = svm.SVC(kernel='linear').fit(X, Y)
    sp_clf = svm.sparse.SVC(kernel='linear').fit(X, Y)

    assert_array_equal(sp_clf.predict(T), true_result)


    assert scipy.sparse.issparse(sp_clf.support_)
    assert_array_almost_equal(clf.support_, sp_clf.support_.todense())

    assert scipy.sparse.issparse (sp_clf.dual_coef_)
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())

    assert scipy.sparse.issparse (sp_clf.coef_)
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T), sp_clf.predict(T))

    # refit with a different dataset
    clf.fit(X2, Y2)
    sp_clf.fit(X2, Y2)
    assert_array_almost_equal(clf.support_, sp_clf.support_.todense())
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T2), sp_clf.predict(T2))


def test_SVC_iris():
    """Test the sparse SVC with the iris dataset"""
    for k in ('linear', 'rbf'):
        sp_clf = svm.sparse.SVC(kernel=k).fit(iris.data, iris.target)
        clf = svm.SVC(kernel=k).fit(iris.data, iris.target)

        assert_array_almost_equal(clf.support_, sp_clf.support_.todense())
        assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
        assert_array_almost_equal(clf.predict(iris.data), sp_clf.predict(iris.data))
        if k == 'linear':
            assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())



def test_error():
    """
    Test that it gives proper exception on deficient input
    """
    # impossible value of C
    assert_raises (ValueError, svm.SVC(C=-1).fit, X, Y)

    # impossible value of nu
    clf = svm.sparse.NuSVC(nu=0.0)
    assert_raises(ValueError, clf.fit, X, Y)

    Y2 = Y[:-1] # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)
    assert_raises(AssertionError, svm.SVC, X, Y2)

    clf = svm.sparse.SVC()
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)


def test_LinearSVC():
    """
    Similar to test_SVC
    """
    clf = svm.LinearSVC().fit(X, Y)
    sp_clf = svm.sparse.LinearSVC().fit(X, Y)

    assert sp_clf.fit_intercept
    
    assert_array_almost_equal (clf.raw_coef_, sp_clf.raw_coef_, decimal=4)

    assert_array_almost_equal (clf.predict(X), sp_clf.predict(X))

    clf.fit(X2, Y2)
    sp_clf.fit(X2, Y2)

    assert_array_almost_equal (clf.raw_coef_, sp_clf.raw_coef_, decimal=4)


def test_LinearSVC_iris():
    """Test the sparse LinearSVC with the iris dataset"""
    iris = datasets.load_iris()
    sp_clf = svm.sparse.LinearSVC().fit(iris.data, iris.target)
    clf = svm.LinearSVC().fit(iris.data, iris.target)

    assert_array_almost_equal(clf.label_, sp_clf.label_)
    assert_equal (clf.fit_intercept, sp_clf.fit_intercept)

    assert_array_almost_equal(clf.raw_coef_, sp_clf.raw_coef_, decimal=1)
    assert_array_almost_equal(clf.predict(iris.data), sp_clf.predict(iris.data))
