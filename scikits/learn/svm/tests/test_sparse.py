import scipy.sparse
from scikits.learn import datasets, svm
from numpy.testing import assert_array_almost_equal, \
                          assert_equal

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


def test_SVC():
    """Check that sparse SVC gives the same result as SVC"""

    clf = svm.SVC(kernel='linear').fit(X, Y)
    sp_clf = svm.sparse.SVC(kernel='linear').fit(X, Y)

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
    iris = datasets.load_iris()
    sp_clf = svm.sparse.SVC(kernel='linear').fit(iris.data, iris.target)
    clf = svm.SVC(kernel='linear').fit(iris.data, iris.target)

    assert_array_almost_equal(clf.support_, sp_clf.support_.todense())
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(iris.data), sp_clf.predict(iris.data))


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
