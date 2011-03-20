import numpy as np
import scipy.sparse
from scikits.learn import datasets, svm, linear_model
from numpy.testing import assert_array_almost_equal, \
     assert_array_equal, assert_equal

from nose.tools import assert_raises
from scikits.learn.datasets.samples_generator import test_dataset_classif
from . import test_svm

# test sample 1
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
Y = [1, 1, 1, 2, 2, 2]
T = np.array([[-1, -1], [2, 2], [3, 2]])
true_result = [1, 2, 2]

# test sample 2
X2 = np.array([[0, 0, 0], [1, 1, 1], [2, 0, 0, ],
               [0, 0, 2], [3, 3, 3]])
Y2 = [1, 2, 2, 2, 3]
T2 = np.array([[-1, -1, -1], [1, 1, 1], [2, 2, 2]])
true_result2 = [1, 2, 3]


iris = datasets.load_iris()
# permute
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
# sparsify
iris.data = scipy.sparse.csr_matrix(iris.data)


def test_SVC():
    """Check that sparse SVC gives the same result as SVC"""

    clf = svm.SVC(kernel='linear').fit(X, Y)
    sp_clf = svm.sparse.SVC(kernel='linear').fit(X, Y)

    assert_array_equal(sp_clf.predict(T), true_result)

    assert scipy.sparse.issparse(sp_clf.support_vectors_)
    assert_array_almost_equal(clf.support_vectors_, sp_clf.support_vectors_.todense())

    assert scipy.sparse.issparse(sp_clf.dual_coef_)
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())

    assert scipy.sparse.issparse(sp_clf.coef_)
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T), sp_clf.predict(T))

    # refit with a different dataset
    clf.fit(X2, Y2)
    sp_clf.fit(X2, Y2)
    assert_array_almost_equal(clf.support_vectors_, sp_clf.support_vectors_.todense())
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T2), sp_clf.predict(T2))


def test_SVC_iris():
    """Test the sparse SVC with the iris dataset"""
    for k in ('linear', 'poly', 'rbf'):
        sp_clf = svm.sparse.SVC(kernel=k).fit(iris.data, iris.target)
        clf = svm.SVC(kernel=k).fit(iris.data.todense(), iris.target)

        assert_array_almost_equal(clf.support_vectors_, sp_clf.support_vectors_.todense())
        assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
        assert_array_almost_equal(
            clf.predict(iris.data.todense()), sp_clf.predict(iris.data))
        if k == 'linear':
            assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())


def test_error():
    """
    Test that it gives proper exception on deficient input
    """
    # impossible value of C
    assert_raises(ValueError, svm.SVC(C=-1).fit, X, Y)

    # impossible value of nu
    clf = svm.sparse.NuSVC(nu=0.0)
    assert_raises(ValueError, clf.fit, X, Y)

    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)
    assert_raises(ValueError, svm.SVC, X, Y2)

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

    assert_array_almost_equal(clf.raw_coef_, sp_clf.raw_coef_, decimal=4)

    assert_array_almost_equal(clf.predict(X), sp_clf.predict(X))

    clf.fit(X2, Y2)
    sp_clf.fit(X2, Y2)

    assert_array_almost_equal(clf.raw_coef_, sp_clf.raw_coef_, decimal=4)


def test_LinearSVC_iris():
    """Test the sparse LinearSVC with the iris dataset"""

    sp_clf = svm.sparse.LinearSVC().fit(iris.data, iris.target)
    clf = svm.LinearSVC().fit(iris.data.todense(), iris.target)

    assert_array_almost_equal(clf.label_, sp_clf.label_)
    assert_equal(clf.fit_intercept, sp_clf.fit_intercept)

    assert_array_almost_equal(clf.raw_coef_, sp_clf.raw_coef_, decimal=1)
    assert_array_almost_equal(
        clf.predict(iris.data.todense()), sp_clf.predict(iris.data))

    # check decision_function
    pred = np.argmax(sp_clf.decision_function(iris.data), 1)
    assert_array_almost_equal(pred, clf.predict(iris.data.todense()))


def test_weight():
    """
    Test class weights
    """

    X_, y_ = test_dataset_classif(n_samples=200, n_features=100, param=[5, 1],
                                  seed=0)
    X_ = scipy.sparse.csr_matrix(X_)
    for clf in (linear_model.sparse.LogisticRegression(),
                svm.sparse.LinearSVC(),
                svm.sparse.SVC()):
        clf.fit(X_[:180], y_[:180], class_weight={0: 5})
        y_pred = clf.predict(X_[180:])
        assert np.sum(y_pred == y_[180:]) >= 11


def test_sample_weights():
    """
    Test weights on individual samples
    """
    clf = svm.sparse.SVC()
    clf.fit(X, Y)
    assert_array_equal(clf.predict(X[2]), [1.])

    sample_weight = [.1] * 3 + [10] * 3
    clf.fit(X, Y, sample_weight=sample_weight)
    assert_array_equal(clf.predict(X[2]), [2.])


def test_sparse_liblinear_intercept_handling():
    """
    Test that sparse liblinear honours intercept_scaling param
    """
    test_svm.test_dense_liblinear_intercept_handling(svm.sparse.LinearSVC)


if __name__ == '__main__':
    import nose
    nose.runmodule()
