import numpy as np
from scipy import linalg
from scipy import sparse
from sklearn import datasets, svm, linear_model
from numpy.testing import assert_array_almost_equal, \
     assert_array_equal, assert_equal

from nose.tools import assert_raises, assert_true
from sklearn.datasets.samples_generator import make_classification
from sklearn.svm.tests import test_svm

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
iris.data = sparse.csr_matrix(iris.data)


def test_SVC():
    """Check that sparse SVC gives the same result as SVC"""

    clf = svm.SVC(kernel='linear').fit(X, Y)
    sp_clf = svm.sparse.SVC(kernel='linear').fit(X, Y)

    assert_array_equal(sp_clf.predict(T), true_result)

    assert sparse.issparse(sp_clf.support_vectors_)
    assert_array_almost_equal(clf.support_vectors_,
            sp_clf.support_vectors_.todense())

    assert sparse.issparse(sp_clf.dual_coef_)
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())

    assert sparse.issparse(sp_clf.coef_)
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T), sp_clf.predict(T))

    # refit with a different dataset
    clf.fit(X2, Y2)
    sp_clf.fit(X2, Y2)
    assert_array_almost_equal(clf.support_vectors_,
            sp_clf.support_vectors_.todense())
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T2), sp_clf.predict(T2))


def test_SVC_iris():
    """Test the sparse SVC with the iris dataset"""
    for k in ('linear', 'poly', 'rbf'):
        sp_clf = svm.sparse.SVC(kernel=k).fit(iris.data, iris.target)
        clf = svm.SVC(kernel=k).fit(iris.data.todense(), iris.target)

        assert_array_almost_equal(clf.support_vectors_,
                sp_clf.support_vectors_.todense())
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

    X_, y_ = make_classification(n_samples=200, n_features=100,
                                 weights=[0.833, 0.167], random_state=0)

    X_ = sparse.csr_matrix(X_)
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


def test_sparse_realdata():
    """
    Test on a subset from the 20newsgroups dataset.

    This catchs some bugs if input is not correctly converted into
    sparse format or weights are not correctly initialized.
    """

    data = np.array([0.03771744,  0.1003567,  0.01174647,  0.027069])
    indices = np.array([6, 5, 35, 31])
    indptr = np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 4])
    X = sparse.csr_matrix((data, indices, indptr))
    y = np.array(
        [1.,  0.,  2.,  2.,  1.,  1.,  1.,  2.,  2.,  0.,  1.,  2.,  2.,
        0.,  2.,  0.,  3.,  0.,  3.,  0.,  1.,  1.,  3.,  2.,  3.,  2.,
        0.,  3.,  1.,  0.,  2.,  1.,  2.,  0.,  1.,  0.,  2.,  3.,  1.,
        3.,  0.,  1.,  0.,  0.,  2.,  0.,  1.,  2.,  2.,  2.,  3.,  2.,
        0.,  3.,  2.,  1.,  2.,  3.,  2.,  2.,  0.,  1.,  0.,  1.,  2.,
        3.,  0.,  0.,  2.,  2.,  1.,  3.,  1.,  1.,  0.,  1.,  2.,  1.,
        1.,  3.])

    clf = svm.SVC(kernel='linear').fit(X.todense(), y)
    sp_clf = svm.sparse.SVC(kernel='linear').fit(X, y)

    assert_array_equal(clf.support_vectors_, sp_clf.support_vectors_.todense())
    assert_array_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())


def test_sparse_scale_C():
    """Check that sparse LibSVM/LibLinear works ok with scaling of C"""

    params = dict(kernel='linear', C=0.1)
    klasses = [(svm.SVC, svm.sparse.SVC, params),
               (svm.SVR, svm.sparse.SVR, params),
               (svm.NuSVR, svm.sparse.NuSVR, params),
               (svm.LinearSVC, svm.sparse.LinearSVC, {}),
               (linear_model.LogisticRegression,
                    linear_model.sparse.LogisticRegression, {})
              ]

    for klass, sparse_klass, this_params in klasses:
        clf = klass(scale_C=True, **this_params).fit(X, Y)
        clf_no_scale = klass(scale_C=False, **this_params).fit(X, Y)
        sp_clf = sparse_klass(scale_C=True, **this_params).fit(X, Y)

        sp_clf_coef_ = sp_clf.coef_
        if sparse.issparse(sp_clf_coef_):
            sp_clf_coef_ = sp_clf_coef_.todense()
        assert_array_almost_equal(clf.coef_, sp_clf_coef_, 5)

        error_with_scale = linalg.norm(clf_no_scale.coef_
                           - sp_clf_coef_) / linalg.norm(clf_no_scale.coef_)
        assert_true(error_with_scale > 1e-3)


if __name__ == '__main__':
    import nose
    nose.runmodule()
