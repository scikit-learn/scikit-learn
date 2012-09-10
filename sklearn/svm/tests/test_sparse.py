import numpy as np
from scipy import sparse
from sklearn import datasets, svm, linear_model, base
from numpy.testing import assert_array_almost_equal, \
     assert_array_equal, assert_equal

from nose.tools import assert_raises, assert_true
from sklearn.datasets.samples_generator import make_classification
from sklearn.svm.tests import test_svm

# test sample 1
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
X_sp = sparse.lil_matrix(X)
Y = [1, 1, 1, 2, 2, 2]
T = np.array([[-1, -1], [2, 2], [3, 2]])
true_result = [1, 2, 2]

# test sample 2
X2 = np.array([[0, 0, 0], [1, 1, 1], [2, 0, 0, ],
               [0, 0, 2], [3, 3, 3]])
X2_sp = sparse.dok_matrix(X2)
Y2 = [1, 2, 2, 2, 3]
T2 = np.array([[-1, -1, -1], [1, 1, 1], [2, 2, 2]])
true_result2 = [1, 2, 3]


iris = datasets.load_iris()
# permute
rng = np.random.RandomState(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]
# sparsify
iris.data = sparse.csr_matrix(iris.data)


def test_svc():
    """Check that sparse SVC gives the same result as SVC"""

    clf = svm.SVC(kernel='linear').fit(X, Y)
    sp_clf = svm.SVC(kernel='linear').fit(X_sp, Y)

    assert_array_equal(sp_clf.predict(T), true_result)

    assert_true(sparse.issparse(sp_clf.support_vectors_))
    assert_array_almost_equal(clf.support_vectors_,
            sp_clf.support_vectors_.todense())

    assert_true(sparse.issparse(sp_clf.dual_coef_))
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())

    assert_true(sparse.issparse(sp_clf.coef_))
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T), sp_clf.predict(T))

    # refit with a different dataset
    clf.fit(X2, Y2)
    sp_clf.fit(X2_sp, Y2)
    assert_array_almost_equal(clf.support_vectors_,
            sp_clf.support_vectors_.todense())
    assert_array_almost_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())
    assert_array_almost_equal(clf.coef_, sp_clf.coef_.todense())
    assert_array_almost_equal(clf.predict(T2), sp_clf.predict(T2))


def test_svc_iris():
    """Test the sparse SVC with the iris dataset"""
    for k in ('linear', 'poly', 'rbf'):
        sp_clf = svm.SVC(kernel=k).fit(iris.data, iris.target)
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
    clf = svm.NuSVC(nu=0.0)
    assert_raises(ValueError, clf.fit, X_sp, Y)

    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X_sp, Y2)

    clf = svm.SVC()
    clf.fit(X_sp, Y)
    assert_array_equal(clf.predict(T), true_result)


def test_linearsvc():
    """
    Similar to test_SVC
    """
    clf = svm.LinearSVC().fit(X, Y)
    sp_clf = svm.LinearSVC().fit(X_sp, Y)

    assert_true(sp_clf.fit_intercept)

    assert_array_almost_equal(clf.raw_coef_, sp_clf.raw_coef_, decimal=4)

    assert_array_almost_equal(clf.predict(X), sp_clf.predict(X_sp))

    clf.fit(X2, Y2)
    sp_clf.fit(X2_sp, Y2)

    assert_array_almost_equal(clf.raw_coef_, sp_clf.raw_coef_, decimal=4)


def test_linearsvc_iris():
    """Test the sparse LinearSVC with the iris dataset"""

    sp_clf = svm.LinearSVC().fit(iris.data, iris.target)
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
    for clf in (linear_model.LogisticRegression(),
                svm.LinearSVC(),
                svm.SVC()):
        clf.set_params(class_weight={0: 5})
        clf.fit(X_[:180], y_[:180])
        y_pred = clf.predict(X_[180:])
        assert_true(np.sum(y_pred == y_[180:]) >= 11)


def test_sample_weights():
    """
    Test weights on individual samples
    """
    clf = svm.SVC()
    clf.fit(X_sp, Y)
    assert_array_equal(clf.predict(X[2]), [1.])

    sample_weight = [.1] * 3 + [10] * 3
    clf.fit(X_sp, Y, sample_weight=sample_weight)
    assert_array_equal(clf.predict(X[2]), [2.])


def test_sparse_liblinear_intercept_handling():
    """
    Test that sparse liblinear honours intercept_scaling param
    """
    test_svm.test_dense_liblinear_intercept_handling(svm.LinearSVC)


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
    sp_clf = svm.SVC(kernel='linear').fit(sparse.coo_matrix(X), y)

    assert_array_equal(clf.support_vectors_, sp_clf.support_vectors_.todense())
    assert_array_equal(clf.dual_coef_, sp_clf.dual_coef_.todense())


def test_sparse_svc_clone_with_callable_kernel():
    # first, test that we raise a value error for "sparse kernels"
    # this test is only relevant for the deprecated sparse.SVC class.
    sp = svm.sparse.SVC(C=1, kernel=lambda x, y: x * y.T, probability=True)
    assert_raises(ValueError, sp.fit, X_sp, Y)

    # Test that the "dense_fit" is called even though we use sparse input
    # meaning that everything works fine.
    a = svm.SVC(C=1, kernel=lambda x, y: x * y.T, probability=True)
    b = base.clone(a)

    b.fit(X_sp, Y)
    pred = b.predict(X_sp)
    b.predict_proba(X_sp)

    dense_svm = svm.SVC(C=1, kernel=lambda x, y: np.dot(x, y.T),
            probability=True)
    pred_dense = dense_svm.fit(X, Y).predict(X)
    assert_array_equal(pred_dense, pred)
    # b.decision_function(X_sp)  # XXX : should be supported


if __name__ == '__main__':
    import nose
    nose.runmodule()
