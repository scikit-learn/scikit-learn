"""
Testing for Support Vector Machine module (sklearn.svm)

TODO: remove hard coded numerical results when possible
"""

import numpy as np
from scipy import linalg
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
                          assert_almost_equal
from nose.tools import assert_raises, assert_true

from sklearn import svm, linear_model, datasets, metrics
from sklearn.datasets.samples_generator import make_classification
from sklearn.utils import check_random_state

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

# also load the iris dataset
iris = datasets.load_iris()
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_libsvm_parameters():
    """
    Test parameters on classes that make use of libsvm.
    """
    clf = svm.SVC(kernel='linear', C=len(X)).fit(X, Y)
    assert_array_equal(clf.dual_coef_, [[0.25, -.25]])
    assert_array_equal(clf.support_, [1, 3])
    assert_array_equal(clf.support_vectors_, (X[1], X[3]))
    assert_array_equal(clf.intercept_, [0.])
    assert_array_equal(clf.predict(X), Y)


def test_libsvm_iris():
    """
    Check consistency on dataset iris.
    """

    # shuffle the dataset so that labels are not ordered
    for k in ('linear', 'rbf'):
        clf = svm.SVC(kernel=k).fit(iris.data, iris.target)
        assert_true(np.mean(clf.predict(iris.data) == iris.target) > 0.9)

    assert_array_equal(clf.label_, np.sort(clf.label_))

    # check also the low-level API
    model = svm.libsvm.fit(iris.data, iris.target.astype(np.float64))
    pred = svm.libsvm.predict(iris.data, *model)
    assert_true(np.mean(pred == iris.target) > .95)

    model = svm.libsvm.fit(iris.data,
            iris.target.astype(np.float64), kernel='linear')
    pred = svm.libsvm.predict(iris.data, *model, kernel='linear')
    assert_true(np.mean(pred == iris.target) > .95)

    pred = svm.libsvm.cross_validation(iris.data,
            iris.target.astype(np.float64), 5, kernel='linear')
    assert_true(np.mean(pred == iris.target) > .95)


def test_single_sample_1d():
    """
    Test whether SVCs work on a single sample given as a 1-d array
    """

    clf = svm.SVC().fit(X, Y)
    clf.predict(X[0])

    clf = svm.LinearSVC().fit(X, Y)
    clf.predict(X[0])


def test_precomputed():
    """
    SVC with a precomputed kernel.

    We test it with a toy dataset and with iris.
    """
    clf = svm.SVC(kernel='precomputed', C=len(X))
    # Gram matrix for train data (square matrix)
    # (we use just a linear kernel)
    K = np.dot(X, np.array(X).T)
    clf.fit(K, Y)
    # Gram matrix for test data (rectangular matrix)
    KT = np.dot(T, np.array(X).T)
    pred = clf.predict(KT)

    assert_array_equal(clf.dual_coef_, [[0.25, -.25]])
    assert_array_equal(clf.support_, [1, 3])
    assert_array_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.support_, [1, 3])
    assert_array_equal(pred, true_result)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    KT = np.zeros_like(KT)
    for i in range(len(T)):
        for j in clf.support_:
            KT[i, j] = np.dot(T[i], X[j])

    pred = clf.predict(KT)
    assert_array_equal(pred, true_result)

    # same as before, but using a callable function instead of the kernel
    # matrix. kernel is just a linear kernel

    kfunc = lambda x, y: np.dot(x, y.T)
    clf = svm.SVC(kernel=kfunc, C=len(X))
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_equal(clf.dual_coef_, [[0.25, -.25]])
    assert_array_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.support_, [1, 3])
    assert_array_equal(pred, true_result)

    # test a precomputed kernel with the iris dataset
    # and check parameters against a linear SVC
    clf = svm.SVC(kernel='precomputed', C=len(iris.data))
    clf2 = svm.SVC(kernel='linear', C=len(iris.data))
    K = np.dot(iris.data, iris.data.T)
    clf.fit(K, iris.target)
    clf2.fit(iris.data, iris.target)
    pred = clf.predict(K)
    assert_array_almost_equal(clf.support_, clf2.support_)
    assert_array_almost_equal(clf.dual_coef_, clf2.dual_coef_)
    assert_array_almost_equal(clf.intercept_, clf2.intercept_)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)

    # Gram matrix for test data but compute KT[i,j]
    # for support vectors j only.
    K = np.zeros_like(K)
    for i in range(len(iris.data)):
        for j in clf.support_:
            K[i, j] = np.dot(iris.data[i], iris.data[j])

    pred = clf.predict(K)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)

    clf = svm.SVC(kernel=kfunc, C=len(X))
    clf.fit(iris.data, iris.target)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)


def test_SVR():
    """
    Test Support Vector Regression
    """

    diabetes = datasets.load_diabetes()
    n_samples = len(diabetes.data)
    for clf in (svm.NuSVR(kernel='linear', nu=.4, C=n_samples),
                svm.NuSVR(kernel='linear', nu=.4, C=10. * n_samples),
                svm.SVR(kernel='linear', C=10. * n_samples),
                svm.sparse.NuSVR(kernel='linear', nu=.4, C=n_samples),
                svm.sparse.NuSVR(kernel='linear', nu=.4, C=10. * n_samples),
                svm.sparse.SVR(kernel='linear', C=10. * n_samples)):
        clf.fit(diabetes.data, diabetes.target)
        assert_true(clf.score(diabetes.data, diabetes.target) > 0.02)


def test_oneclass():
    """
    Test OneClassSVM
    """
    clf = svm.OneClassSVM()
    clf.fit(X)
    pred = clf.predict(T)

    assert_array_almost_equal(pred, [-1, -1, -1])
    assert_array_almost_equal(clf.intercept_, [-1.008], decimal=3)
    assert_array_almost_equal(clf.dual_coef_,
                              [[0.632, 0.233, 0.633, 0.234, 0.632, 0.633]],
                              decimal=3)
    assert_raises(ValueError, lambda: clf.coef_)


def test_oneclass_decision_function():
    """
    Test OneClassSVM decision function
    """
    clf = svm.OneClassSVM()
    rnd = check_random_state(0)

    # Generate train data
    X = 0.3 * rnd.randn(100, 2)
    X_train = np.r_[X + 2, X - 2]

    # Generate some regular novel observations
    X = 0.3 * rnd.randn(20, 2)
    X_test = np.r_[X + 2, X - 2]
    # Generate some abnormal novel observations
    X_outliers = np.random.uniform(low=-4, high=4, size=(20, 2))

    # fit the model
    clf = svm.OneClassSVM(nu=0.1, kernel="rbf", gamma=0.1)
    clf.fit(X_train)

    # predict things
    y_pred_test = clf.predict(X_test)
    assert_true(np.mean(y_pred_test == 1) > .9)
    y_pred_outliers = clf.predict(X_outliers)
    assert_true(np.mean(y_pred_outliers == -1) > .9)
    dec_func_test = clf.decision_function(X_test)
    assert_array_equal((dec_func_test > 0).ravel(), y_pred_test == 1)
    dec_func_outliers = clf.decision_function(X_outliers)
    assert_array_equal((dec_func_outliers > 0).ravel(),  y_pred_outliers == 1)


def test_tweak_params():
    """
    Make sure some tweaking of parameters works.

    We change clf.dual_coef_ at run time and expect .predict() to change
    accordingly. Notice that this is not trivial since it involves a lot
    of C/Python copying in the libsvm bindings.

    The success of this test ensures that the mapping between libsvm and
    the python classifier is complete.
    """
    clf = svm.SVC(kernel='linear', C=len(X))
    clf.fit(X, Y)
    assert_array_equal(clf.dual_coef_, [[.25, -.25]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    clf.dual_coef_ = np.array([[.0, 1.]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [2])


def test_probability():
    """
    Predict probabilities using SVC

    This uses cross validation, so we use a slightly bigger testing set.
    """

    for clf in (
        svm.SVC(probability=True, C=len(X)),
        svm.NuSVC(probability=True),
        svm.sparse.SVC(probability=True, C=len(X)),
        svm.sparse.NuSVC(probability=True)):

        clf.fit(iris.data, iris.target)

        prob_predict = clf.predict_proba(iris.data)
        assert_array_almost_equal(
            np.sum(prob_predict, 1), np.ones(iris.data.shape[0]))
        assert_true(np.mean(np.argmax(prob_predict, 1)
                       == clf.predict(iris.data)) > 0.9)

        assert_almost_equal(clf.predict_proba(iris.data),
                            np.exp(clf.predict_log_proba(iris.data)), 8)


def test_decision_function():
    """
    Test decision_function

    Sanity check, test that decision_function implemented in python
    returns the same as the one in libsvm

    """
    # multi class:
    clf = svm.SVC(kernel='linear').fit(iris.data, iris.target)

    dec = np.dot(iris.data, clf.coef_.T) + clf.intercept_

    assert_array_almost_equal(dec, clf.decision_function(iris.data))

    # binary:
    clf.fit(X, Y)
    dec = np.dot(X, clf.coef_.T) + clf.intercept_
    prediction = clf.predict(X)
    assert_array_almost_equal(dec, clf.decision_function(X))
    assert_array_almost_equal(prediction, clf.label_[(clf.decision_function(X)
        > 0).astype(np.int).ravel()])
    expected = np.array([[-1.], [-0.66], [-1.], [0.66], [1.], [1.]])
    assert_array_almost_equal(clf.decision_function(X), expected, 2)


def test_weight():
    """
    Test class weights
    """
    clf = svm.SVC(class_weight={1: 0.1}, C=len(X))
    # we give a small weights to class 1
    clf.fit(X, Y)
    # so all predicted values belong to class 2
    assert_array_almost_equal(clf.predict(X), [2] * 6)

    X_, y_ = make_classification(n_samples=200, n_features=100,
                                 weights=[0.833, 0.167], random_state=0)

    for clf in (linear_model.LogisticRegression(C=180),
            svm.LinearSVC(C=len(X)), svm.SVC(C=len(X))):
        clf.set_params(class_weight={0: 5})
        clf.fit(X_[: 180], y_[: 180])
        y_pred = clf.predict(X_[180:])
        assert_true(np.sum(y_pred == y_[180:]) >= 11)


def test_sample_weights():
    """
    Test weights on individual samples
    """
    # TODO: check on NuSVR, OneClass, etc.
    clf = svm.SVC()
    clf.fit(X, Y)
    assert_array_equal(clf.predict(X[2]), [1.])

    sample_weight = [.1] * 3 + [10] * 3
    clf.fit(X, Y, sample_weight=sample_weight)
    assert_array_equal(clf.predict(X[2]), [2.])


def test_auto_weight():
    """Test class weights for imbalanced data"""
    from sklearn.linear_model import LogisticRegression
    # we take as dataset a the two-dimensional projection of iris so
    # that it is not separable and remove half of predictors from
    # class 1
    from sklearn.svm.base import _get_class_weight
    X, y = iris.data[:, :2], iris.target
    unbalanced = np.delete(np.arange(y.size), np.where(y > 1)[0][::2])

    assert_true(np.argmax(_get_class_weight('auto', y[unbalanced])[0]) == 2)

    C = len(X)
    for clf in (svm.SVC(kernel='linear', C=C),
            svm.LinearSVC(C=C), LogisticRegression(C=C)):
        # check that score is better when class='auto' is set.
        y_pred = clf.fit(X[unbalanced], y[unbalanced]).predict(X)
        clf.set_params(class_weight='auto')
        y_pred_balanced = clf.fit(X[unbalanced], y[unbalanced],).predict(X)
        assert_true(metrics.f1_score(y, y_pred) <=
                metrics.f1_score(y, y_pred_balanced))


def test_bad_input():
    """
    Test that it gives proper exception on deficient input
    """
    # impossible value of C
    assert_raises(ValueError, svm.SVC(C=-1).fit, X, Y)

    # impossible value of nu
    clf = svm.NuSVC(nu=0.0)
    assert_raises(ValueError, clf.fit, X, Y)

    Y2 = Y[:-1]  # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)

    # Test with arrays that are non-contiguous.
    for clf in (svm.SVC(), svm.LinearSVC(), svm.sparse.SVC()):
        Xf = np.asfortranarray(X)
        assert_true(Xf.flags['C_CONTIGUOUS'] == False)
        yf = np.ascontiguousarray(np.tile(Y, (2, 1)).T)
        yf = yf[:, -1]
        assert_true(yf.flags['F_CONTIGUOUS'] == False)
        assert_true(yf.flags['C_CONTIGUOUS'] == False)
        clf.fit(Xf, yf)
        assert_array_equal(clf.predict(T), true_result)

    # error for precomputed kernelsx
    clf = svm.SVC(kernel='precomputed')
    assert_raises(ValueError, clf.fit, X, Y)

    Xt = np.array(X).T
    clf = svm.SVC(kernel='precomputed')
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = svm.SVC()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)


def test_LinearSVC_parameters():
    """
    Test possible parameter combinations in LinearSVC
    """
    # generate list of possible parameter combinations
    params = [(dual, loss, penalty) for dual in [True, False]
            for loss in ['l1', 'l2', 'lr'] for penalty in ['l1', 'l2']]

    for dual, loss, penalty in params:
            if loss == 'l1' and penalty == 'l1':
                assert_raises(ValueError, svm.LinearSVC, penalty=penalty,
                        loss=loss, dual=dual)
            elif loss == 'l1' and penalty == 'l2' and dual == False:
                assert_raises(ValueError, svm.LinearSVC, penalty=penalty,
                        loss=loss, dual=dual)
            elif penalty == 'l1' and dual == True:
                assert_raises(ValueError, svm.LinearSVC, penalty=penalty,
                        loss=loss, dual=dual)
            else:
                svm.LinearSVC(penalty=penalty, loss=loss, dual=dual)


def test_LinearSVC():
    """
    Test basic routines using LinearSVC
    """
    clf = svm.LinearSVC().fit(X, Y)

    # by default should have intercept
    assert_true(clf.fit_intercept)

    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.intercept_, [0], decimal=3)

    # the same with l1 penalty
    clf = svm.LinearSVC(penalty='l1', dual=False).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty with dual formulation
    clf = svm.LinearSVC(penalty='l2', dual=True).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty, l1 loss
    clf = svm.LinearSVC(penalty='l2', loss='l1', dual=True).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # test also decision function
    dec = clf.decision_function(T).ravel()
    res = (dec > 0).astype(np.int) + 1
    assert_array_equal(res, true_result)


def test_LinearSVC_iris():
    """
    Test that LinearSVC gives plausible predictions on the iris dataset
    """
    clf = svm.LinearSVC(C=len(iris.data)).fit(iris.data, iris.target)
    assert_true(np.mean(clf.predict(iris.data) == iris.target) > 0.8)

    dec = clf.decision_function(iris.data)
    pred = np.argmax(dec, 1)
    assert_array_equal(pred, clf.predict(iris.data))


def test_dense_liblinear_intercept_handling(classifier=svm.LinearSVC):
    """
    Test that dense liblinear honours intercept_scaling param
    """
    X = [[2, 1],
         [3, 1],
         [1, 3],
         [2, 3]]
    y = [0, 0, 1, 1]
    clf = classifier(fit_intercept=True, penalty='l1', loss='l2',
                     dual=False, C=4, tol=1e-7)
    assert_true(clf.intercept_scaling == 1, clf.intercept_scaling)
    assert_true(clf.fit_intercept)

    # when intercept_scaling is low the intercept value is highly "penalized"
    # by regularization
    clf.intercept_scaling = 1
    clf.fit(X, y)
    assert_almost_equal(clf.intercept_, 0, decimal=5)

    # when intercept_scaling is sufficiently high, the intercept value
    # is not affected by regularization
    clf.intercept_scaling = 100
    clf.fit(X, y)
    intercept1 = clf.intercept_
    assert_true(intercept1 < -1)

    # when intercept_scaling is sufficiently high, the intercept value
    # doesn't depend on intercept_scaling value
    clf.intercept_scaling = 1000
    clf.fit(X, y)
    intercept2 = clf.intercept_
    assert_array_almost_equal(intercept1, intercept2, decimal=2)


def test_liblinear_predict():
    """
    Test liblinear predict

    Sanity check, test that predict implemented in python
    returns the same as the one in libliblinear

    """
    # multi-class case
    clf = svm.LinearSVC().fit(iris.data, iris.target)
    weights = clf.coef_.T
    bias = clf.intercept_
    H = np.dot(iris.data, weights) + bias
    assert_array_equal(clf.predict(iris.data), H.argmax(axis=1))

    # binary-class case
    X = [[2, 1],
         [3, 1],
         [1, 3],
         [2, 3]]
    y = [0, 0, 1, 1]

    clf = svm.LinearSVC().fit(X, y)
    weights = np.ravel(clf.coef_)
    bias = clf.intercept_
    H = np.dot(X, weights) + bias
    assert_array_equal(clf.predict(X), (H > 0).astype(int))


def test_liblinear_set_coef():
    # multi-class case
    clf = svm.LinearSVC().fit(iris.data, iris.target)
    values = clf.decision_function(iris.data)
    clf.coef_ = clf.coef_.copy()
    clf.intercept_ = clf.intercept_.copy()
    values2 = clf.decision_function(iris.data)
    assert_array_equal(values, values2)

    # binary-class case
    X = [[2, 1],
         [3, 1],
         [1, 3],
         [2, 3]]
    y = [0, 0, 1, 1]

    clf = svm.LinearSVC().fit(X, y)
    values = clf.decision_function(X)
    clf.coef_ = clf.coef_.copy()
    clf.intercept_ = clf.intercept_.copy()
    values2 = clf.decision_function(X)
    assert_array_equal(values, values2)


def test_c_samples_scaling():
    """Test C scaling by n_samples
    """
    X = iris.data[iris.target != 2]
    y = iris.target[iris.target != 2]
    X2 = np.r_[X, X]
    y2 = np.r_[y, y]

    clfs = [svm.SVC(tol=1e-6, kernel='linear', C=0.1),
            svm.SVR(tol=1e-6, kernel='linear', C=100),
            svm.LinearSVC(tol=1e-6, C=0.1),
            linear_model.LogisticRegression(penalty='l1', tol=1e-6, C=100),
            linear_model.LogisticRegression(penalty='l2', tol=1e-6),
            svm.NuSVR(tol=1e-6, kernel='linear')]

    for clf in clfs:
        clf.set_params(scale_C=False)
        coef_ = clf.fit(X, y).coef_
        coef2_ = clf.fit(X2, y2).coef_
        error_no_scale = linalg.norm(coef2_ - coef_) / linalg.norm(coef_)
        assert_true(error_no_scale > 1e-3)

        clf.set_params(scale_C=True)
        coef_ = clf.fit(X, y).coef_
        coef2_ = clf.fit(X2, y2).coef_
        error_with_scale = linalg.norm(coef2_ - coef_) / linalg.norm(coef_)
        assert_true(error_with_scale < 1e-5)


def test_nu_svc_samples_scaling():
    """Test NuSVC scaling by n_samples
    """
    X = iris.data[iris.target != 2]
    y = iris.target[iris.target != 2]
    X2 = np.r_[X, X]
    y2 = np.r_[y, y]

    clfs = [svm.NuSVC(tol=1e-6, kernel='linear')]

    for clf in clfs:
        coef_ = clf.fit(X, y).coef_
        coef2_ = clf.fit(X2, y2).coef_
        error_with_scale = linalg.norm(coef2_ - coef_) / linalg.norm(coef_)
        assert_true(error_with_scale < 1e-5)


def test_immutable_coef_property():
    """Check that primal coef modification are not silently ignored"""
    svms = [
        svm.SVC(kernel='linear').fit(iris.data, iris.target),
        svm.NuSVC(kernel='linear').fit(iris.data, iris.target),
        svm.SVR(kernel='linear').fit(iris.data, iris.target),
        svm.NuSVR(kernel='linear').fit(iris.data, iris.target),
        svm.OneClassSVM(kernel='linear').fit(iris.data),
        svm.sparse.SVC(kernel='linear').fit(iris.data, iris.target),
        svm.sparse.NuSVC(kernel='linear').fit(iris.data, iris.target),
        svm.sparse.SVR(kernel='linear').fit(iris.data, iris.target),
        svm.sparse.NuSVR(kernel='linear').fit(iris.data, iris.target),
    ]
    for clf in svms:
        assert_raises(AttributeError, clf.__setattr__, 'coef_', np.arange(3))
        assert_raises(RuntimeError, clf.coef_.__setitem__, (0, 0), 0)


def test_inheritance():
    # check that SVC classes can do inheritance
    class ChildSVC(svm.SVC):
        def __init__(self, foo=0):
            self.foo = foo
            svm.SVC.__init__(self)

    clf = ChildSVC()
    clf.fit(iris.data, iris.target)
    clf.predict(iris.data[-1])
    clf.decision_function(iris.data[-1])


if __name__ == '__main__':
    import nose
    nose.runmodule()
