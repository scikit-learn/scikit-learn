"""
Testing for Support Vector Machine module (scikits.learn.svm)

TODO: remove hard coded numerical results when possible
"""

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
                          assert_almost_equal
from nose.tools import assert_raises

from scikits.learn import svm, datasets

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [1, 1, 1, 2, 2, 2]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [1, 2, 2]

# also load the iris dataset
iris = datasets.load_iris()


def test_libsvm_parameters():
    """
    Test parameters on classes that make use of libsvm.
    """

    clf = svm.SVC(kernel='linear').fit(X, Y)
    assert_array_equal(clf.dual_coef_, [[ 0.25, -.25]])
    assert_array_equal(clf.support_, [[-1, -1], [1, 1]])
    assert_array_equal(clf.intercept_, [0.])
    assert_array_equal(clf.predict(X), Y)


def test_libsvm_iris():
    """
    Check consistency on dataset iris.
    """
    for k in ('linear', 'rbf'):
        clf = svm.SVC(kernel=k).fit(iris.data, iris.target)
        assert np.mean(clf.predict(iris.data) == iris.target) > 0.9
    

def test_precomputed():
    """
    SVC with a precomputed kernel.

    We test it with a toy dataset and with iris.
    """
    clf = svm.SVC(kernel='precomputed')
    # we use just a linear kernel
    K = np.dot(X, np.array(X).T)
    clf.fit(K, Y)
    # KT is the Gram matrix
    KT = np.dot(T, np.array(X).T)
    pred = clf.predict(KT)

    assert_array_equal(clf.dual_coef_, [[0.25, -.25]])
    assert_array_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.support_, [[2], [4]])
    assert_array_equal(pred, true_result)

    # same as before, but using a callable function instead of the kernel
    # matrix. kernel is just a linear kernel

    kfunc = lambda x, y: np.dot(x, y.T)
    clf = svm.SVC(kernel=kfunc)
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_equal(clf.dual_coef_, [[0.25, -.25]])
    assert_array_equal(clf.intercept_, [0])
    assert_array_almost_equal(clf.support_, [[2], [4]])
    assert_array_equal(pred, true_result)

    # test a precomputed kernel with the iris dataset
    clf = svm.SVC(kernel='precomputed')
    K = np.dot(iris.data, iris.data.T)
    clf.fit(K, iris.target)
    pred = clf.predict(K)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)

    clf = svm.SVC(kernel=kfunc)
    clf.fit(iris.data, iris.target)
    assert_almost_equal(np.mean(pred == iris.target), .99, decimal=2)


def test_SVR():
    """
    Test Support Vector Regression
    """

    clf = svm.SVR(kernel='linear')
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_almost_equal(clf.dual_coef_, [[-0.1, 0.1]])
    assert_array_almost_equal(clf.coef_, [[0.2, 0.2]])
    assert_array_almost_equal(clf.support_, [[-1, -1], [1, 1]])
    assert_array_almost_equal(clf.intercept_, [1.5])
    assert_array_almost_equal(pred, [1.1, 2.3, 2.5])

    # the same with kernel='rbf'
    clf = svm.SVR(kernel='rbf')
    clf.fit(X, Y)
    pred = clf.predict(T)

    assert_array_almost_equal(clf.dual_coef_,
                              [[-0.014, -0.515, -0.013, 0.515, 0.013, 0.013]],
                              decimal=3)
    assert_raises(NotImplementedError, lambda: clf.coef_)
    assert_array_almost_equal(clf.support_, X)
    assert_array_almost_equal(clf.intercept_, [ 1.49997261])
    assert_array_almost_equal(pred, [ 1.10001274, 1.86682485, 1.73300377])


def test_oneclass():
    """
    Test OneClassSVM
    """
    clf = svm.OneClassSVM()
    clf.fit(X)
    pred = clf.predict(T)

    assert_array_almost_equal(pred, [1, -1, -1])
    assert_array_almost_equal(clf.intercept_, [-1.351], decimal=3)
    assert_array_almost_equal(clf.dual_coef_, [[ 0.750, 0.749, 0.749, 0.750]],
                              decimal=3)
    assert_raises(NotImplementedError, lambda: clf.coef_)


def test_tweak_params():
    """
    Make sure some tweaking of parameters works.

    We change clf.dual_coef_ at run time and expect .predict() to change
    accordingly. Notice that this is not trivial since it involves a lot
    of C/Python copying in the libsvm bindings.

    The success of this test ensures that the mapping between libsvm and
    the python classifier is complete.
    """
    clf = svm.SVC(kernel='linear')
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

    clf = svm.SVC(probability=True)
    clf.fit(iris.data, iris.target)

    # predict on a simple dataset
    T = [[0, 0, 0, 0],
         [2, 2, 2, 2]]
    assert_array_almost_equal(clf.predict_proba(T),
                [[ 0.993, 0.003, 0.002],
                 [ 0.740, 0.223  , 0.035]],
                 decimal=2)

    # make sure probabilities sum to one
    pprob = clf.predict_proba(X)
    assert_array_almost_equal(pprob.sum(axis=1),
                               np.ones(len(X)))


def test_margin():
    """
    Test predict_margin

    We create a set of points lying in two lines, so that margin is easily
    calculated in a linear kernel.

    TODO: distance should be sqrt(2)/2, but libsvm returns 1.
    """
    X = [(i, i) for i in range(-4, 6)]
    X += [(i, i+2) for i in range(-4, 6)]
    Y = [0]*10 + [1]*10
    T = [[1]]*10 + [[-1]]*10
    clf = svm.SVC(kernel='linear').fit(X, Y)
    assert_array_almost_equal(clf.predict_margin(X), T)

    # the same using a callable kernel
    kfunc = lambda x, y: np.dot(x, y.T)
    clf = svm.SVC(kernel=kfunc).fit(X, Y)
    assert_array_almost_equal(clf.predict_margin(X), T)

    # failing test
    # assert_array_almost_equal (clf.predict_margin(iris.data),
    #                            np.dot(clf.coef_.T, iris.data) + \
    #                            clf.intercept_)


def test_weight():
    """
    Test class weights
    """
    clf = svm.SVC()
    # we give a small weights to class 1
    clf.fit(X, Y, {1: 0.1})
    # so all predicted values belong to class 2
    assert_array_almost_equal(clf.predict(X), [2] * 6)


def test_error():
    """
    Test that it gives proper exception on deficient input
    """
    # impossible value of C
    assert_raises (ValueError, svm.SVC(C=-1).fit, X, Y)

    # impossible value of nu
    clf = svm.NuSVC(nu=0.0)
    assert_raises(ValueError, clf.fit, X, Y)

    Y2 = Y[:-1] # wrong dimensions for labels
    assert_raises(ValueError, clf.fit, X, Y2)
    assert_raises(AssertionError, svm.SVC, X, Y2)

    # Test with arrays that are non-contiguous.
    Xf = np.asfortranarray(X)
    clf = svm.SVC()
    clf.fit(Xf, Y)
    assert_array_equal(clf.predict(T), true_result)


def test_LinearSVC():
    """
    Test basic routines using LinearSVC
    """
    clf = svm.LinearSVC().fit(X, Y)

    # by default should have intercept
    assert clf.fit_intercept

    assert_array_equal(clf.predict(T), true_result)
    assert_array_almost_equal(clf.intercept_, [0], decimal=5)

    # the same with l1 penalty
    clf = svm.LinearSVC(penalty='l1', dual=False).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty with dual formulation
    clf = svm.LinearSVC(penalty='l2', dual=True).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty, l1 loss
    clf = svm.LinearSVC(penalty='l2', loss='l1', dual=True).fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)


def test_LinearSVC_iris():
    """
    Test that LinearSVC gives plausible predictions on the iris dataset
    """
    clf = svm.LinearSVC().fit(iris.data, iris.target)
    assert np.mean(clf.predict(iris.data) == iris.target) > 0.95

