import numpy as np
from scikits.learn import svm
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

# Data is just 6 separable points in the plane
X = np.array( [[-2,-1], [-1, -1], [-1, -2], [1,1], [1,2], [2, 1]])
Y = np.array( [1, 1, 1, 2, 2, 2])
T = np.array( [[-1,-1], [2, 2], [3, 2]] )
true_result = [1, 2, 2]

def test_svm_params():
    """
    C_SVC algorithm and linear kernel.

    This checks that we retrieve the correct parameters.
    """

    clf =  svm.SVC(kernel='linear')
    clf.fit(X, Y)

    assert_array_equal(clf.coef_, [[ 0.25, -.25]])
    assert_array_equal(clf.support_, [[-1,-1], [1, 1]])
    assert_array_equal(clf.rho_, [0.])

def test_fit():
    clf = svm.SVC()
    clf.fit([[1,2]], [0])
    assert_array_equal(clf.predict([[-1, -1]]), [0.])

def test_tweak_params():
    """
    Make sure some tweaking of parameters works.

    We change clf.coef_ at run time and expect .predict() to change
    accordingly. Notice that this is not trivial since it involves a lot
    of C/Python copying in the libsvm bindings.

    The success of this test ensures that the mapping between libsvm and
    the python classifier is complete.
    """
    clf = svm.SVC(kernel='linear')
    clf.fit(X, Y)
    assert_array_equal(clf.coef_, [[.25, -.25]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    clf.coef_ = np.array([[.0, 1.]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [2])

def test_error():
    """
    Test that it gives proper exception on deficient input
    """
    # impossible value of nu
    clf = svm.SVC(impl='nu_svc', kernel='linear', nu=0.0)
    assert_raises(ValueError, clf.fit, X, Y)

    Y2 = Y[:-1] # wrong dimensions for labels
    assert_raises(ValueError, svm.SVC, X, Y2)

def test_predict():
    clf = svm.SVC()
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)


def test_probability():
    """
    predict probabilities
    """
    clf = svm.SVC(probability=True)
    clf.fit(X, Y)
    assert_array_almost_equal(clf.predict_proba(T),
                              [[ 0.81997622,  0.18002378],
                               [ 0.18863819,  0.81136181],
                               [ 0.27647055,  0.72352945]])

def test_noncontiguous():
    """
    Test with arrays that are non-contiguous.
    """
    Xt = X.transpose()
    Yt = [1, 2]
    clf = svm.SVC()
    clf.fit(Xt, Yt)
    assert_array_equal(clf.predict(T), [1, 2, 2])

def test_predict_multiclass():
    """
    this example is from libsvm/README
    """
    X  =   [[0,	  0.1,	  0.2,	  0,	  0],
            [ 0,  0.1,	  0.3,	 -1.2,	  0],
            [0.4,  0,	  0,	  0,	  0],
            [0,	  0.1,	  0,	  1.4,	  0.5],
            [-0.1,-0.2,	  0.1,	  1.1,	  0.1]]
    Y = [1,2,1,2,3]
    test = [[0, 1, 0, -1, 0]]
    clf = svm.SVC()
    clf.fit(X, Y)
    result = clf.predict(test)
    assert_array_equal(result, [2])

def test_regression():
    """
    Test SVM regression
    TODO: simplify this. btw, is it correct ?
    """
    clf = svm.SVR()
    clf.fit([[0,0], [1, 1]], [0, 1])
    assert_array_almost_equal(clf.predict([[0,0], [1, 1]]), [.099999, .9])
    

def test_oneclass():
    """
    FIXME: this does nothing
    """
    clf = svm.OneClassSVM()
    clf.fit(X, Y)
    assert_array_equal(Y, [1, 1, 1, 2, 2, 2])

def test_LinearSVC():
    clf = svm.LinearSVC()
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # the same with l1 penalty
    clf = svm.LinearSVC(penalty='l1')
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)

    # l2 penalty with dual formulation
    clf = svm.LinearSVC(penalty='l2', dual=True)
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)
    
    # 
    clf = svm.LinearSVC(penalty='l2', loss='l1', dual=True)
    clf.fit(X, Y)
    assert_array_equal(clf.predict(T), true_result)
