from scikits.learn import svm
import numpy as np
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

# define our data
X = np.array( [[-2,-1], [-1, -1], [-1, -2], [1,1], [1,2], [2, 1]])
y = np.array( [1, 1, 1, 2, 2, 2])
T = np.array( [[-1,-1], [2, 2], [3, 2]] )

def test_svm_params():
    """
    C_SVC algorithm and linear kernel.

    This checks that we retrieve the correct parameters.

    Data is just 6 separable points in the plane. We use linear kernel
    """

    clf =  svm.SVM(kernel_type='linear')
    clf.fit(X, y)

    assert_array_equal(clf.coef_, [[ 0.25, -.25]])
    assert_array_equal(clf.support_, [[-1,-1], [1, 1]])
    assert_array_equal(clf.rho_, [0.])

def test_tweak_params():
    """
    Make sure some tweaking of parameters works.
    Currently this will fail because communication between fit/predict
    is done via a pointer and not copying the arguments.

    We change clf.coef_ at run time and expect .predict() to change
    accordingly. Notice that this is not trivial since it involves a lot
    of C/Python copying in the libsvm bindings.
    """
    clf = svm.SVM(kernel_type='linear')
    clf.fit(X, y)
    assert_array_equal(clf.coef_, [[.25, -.25]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [1])
    clf.coef_ = np.array([[.0, 1.]])
    assert_array_equal(clf.predict([[-.1, -.1]]), [2])

def test_error():
    """
    test that it gives proper exception on deficient input
    """
    # impossible value of nu
    clf = svm.SVM(svm_type='nu_svc', kernel_type='linear', nu=0.0)
    assert_raises(ValueError, clf.fit, X, y)

    y2 = y[:-1] # wrong dimensions for labels
    assert_raises(ValueError, svm.SVM, X, y2)

def test_predict():
    true_result = [1, 2, 2]
    assert_array_equal(svm.predict(X, y, T) , true_result)
    # the same, but using SVM object
    clf = svm.SVM()
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

def test_noncontiguous():
    """
    Test with arrays that are non-contiguous.

    """
    Xt = X.transpose()
    yt = [1, 2]
    assert_array_equal(svm.predict(Xt, yt, T), [1, 2, 2])

def test_dimension_mismatch():
    """
    Test with data that in which dimensions of data space and labels do not
    match
    """
    Y2 = y[:-1]
    assert_raises(AssertionError, svm.predict, X, Y2, T)

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
    result = svm.predict(X, Y, test)
    assert_array_equal(result, [2])

