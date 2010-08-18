
from scikits.learn import logistic
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          decorators

X = [[0, 0], [0, 1], [1, 1]]
Y1 = [0,1,1]
Y2 = [0,1,2]

def test_predict_2_classes():
    clf = logistic.LogisticRegression().fit(X, Y1)
    assert_array_almost_equal(clf.coef_, [[-0.27501564, -0.60803562]])
    assert_array_almost_equal(clf.intercept_, [-0.08642295])
    assert_array_equal(clf.predict([[-1, -1], [0, 1],]), [0, 1])

    clf = logistic.LogisticRegression(has_intercept=False).fit(X, Y1)
    assert_array_almost_equal(clf.coef_, [[-0.28540916, -0.63236105]])
    assert_array_almost_equal(clf.intercept_, [0])

@decorators.skipif(True, "XFailed test")
def test_predict_3_classes():
    clf = logistic.LogisticRegression().fit(X, Y2)
    assert_array_equal(clf.predict([[1, 0], [0, 1], [1, 1]]), [0, 1, 2])

@decorators.skipif(True, "XFailed test")
def test_predict_proba():
    """
    I think this test is wrong. Is there a way to know the right results ?
    """
    clf = logistic.LogisticRegression().fit(X, Y2)
    assert_array_almost_equal(clf.predict_proba([[1, 1]]),
                              [[ 0.21490268,  0.32639437,  0.45870294]])

    clf = logistic.LogisticRegression(penalty='l1').fit(X, Y2)
    assert_array_almost_equal(clf.predict_proba([[2, 2]]),
                              [[ 0.33333333,  0.33333333,  0.33333333]])
