import numpy as np
from scikits.learn import logistic
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

X = [[0, 0], [0, 1], [1, 1]]
Y1 = [0,1,1]
Y2 = [0,1,2]


def test_predict_2_classes():
    clf = logistic.LogisticRegression()
    clf.fit(X, Y1)
    assert_array_almost_equal(clf.coef_, [[-0.27501564, -0.60803562]])
    assert_array_almost_equal(clf.intercept_, [-0.08642295])
    assert_array_equal(clf.predict([[-1, -1], [0, 1],]), [0, 1])

    clf = logistic.LogisticRegression(intercept=False)
    clf.fit(X, Y1)
    assert_array_almost_equal(clf.coef_, [[-0.28540916, -0.63236105]])
    assert_array_almost_equal(clf.intercept_, [0])

def test_predict_3_classes():
    clf = logistic.LogisticRegression()
    clf.fit(X, Y2)
    assert_array_equal(clf.predict([[1, 0], [0, 1], [1, 1]]), [2, 1, 2])

def test_predict_proba():
    clf = logistic.LogisticRegression()
    clf.fit(X, Y2)
    assert_array_almost_equal(clf.predict_proba([[1, 1]]),
                              [[ 0.21490268,  0.32639437,  0.45870294]])

    clf = logistic.LogisticRegression(penalty='l1')
    clf.fit(X, Y2)
    assert_array_almost_equal(clf.predict_proba([[2, 2]]),
                              [[ 0.33333333,  0.33333333,  0.33333333]])
