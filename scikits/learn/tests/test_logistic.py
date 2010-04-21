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
    assert_array_equal(clf.predict([[-1, -1], [0, 1],]), [0, 1])

def test_predict_3_classes():
    clf = logistic.LogisticRegression()
    clf.fit(X, Y2)
    assert_array_equal(clf.predict([[1, 0], [0, 1], [1, 1]]), [2, 1, 2])

def test_predict_proba():
    clf = logistic.LogisticRegression()
    clf.fit(X, Y2)
    assert_array_almost_equal(clf.predict_proba([[1, 1]]),
                              [[ 0.23148573,  0.31760051,  0.45091376]])

    clf = logistic.LogisticRegression(penalty='l1')
    clf.fit(X, Y2)
    assert_array_almost_equal(clf.predict_proba([[2, 2]]),
                              [[ 0.33333333,  0.33333333,  0.33333333]])
