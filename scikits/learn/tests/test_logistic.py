import numpy as np
from scikits.learn import logistic
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

def test_predict_2_features():
    clf = logistic.LogisticRegression()
    clf.fit([[-1, 0], [0, 1]], [0, 1])
    assert_array_equal(clf.predict([[-1, 0], [0, 1],]), [0, 1])

def test_predict_3_features():
    clf = logistic.LogisticRegression()
    clf.fit([[1, 0], [0, 1], [1, 1]], [0, 1, 2])
    assert_array_equal(clf.predict([[1, 0], [0, 1], [1, 1]]), [0, 1, 2])
