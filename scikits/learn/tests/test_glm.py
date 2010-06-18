import numpy as np
from numpy.testing import *

from scikits.learn import glm

X = [[1, 0, -1.],
     [0, 0, 0],
     [0, 1, 1.]]
Y = [1, 0, -1]

def test_1():
    """
    Very simple test
    """
    clf = glm.LeastAngleRegression().fit(X, Y)
    assert_array_almost_equal(clf.coef_, [0, 0, -1.4142], decimal=4)
    assert_array_almost_equal(clf.alphas_.shape, clf.coef_path_.shape[1])
    assert_array_almost_equal(clf.predict(X), np.sqrt(2)* np.array(Y))

def test_2():

    clf = glm.LeastAngleRegression().fit(X, Y, n_features=1)
    assert_array_almost_equal(clf.coef_, [0, 0, -1.4142], decimal=4)
    print clf.coef_path_
    assert_array_almost_equal(clf.alphas_.shape, clf.coef_path_.shape[1])
    


