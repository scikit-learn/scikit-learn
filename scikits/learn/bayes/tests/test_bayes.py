import numpy as np
from scikits.learn.bayes.bayes import bayesian_regression
from numpy.testing import assert_array_equal

X = np.array([[1], [2]])
Y = np.array([1, 2])

def test_toy():
    w = bayesian_regression(X, Y)
    assert_array_equal(w, [1])
