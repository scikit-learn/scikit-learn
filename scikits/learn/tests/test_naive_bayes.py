import numpy as np
from numpy.testing import assert_array_equal

from .. import naive_bayes

# Data is just 6 separable points in the plane
X = np.array( [[-2,-1], [-1, -1], [-1, -2], [1,1], [1,2], [2, 1]])
y = np.array( [1, 1, 1, 2, 2, 2])

def test_gnb():
    """
    Gaussian Naive Bayes classification.

    This checks that GNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf =  naive_bayes.GNB()
    y_pred = clf.fit(X, y).predict(X)
    assert_array_equal(y_pred, y)

