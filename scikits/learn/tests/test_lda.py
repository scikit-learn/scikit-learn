import numpy as np
from scikits.learn import lda
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises

# Data is just 6 separable points in the plane
X = np.array( [[-2,-1], [-1, -1], [-1, -2], [1,1], [1,2], [2, 1]])
y = np.array( [1, 1, 1, 2, 2, 2])

def test_lda():
    """
    LDA classification.

    This checks that LDA implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf =  lda.LDA()
    y_pred = clf.fit(X, y).predict(X)

    assert_array_equal(y_pred, y)

