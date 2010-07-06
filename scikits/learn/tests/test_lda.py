import numpy as np
from scikits.learn import lda
from numpy.testing import assert_array_equal, \
                          assert_array_almost_equal, \
                          assert_raises
from nose.tools import assert_true

# Data is just 6 separable points in the plane
X = np.array( [[-2,-1], [-1, -1], [-1, -2], [1,1], [1,2], [2, 1]])
y = np.array( [1, 1, 1, 2, 2, 2])
y3 = np.array([1, 1, 2, 2, 3, 3])

# Degenerate data with 1 feature (still should be separable)
X1 = np.array( [[-2,], [-1,], [-1,], [1,], [1,], [2,]])

def test_lda():
    """
    LDA classification.

    This checks that LDA implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf =  lda.LDA()
    y_pred = clf.fit(X, y).predict(X)

    assert_array_equal(y_pred, y)

    # Assure that it works with 1D data
    y_pred1 = clf.fit(X1, y).predict(X1)
    assert_array_equal(y_pred1, y)

    # Primarily test for commit 2f34950 -- "reuse" of priors
    y_pred3 = clf.fit(X, y3).predict(X)
    # LDA shouldn't be able to separate those
    assert_true(np.any(y_pred3 != y3))

