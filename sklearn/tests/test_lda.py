import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_equal, assert_true

from .. import lda

# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
y = np.array([1, 1, 1, 2, 2, 2])
y3 = np.array([1, 1, 2, 2, 3, 3])

# Degenerate data with 1 feature (still should be separable)
X1 = np.array([[-2, ], [-1, ], [-1, ], [1, ], [1, ], [2, ]])


def test_lda_predict():
    """
    LDA classification.

    This checks that LDA implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf = lda.LDA()
    y_pred = clf.fit(X, y).predict(X)

    assert_array_equal(y_pred, y)

    # Assure that it works with 1D data
    y_pred1 = clf.fit(X1, y).predict(X1)
    assert_array_equal(y_pred1, y)

    # Test probas estimates
    y_proba_pred1 = clf.predict_proba(X1)
    assert_array_equal((y_proba_pred1[:, 1] > 0.5) + 1, y)
    y_log_proba_pred1 = clf.predict_log_proba(X1)
    assert_array_almost_equal(np.exp(y_log_proba_pred1), y_proba_pred1, 8)

    # Primarily test for commit 2f34950 -- "reuse" of priors
    y_pred3 = clf.fit(X, y3).predict(X)
    # LDA shouldn't be able to separate those
    assert_true(np.any(y_pred3 != y3))


def test_lda_transform():
    clf = lda.LDA()
    X_transformed = clf.fit(X, y).transform(X)
    assert_equal(X_transformed.shape[1], 1)
