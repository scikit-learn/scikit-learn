import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_raises

from .. import slda

# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
y1 = np.array([1, 1, 1, 1, 1, 1])
y2 = np.array([1, 1, 1, 2, 2, 2])
y3 = np.array([1, 1, 2, 2, 3, 3])

# Degenerate data with 1 feature (still should be separable)
X1 = np.array([[-2, ], [-1, ], [-1, ], [1, ], [1, ], [2, ]])


def test_lda_predict():
    """
    SLDA classification.

    This checks that SLDA implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    # Test invalid skrinkage
    clf = slda.SLDA(shrinkage='000')
    y_pred = clf.fit(X, y2).predict(X)
    assert_array_equal(y_pred, y2)

    # Test without shrinkage
    clf = slda.SLDA(shrinkage=None)
    y_pred = clf.fit(X, y2).predict(X)
    assert_array_equal(y_pred, y2)

    # Test with shrinkage
    clf = slda.SLDA(shrinkage='auto')
    y_pred = clf.fit(X, y2).predict(X)
    assert_array_equal(y_pred, y2)

    # Assure that it works with 1D data
    y_pred1 = clf.fit(X1, y2).predict(X1)
    assert_array_equal(y_pred1, y2)

    # Test decision function (two-class case)
    clf = slda.SLDA().fit(X, y2)
    y_decf = clf.decision_function(X)
    y_logp = clf.predict_log_proba(X)
    assert_array_almost_equal(y_decf, y_logp[:, 1] - y_logp[:, 0], 8)

    # Test probability estimates
    y_proba_pred1 = clf.predict_proba(X1)
    assert_array_equal((y_proba_pred1[:, 1] > 0.5) + 1, y2)
    y_log_proba_pred1 = clf.predict_log_proba(X1)
    assert_array_almost_equal(np.exp(y_log_proba_pred1), y_proba_pred1, 8)

    # Primarily test for commit 2f34950 -- "reuse" of priors
    y_pred3 = clf.fit(X, y3).predict(X)
    # LDA shouldn't be able to separate those
    assert_true(np.any(y_pred3 != y3))

    # Test priors
    assert_raises(ValueError, slda.SLDA, priors=[2, -1])
    clf = slda.SLDA(priors=[2, 1])
    y_pred = clf.fit(X, y2).predict(X)
    assert_array_equal(y_pred, y2)

    # Test one class
    assert_raises(ValueError, clf.fit, X, y1)