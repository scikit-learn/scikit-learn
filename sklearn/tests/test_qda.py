import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_true

from sklearn import qda

# Data is just 6 separable points in the plane
X = np.array([[0, 0], [-2, -2], [-2, -1], [-1, -1], [-1, -2],
              [1, 3], [1, 2], [2, 1], [2, 2]])
y = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2])
y3 = np.array([1, 2, 3, 2, 3, 1, 2, 3, 1])

# Degenerate data with 1 feature (still should be separable)
X1 = np.array([[-3, ], [-2, ], [-1, ], [-1, ], [0, ], [1, ], [1, ],
               [2, ], [3, ]])


def test_qda():
    """
    QDA classification.

    This checks that QDA implements fit and predict and returns
    correct values for a simple toy dataset.
    """
    clf = qda.QDA()
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

    y_pred3 = clf.fit(X, y3).predict(X)
    # QDA shouldn't be able to separate those
    assert_true(np.any(y_pred3 != y3))


def test_qda_priors():
    clf = qda.QDA(priors=np.array([0.0, 1.0]))
    y_pred = clf.fit(X, y).predict(X)
    assert (y_pred == 2).all()


def test_qda_store_covariances():
    # The default is to not set the covariances_ attribute
    clf = qda.QDA().fit(X, y)
    assert_true(not hasattr(clf, 'covariances_'))

    # Test the actual attribute:
    clf = qda.QDA().fit(X, y, store_covariances=True)
    assert_true(hasattr(clf, 'covariances_'))

    assert_array_almost_equal(
        clf.covariances_[0],
        np.array([[0.7, 0.45], [0.45, 0.7]])
    )

    assert_array_almost_equal(
        clf.covariances_[1],
        np.array([[0.33333333, -0.33333333], [-0.33333333, 0.66666667]])
    )
