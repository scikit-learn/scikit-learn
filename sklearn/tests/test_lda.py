import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true

from sklearn import lda

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


def test_lda_orthogonality():
    # arrange four classes with their means in a kite-shaped pattern
    # the longer distance should be transformed to the first component, and
    # the shorter distance to the second component.
    means = np.array([[0, 0, -1], [0, 2, 0], [0, -2, 0], [0, 0, 5]])

    # We construct perfectly symmetric distributions, so the LDA can estimate
    # precise means.
    scatter = np.array([[0.1, 0, 0], [-0.1, 0, 0], [0, 0.1, 0], [0, -0.1, 0],
                        [0, 0, 0.1], [0, 0, -0.1]])

    X = (means[:, np.newaxis, :] + scatter[np.newaxis, :, :]).reshape((-1, 3))
    y = np.repeat(np.arange(means.shape[0]), scatter.shape[0])

    # Fit LDA and transform the means
    clf = lda.LDA().fit(X, y)
    means_transformed = clf.transform(means)

    d1 = means_transformed[3] - means_transformed[0]
    d2 = means_transformed[2] - means_transformed[1]
    d1 /= np.sqrt(np.sum(d1**2))
    d2 /= np.sqrt(np.sum(d2**2))

    # the transformed within-class covariance should be the identity matrix
    assert_almost_equal(np.cov(clf.transform(scatter).T), np.eye(2))

    # the means of classes 0 and 3 should lie on the first component
    assert_almost_equal(np.abs(np.dot(d1[:2], [1, 0])), 1.0)

    # the means of classes 1 and 2 should lie on the second component
    assert_almost_equal(np.abs(np.dot(d2[:2], [0, 1])), 1.0)
