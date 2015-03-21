import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message

from sklearn.datasets import make_blobs
from sklearn import lda


# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]], dtype='f')
y = np.array([1, 1, 1, 2, 2, 2])
y3 = np.array([1, 1, 2, 2, 3, 3])

# Degenerate data with only one feature (still should be separable)
X1 = np.array([[-2, ], [-1, ], [-1, ], [1, ], [1, ], [2, ]], dtype='f')

solver_shrinkage = [('svd', None), ('lsqr', None), ('eigen', None),
                    ('lsqr', 'auto'), ('lsqr', 0), ('lsqr', 0.43),
                    ('eigen', 'auto'), ('eigen', 0), ('eigen', 0.43)]


def test_lda_predict():
    # Test LDA classification.
    # This checks that LDA implements fit and predict and returns correct values
    # for simple toy data.
    for test_case in solver_shrinkage:
        solver, shrinkage = test_case
        clf = lda.LDA(solver=solver, shrinkage=shrinkage)
        y_pred = clf.fit(X, y).predict(X)
        assert_array_equal(y_pred, y, 'solver %s' % solver)

        # Assert that it works with 1D data
        y_pred1 = clf.fit(X1, y).predict(X1)
        assert_array_equal(y_pred1, y, 'solver %s' % solver)

        # Test probability estimates
        y_proba_pred1 = clf.predict_proba(X1)
        assert_array_equal((y_proba_pred1[:, 1] > 0.5) + 1, y,
                           'solver %s' % solver)
        y_log_proba_pred1 = clf.predict_log_proba(X1)
        assert_array_almost_equal(np.exp(y_log_proba_pred1), y_proba_pred1,
                                  8, 'solver %s' % solver)

        # Primarily test for commit 2f34950 -- "reuse" of priors
        y_pred3 = clf.fit(X, y3).predict(X)
        # LDA shouldn't be able to separate those
        assert_true(np.any(y_pred3 != y3), 'solver %s' % solver)

    # Test invalid shrinkages
    clf = lda.LDA(solver="lsqr", shrinkage=-0.2231)
    assert_raises(ValueError, clf.fit, X, y)
    clf = lda.LDA(solver="eigen", shrinkage="dummy")
    assert_raises(ValueError, clf.fit, X, y)
    clf = lda.LDA(solver="svd", shrinkage="auto")
    assert_raises(NotImplementedError, clf.fit, X, y)
    # Test unknown solver
    clf = lda.LDA(solver="dummy")
    assert_raises(ValueError, clf.fit, X, y)


def test_lda_coefs():
    # Test if the coefficients of the solvers are approximately the same.
    n_features = 2
    n_classes = 2
    n_samples = 1000
    X, y = make_blobs(n_samples=n_samples, n_features=n_features,
                      centers=n_classes, random_state=11)

    clf_lda_svd = lda.LDA(solver="svd")
    clf_lda_lsqr = lda.LDA(solver="lsqr")
    clf_lda_eigen = lda.LDA(solver="eigen")

    clf_lda_svd.fit(X, y)
    clf_lda_lsqr.fit(X, y)
    clf_lda_eigen.fit(X, y)

    assert_array_almost_equal(clf_lda_svd.coef_, clf_lda_lsqr.coef_, 1)
    assert_array_almost_equal(clf_lda_svd.coef_, clf_lda_eigen.coef_, 1)
    assert_array_almost_equal(clf_lda_eigen.coef_, clf_lda_lsqr.coef_, 1)


def test_lda_transform():
    # Test LDA transform.
    clf = lda.LDA(solver="svd", n_components=1)
    X_transformed = clf.fit(X, y).transform(X)
    assert_equal(X_transformed.shape[1], 1)
    clf = lda.LDA(solver="eigen", n_components=1)
    X_transformed = clf.fit(X, y).transform(X)
    assert_equal(X_transformed.shape[1], 1)

    clf = lda.LDA(solver="lsqr", n_components=1)
    clf.fit(X, y)
    msg = "transform not implemented for 'lsqr'"
    assert_raise_message(NotImplementedError, msg, clf.transform, X)


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
    clf = lda.LDA(solver="svd").fit(X, y)
    means_transformed = clf.transform(means)

    d1 = means_transformed[3] - means_transformed[0]
    d2 = means_transformed[2] - means_transformed[1]
    d1 /= np.sqrt(np.sum(d1 ** 2))
    d2 /= np.sqrt(np.sum(d2 ** 2))

    # the transformed within-class covariance should be the identity matrix
    assert_almost_equal(np.cov(clf.transform(scatter).T), np.eye(2))

    # the means of classes 0 and 3 should lie on the first component
    assert_almost_equal(np.abs(np.dot(d1[:2], [1, 0])), 1.0)

    # the means of classes 1 and 2 should lie on the second component
    assert_almost_equal(np.abs(np.dot(d2[:2], [0, 1])), 1.0)


def test_lda_scaling():
    # Test if classification works correctly with differently scaled features.
    n = 100
    rng = np.random.RandomState(1234)
    # use uniform distribution of features to make sure there is absolutely no
    # overlap between classes.
    x1 = rng.uniform(-1, 1, (n, 3)) + [-10, 0, 0]
    x2 = rng.uniform(-1, 1, (n, 3)) + [10, 0, 0]
    x = np.vstack((x1, x2)) * [1, 100, 10000]
    y = [-1] * n + [1] * n

    for solver in ('svd', 'lsqr', 'eigen'):
        clf = lda.LDA(solver=solver)
        # should be able to separate the data perfectly
        assert_equal(clf.fit(x, y).score(x, y), 1.0,
                     'using covariance: %s' % solver)
