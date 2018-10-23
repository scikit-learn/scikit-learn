import numpy as np

import pytest

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_false
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import ignore_warnings

from sklearn.datasets import make_blobs
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
from sklearn.discriminant_analysis import _cov


# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]], dtype='f')
y = np.array([1, 1, 1, 2, 2, 2])
y3 = np.array([1, 1, 2, 2, 3, 3])

# Degenerate data with only one feature (still should be separable)
X1 = np.array([[-2, ], [-1, ], [-1, ], [1, ], [1, ], [2, ]], dtype='f')

# Data is just 9 separable points in the plane
X6 = np.array([[0, 0], [-2, -2], [-2, -1], [-1, -1], [-1, -2],
               [1, 3], [1, 2], [2, 1], [2, 2]])
y6 = np.array([1, 1, 1, 1, 1, 2, 2, 2, 2])
y7 = np.array([1, 2, 3, 2, 3, 1, 2, 3, 1])

# Degenerate data with 1 feature (still should be separable)
X7 = np.array([[-3, ], [-2, ], [-1, ], [-1, ], [0, ], [1, ], [1, ],
               [2, ], [3, ]])

# Data that has zero variance in one dimension and needs regularization
X2 = np.array([[-3, 0], [-2, 0], [-1, 0], [-1, 0], [0, 0], [1, 0], [1, 0],
               [2, 0], [3, 0]])

# One element class
y4 = np.array([1, 1, 1, 1, 1, 1, 1, 1, 2])

# Data with less samples in a class than n_features
X5 = np.c_[np.arange(8), np.zeros((8, 3))]
y5 = np.array([0, 0, 0, 0, 0, 1, 1, 1])

solver_shrinkage = [('svd', None), ('lsqr', None), ('eigen', None),
                    ('lsqr', 'auto'), ('lsqr', 0), ('lsqr', 0.43),
                    ('eigen', 'auto'), ('eigen', 0), ('eigen', 0.43)]


def test_lda_predict():
    # Test LDA classification.
    # This checks that LDA implements fit and predict and returns correct
    # values for simple toy data.
    for test_case in solver_shrinkage:
        solver, shrinkage = test_case
        clf = LinearDiscriminantAnalysis(solver=solver, shrinkage=shrinkage)
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
        assert np.any(y_pred3 != y3), 'solver %s' % solver

    # Test invalid shrinkages
    clf = LinearDiscriminantAnalysis(solver="lsqr", shrinkage=-0.2231)
    assert_raises(ValueError, clf.fit, X, y)
    clf = LinearDiscriminantAnalysis(solver="eigen", shrinkage="dummy")
    assert_raises(ValueError, clf.fit, X, y)
    clf = LinearDiscriminantAnalysis(solver="svd", shrinkage="auto")
    assert_raises(NotImplementedError, clf.fit, X, y)
    # Test unknown solver
    clf = LinearDiscriminantAnalysis(solver="dummy")
    assert_raises(ValueError, clf.fit, X, y)


def test_lda_priors():
    # Test priors (negative priors)
    priors = np.array([0.5, -0.5])
    clf = LinearDiscriminantAnalysis(priors=priors)
    msg = "priors must be non-negative"
    assert_raise_message(ValueError, msg, clf.fit, X, y)

    # Test that priors passed as a list are correctly handled (run to see if
    # failure)
    clf = LinearDiscriminantAnalysis(priors=[0.5, 0.5])
    clf.fit(X, y)

    # Test that priors always sum to 1
    priors = np.array([0.5, 0.6])
    prior_norm = np.array([0.45, 0.55])
    clf = LinearDiscriminantAnalysis(priors=priors)
    assert_warns(UserWarning, clf.fit, X, y)
    assert_array_almost_equal(clf.priors_, prior_norm, 2)


def test_lda_coefs():
    # Test if the coefficients of the solvers are approximately the same.
    n_features = 2
    n_classes = 2
    n_samples = 1000
    X, y = make_blobs(n_samples=n_samples, n_features=n_features,
                      centers=n_classes, random_state=11)

    clf_lda_svd = LinearDiscriminantAnalysis(solver="svd")
    clf_lda_lsqr = LinearDiscriminantAnalysis(solver="lsqr")
    clf_lda_eigen = LinearDiscriminantAnalysis(solver="eigen")

    clf_lda_svd.fit(X, y)
    clf_lda_lsqr.fit(X, y)
    clf_lda_eigen.fit(X, y)

    assert_array_almost_equal(clf_lda_svd.coef_, clf_lda_lsqr.coef_, 1)
    assert_array_almost_equal(clf_lda_svd.coef_, clf_lda_eigen.coef_, 1)
    assert_array_almost_equal(clf_lda_eigen.coef_, clf_lda_lsqr.coef_, 1)


def test_lda_transform():
    # Test LDA transform.
    clf = LinearDiscriminantAnalysis(solver="svd", n_components=1)
    X_transformed = clf.fit(X, y).transform(X)
    assert_equal(X_transformed.shape[1], 1)
    clf = LinearDiscriminantAnalysis(solver="eigen", n_components=1)
    X_transformed = clf.fit(X, y).transform(X)
    assert_equal(X_transformed.shape[1], 1)

    clf = LinearDiscriminantAnalysis(solver="lsqr", n_components=1)
    clf.fit(X, y)
    msg = "transform not implemented for 'lsqr'"
    assert_raise_message(NotImplementedError, msg, clf.transform, X)


def test_lda_explained_variance_ratio():
    # Test if the sum of the normalized eigen vectors values equals 1,
    # Also tests whether the explained_variance_ratio_ formed by the
    # eigen solver is the same as the explained_variance_ratio_ formed
    # by the svd solver

    state = np.random.RandomState(0)
    X = state.normal(loc=0, scale=100, size=(40, 20))
    y = state.randint(0, 3, size=(40,))

    clf_lda_eigen = LinearDiscriminantAnalysis(solver="eigen")
    clf_lda_eigen.fit(X, y)
    assert_almost_equal(clf_lda_eigen.explained_variance_ratio_.sum(), 1.0, 3)
    assert_equal(clf_lda_eigen.explained_variance_ratio_.shape, (2,),
                 "Unexpected length for explained_variance_ratio_")

    clf_lda_svd = LinearDiscriminantAnalysis(solver="svd")
    clf_lda_svd.fit(X, y)
    assert_almost_equal(clf_lda_svd.explained_variance_ratio_.sum(), 1.0, 3)
    assert_equal(clf_lda_svd.explained_variance_ratio_.shape, (2,),
                 "Unexpected length for explained_variance_ratio_")

    assert_array_almost_equal(clf_lda_svd.explained_variance_ratio_,
                              clf_lda_eigen.explained_variance_ratio_)


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
    clf = LinearDiscriminantAnalysis(solver="svd").fit(X, y)
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
        clf = LinearDiscriminantAnalysis(solver=solver)
        # should be able to separate the data perfectly
        assert_equal(clf.fit(x, y).score(x, y), 1.0,
                     'using covariance: %s' % solver)


def test_lda_store_covariance():
    # Test for slover 'lsqr' and 'eigen'
    # 'store_covariance' has no effect on 'lsqr' and 'eigen' solvers
    for solver in ('lsqr', 'eigen'):
        clf = LinearDiscriminantAnalysis(solver=solver).fit(X6, y6)
        assert hasattr(clf, 'covariance_')

        # Test the actual attribute:
        clf = LinearDiscriminantAnalysis(solver=solver,
                                         store_covariance=True).fit(X6, y6)
        assert hasattr(clf, 'covariance_')

        assert_array_almost_equal(
            clf.covariance_,
            np.array([[0.422222, 0.088889], [0.088889, 0.533333]])
        )

    # Test for SVD slover, the default is to not set the covariances_ attribute
    clf = LinearDiscriminantAnalysis(solver='svd').fit(X6, y6)
    assert_false(hasattr(clf, 'covariance_'))

    # Test the actual attribute:
    clf = LinearDiscriminantAnalysis(solver=solver,
                                     store_covariance=True).fit(X6, y6)
    assert hasattr(clf, 'covariance_')

    assert_array_almost_equal(
        clf.covariance_,
        np.array([[0.422222, 0.088889], [0.088889, 0.533333]])
    )


def test_qda():
    # QDA classification.
    # This checks that QDA implements fit and predict and returns
    # correct values for a simple toy dataset.
    clf = QuadraticDiscriminantAnalysis()
    y_pred = clf.fit(X6, y6).predict(X6)
    assert_array_equal(y_pred, y6)

    # Assure that it works with 1D data
    y_pred1 = clf.fit(X7, y6).predict(X7)
    assert_array_equal(y_pred1, y6)

    # Test probas estimates
    y_proba_pred1 = clf.predict_proba(X7)
    assert_array_equal((y_proba_pred1[:, 1] > 0.5) + 1, y6)
    y_log_proba_pred1 = clf.predict_log_proba(X7)
    assert_array_almost_equal(np.exp(y_log_proba_pred1), y_proba_pred1, 8)

    y_pred3 = clf.fit(X6, y7).predict(X6)
    # QDA shouldn't be able to separate those
    assert np.any(y_pred3 != y7)

    # Classes should have at least 2 elements
    assert_raises(ValueError, clf.fit, X6, y4)


def test_qda_priors():
    clf = QuadraticDiscriminantAnalysis()
    y_pred = clf.fit(X6, y6).predict(X6)
    n_pos = np.sum(y_pred == 2)

    neg = 1e-10
    clf = QuadraticDiscriminantAnalysis(priors=np.array([neg, 1 - neg]))
    y_pred = clf.fit(X6, y6).predict(X6)
    n_pos2 = np.sum(y_pred == 2)

    assert_greater(n_pos2, n_pos)


def test_qda_store_covariance():
    # The default is to not set the covariances_ attribute
    clf = QuadraticDiscriminantAnalysis().fit(X6, y6)
    assert_false(hasattr(clf, 'covariance_'))

    # Test the actual attribute:
    clf = QuadraticDiscriminantAnalysis(store_covariance=True).fit(X6, y6)
    assert hasattr(clf, 'covariance_')

    assert_array_almost_equal(
        clf.covariance_[0],
        np.array([[0.7, 0.45], [0.45, 0.7]])
    )

    assert_array_almost_equal(
        clf.covariance_[1],
        np.array([[0.33333333, -0.33333333], [-0.33333333, 0.66666667]])
    )


def test_qda_deprecation():
    # Test the deprecation
    clf = QuadraticDiscriminantAnalysis(store_covariances=True)
    assert_warns_message(DeprecationWarning, "'store_covariances' was renamed"
                         " to store_covariance in version 0.19 and will be "
                         "removed in 0.21.", clf.fit, X, y)

    # check that covariance_ (and covariances_ with warning) is stored
    assert_warns_message(DeprecationWarning, "Attribute ``covariances_`` was "
                         "deprecated in version 0.19 and will be removed "
                         "in 0.21. Use ``covariance_`` instead", getattr, clf,
                         'covariances_')


def test_qda_regularization():
    # the default is reg_param=0. and will cause issues
    # when there is a constant variable
    clf = QuadraticDiscriminantAnalysis()
    with ignore_warnings():
        y_pred = clf.fit(X2, y6).predict(X2)
    assert np.any(y_pred != y6)

    # adding a little regularization fixes the problem
    clf = QuadraticDiscriminantAnalysis(reg_param=0.01)
    with ignore_warnings():
        clf.fit(X2, y6)
    y_pred = clf.predict(X2)
    assert_array_equal(y_pred, y6)

    # Case n_samples_in_a_class < n_features
    clf = QuadraticDiscriminantAnalysis(reg_param=0.1)
    with ignore_warnings():
        clf.fit(X5, y5)
    y_pred5 = clf.predict(X5)
    assert_array_equal(y_pred5, y5)


def test_covariance():
    x, y = make_blobs(n_samples=100, n_features=5,
                      centers=1, random_state=42)

    # make features correlated
    x = np.dot(x, np.arange(x.shape[1] ** 2).reshape(x.shape[1], x.shape[1]))

    c_e = _cov(x, 'empirical')
    assert_almost_equal(c_e, c_e.T)

    c_s = _cov(x, 'auto')
    assert_almost_equal(c_s, c_s.T)


@pytest.mark.parametrize("solver", ['svd, lsqr', 'eigen'])
def test_raises_value_error_on_same_number_of_classes_and_samples(solver):
    """
    Tests that if the number of samples equals the number
    of classes, a ValueError is raised.
    """
    X = np.array([[0.5, 0.6], [0.6, 0.5]])
    y = np.array(["a", "b"])
    clf = LinearDiscriminantAnalysis(solver=solver)
    with pytest.raises(ValueError, match="The number of samples must be more"):
        clf.fit(X, y)
