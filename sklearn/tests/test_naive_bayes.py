import pickle
from io import BytesIO
import numpy as np
import scipy.sparse

import warnings

from sklearn.datasets import load_digits
from sklearn.cross_validation import cross_val_score

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater

from sklearn.naive_bayes import GaussianNB, BernoulliNB, MultinomialNB

# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
y = np.array([1, 1, 1, 2, 2, 2])

# A bit more random tests
rng = np.random.RandomState(0)
X1 = rng.normal(size=(10, 3))
y1 = (rng.normal(size=(10)) > 0).astype(np.int)

# Data is 6 random integer points in a 100 dimensional space classified to
# three classes.
X2 = rng.randint(5, size=(6, 100))
y2 = np.array([1, 1, 2, 2, 3, 3])


def test_gnb():
    """
    Gaussian Naive Bayes classification.

    This checks that GaussianNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf = GaussianNB()
    y_pred = clf.fit(X, y).predict(X)
    assert_array_equal(y_pred, y)

    y_pred_proba = clf.predict_proba(X)
    y_pred_log_proba = clf.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)


def test_gnb_prior():
    """Test whether class priors are properly set. """
    clf = GaussianNB().fit(X, y)
    assert_array_almost_equal(np.array([3, 3]) / 6.0,
                              clf.class_prior_, 8)
    clf.fit(X1, y1)
    # Check that the class priors sum to 1
    assert_array_almost_equal(clf.class_prior_.sum(), 1)


def test_discrete_prior():
    """Test whether class priors are properly set. """
    for cls in [BernoulliNB, MultinomialNB]:
        clf = cls().fit(X2, y2)
        assert_array_almost_equal(np.log(np.array([2, 2, 2]) / 6.0),
                                  clf.class_log_prior_, 8)


def test_mnnb():
    """Test Multinomial Naive Bayes classification.

    This checks that MultinomialNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    for X in [X2, scipy.sparse.csr_matrix(X2)]:
        # Check the ability to predict the learning set.
        clf = MultinomialNB()
        assert_raises(ValueError, clf.fit, -X, y2)
        y_pred = clf.fit(X, y2).predict(X)

        assert_array_equal(y_pred, y2)

        # Verify that np.log(clf.predict_proba(X)) gives the same results as
        # clf.predict_log_proba(X)
        y_pred_proba = clf.predict_proba(X)
        y_pred_log_proba = clf.predict_log_proba(X)
        assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)

        # Check that incremental fitting yields the same results
        clf2 = MultinomialNB()
        clf2.partial_fit(X[:2], y2[:2], classes=np.unique(y2))
        clf2.partial_fit(X[2:5], y2[2:5])
        clf2.partial_fit(X[5:], y2[5:])

        y_pred2 = clf2.predict(X)
        assert_array_equal(y_pred2, y2)

        y_pred_proba2 = clf2.predict_proba(X)
        y_pred_log_proba2 = clf2.predict_log_proba(X)
        assert_array_almost_equal(np.log(y_pred_proba2), y_pred_log_proba2, 8)
        assert_array_almost_equal(y_pred_proba2, y_pred_proba)
        assert_array_almost_equal(y_pred_log_proba2, y_pred_log_proba)

        # Partial fit on the whole data at once should be the same as fit too
        clf3 = MultinomialNB()
        clf3.partial_fit(X, y2, classes=np.unique(y2))

        y_pred3 = clf3.predict(X)
        assert_array_equal(y_pred3, y2)
        y_pred_proba3 = clf3.predict_proba(X)
        y_pred_log_proba3 = clf3.predict_log_proba(X)
        assert_array_almost_equal(np.log(y_pred_proba3), y_pred_log_proba3, 8)
        assert_array_almost_equal(y_pred_proba3, y_pred_proba)
        assert_array_almost_equal(y_pred_log_proba3, y_pred_log_proba)


def check_partial_fit(cls):
    clf1 = cls()
    clf1.fit([[0, 1], [1, 0]], [0, 1])

    clf2 = cls()
    clf2.partial_fit([[0, 1], [1, 0]], [0, 1], classes=[0, 1])
    assert_array_equal(clf1.class_count_, clf2.class_count_)
    assert_array_equal(clf1.feature_count_, clf2.feature_count_)

    clf3 = cls()
    clf3.partial_fit([[0, 1]], [0], classes=[0, 1])
    clf3.partial_fit([[1, 0]], [1])
    assert_array_equal(clf1.class_count_, clf3.class_count_)
    assert_array_equal(clf1.feature_count_, clf3.feature_count_)


def test_discretenb_partial_fit():
    for cls in [MultinomialNB, BernoulliNB]:
        yield check_partial_fit, cls


def test_discretenb_pickle():
    """Test picklability of discrete naive Bayes classifiers"""

    for cls in [BernoulliNB, MultinomialNB, GaussianNB]:
        clf = cls().fit(X2, y2)
        y_pred = clf.predict(X2)

        store = BytesIO()
        pickle.dump(clf, store)
        clf = pickle.load(BytesIO(store.getvalue()))

        assert_array_equal(y_pred, clf.predict(X2))

        if cls is not GaussianNB:
            # TODO re-enable me when partial_fit is implemented for GaussianNB

            # Test pickling of estimator trained with partial_fit
            clf2 = cls().partial_fit(X2[:3], y2[:3], classes=np.unique(y2))
            clf2.partial_fit(X2[3:], y2[3:])
            store = BytesIO()
            pickle.dump(clf2, store)
            clf2 = pickle.load(BytesIO(store.getvalue()))
            assert_array_equal(y_pred, clf2.predict(X2))


def test_input_check_fit():
    """Test input checks for the fit method"""
    for cls in [BernoulliNB, MultinomialNB, GaussianNB]:
        # check shape consistency for number of samples at fit time
        assert_raises(ValueError, cls().fit, X2, y2[:-1])

        # check shape consistency for number of input features at predict time
        clf = cls().fit(X2, y2)
        assert_raises(ValueError, clf.predict, X2[:, :-1])


def test_input_check_partial_fit():
    for cls in [BernoulliNB, MultinomialNB]:
        # check shape consistency
        assert_raises(ValueError, cls().partial_fit, X2, y2[:-1],
                      classes=np.unique(y2))

        # classes is required for first call to partial fit
        assert_raises(ValueError, cls().partial_fit, X2, y2)

        # check consistency of consecutive classes values
        clf = cls()
        clf.partial_fit(X2, y2, classes=np.unique(y2))
        assert_raises(ValueError, clf.partial_fit, X2, y2,
                      classes=np.arange(42))

        # check consistency of input shape for partial_fit
        assert_raises(ValueError, clf.partial_fit, X2[:, :-1], y2)

        # check consistency of input shape for predict
        assert_raises(ValueError, clf.predict, X2[:, :-1])


def test_discretenb_predict_proba():
    """Test discrete NB classes' probability scores"""

    # The 100s below distinguish Bernoulli from multinomial.
    # FIXME: write a test to show this.
    X_bernoulli = [[1, 100, 0], [0, 1, 0], [0, 100, 1]]
    X_multinomial = [[0, 1], [1, 3], [4, 0]]

    # test binary case (1-d output)
    y = [0, 0, 2]   # 2 is regression test for binary case, 02e673
    for cls, X in zip([BernoulliNB, MultinomialNB],
                      [X_bernoulli, X_multinomial]):
        clf = cls().fit(X, y)
        assert_equal(clf.predict(X[-1]), 2)
        assert_equal(clf.predict_proba(X[0]).shape, (1, 2))
        assert_array_almost_equal(clf.predict_proba(X[:2]).sum(axis=1),
                                  np.array([1., 1.]), 6)

    # test multiclass case (2-d output, must sum to one)
    y = [0, 1, 2]
    for cls, X in zip([BernoulliNB, MultinomialNB],
                      [X_bernoulli, X_multinomial]):
        clf = cls().fit(X, y)
        assert_equal(clf.predict_proba(X[0]).shape, (1, 3))
        assert_equal(clf.predict_proba(X[:2]).shape, (2, 3))
        assert_almost_equal(np.sum(clf.predict_proba(X[1])), 1)
        assert_almost_equal(np.sum(clf.predict_proba(X[-1])), 1)
        assert_almost_equal(np.sum(np.exp(clf.class_log_prior_)), 1)
        assert_almost_equal(np.sum(np.exp(clf.intercept_)), 1)


def test_discretenb_uniform_prior():
    """Test whether discrete NB classes fit a uniform prior
       when fit_prior=False and class_prior=None"""

    for cls in [BernoulliNB, MultinomialNB]:
        clf = cls()
        clf.set_params(fit_prior=False)
        clf.fit([[0], [0], [1]], [0, 0, 1])
        prior = np.exp(clf.class_log_prior_)
        assert_array_equal(prior, np.array([.5, .5]))


def test_discretenb_provide_prior():
    """Test whether discrete NB classes use provided prior"""

    for cls in [BernoulliNB, MultinomialNB]:
        clf = cls(class_prior=[0.5, 0.5])
        clf.fit([[0], [0], [1]], [0, 0, 1])
        prior = np.exp(clf.class_log_prior_)
        assert_array_equal(prior, np.array([.5, .5]))

        # Inconsistent number of classes with prior
        assert_raises(ValueError, clf.fit, [[0], [1], [2]], [0, 1, 2])
        assert_raises(ValueError, clf.partial_fit, [[0], [1]], [0, 1],
                      classes=[0, 1, 1])


def test_deprecated_fit_param():
    warnings.simplefilter("always", DeprecationWarning)
    try:
        for cls in [BernoulliNB, MultinomialNB]:
            clf = cls()
            with warnings.catch_warnings(record=True) as w:
                clf.fit([[0], [1], [2]], [0, 1, 1], class_prior=[0.5, 0.5])

            # Passing class_prior as a fit param should raise a deprecation
            # warning
            assert_equal(len(w), 1)
            assert_equal(w[0].category, DeprecationWarning)

            with warnings.catch_warnings(record=True):
                # Inconsistent number of classes with prior
                assert_raises(ValueError, clf.fit, [[0], [1], [2]], [0, 1, 2],
                              class_prior=[0.5, 0.5])
    finally:
        warnings.filters.pop(0)


def test_sample_weight_multiclass():
    for cls in [BernoulliNB, MultinomialNB]:
        # check shape consistency for number of samples at fit time
        yield check_sample_weight_multiclass, cls


def check_sample_weight_multiclass(cls):
    X = [
        [0, 0, 1],
        [0, 1, 1],
        [0, 1, 1],
        [1, 0, 0],
    ]
    y = [0, 0, 1, 2]
    sample_weight = np.array([1, 1, 2, 2], dtype=np.float)
    sample_weight /= sample_weight.sum()
    clf = cls().fit(X, y, sample_weight=sample_weight)
    assert_array_equal(clf.predict(X), [0, 1, 1, 2])

    # Check wample weight using the partial_fit method
    clf = cls()
    clf.partial_fit(X[:2], y[:2], classes=[0, 1, 2],
                    sample_weight=sample_weight[:2])
    clf.partial_fit(X[2:3], y[2:3], sample_weight=sample_weight[2:3])
    clf.partial_fit(X[3:], y[3:], sample_weight=sample_weight[3:])
    assert_array_equal(clf.predict(X), [0, 1, 1, 2])


def test_sample_weight_mnb():
    clf = MultinomialNB()
    clf.fit([[1, 2], [1, 2], [1, 0]],
            [0, 0, 1],
            sample_weight=[1, 1, 4])
    assert_array_equal(clf.predict([1, 0]), [1])
    positive_prior = np.exp(clf.intercept_[0])
    assert_array_almost_equal([1 - positive_prior, positive_prior],
                              [1 / 3., 2 / 3.])


def test_coef_intercept_shape():
    """coef_ and intercept_ should have shapes as in other linear models.

    Non-regression test for issue #2127.
    """
    X = [[1, 0, 0], [1, 1, 1]]
    y = [1, 2]  # binary classification

    for clf in [MultinomialNB(), BernoulliNB()]:
        clf.fit(X, y)
        assert_equal(clf.coef_.shape, (1, 3))
        assert_equal(clf.intercept_.shape, (1,))


def test_check_accuracy_on_digits():
    # Non regression test to make sure that any further refactoring / optim
    # of the NB models do not harm the performance on a non linearly separable
    # dataset
    digits = load_digits()
    X, y = digits.data, digits.target
    binary_3v8 = np.logical_or(digits.target == 3, digits.target == 8)
    X_3v8, y_3v8 = X[binary_3v8], y[binary_3v8]

    # Multinomial NB
    scores = cross_val_score(MultinomialNB(alpha=10), X, y, cv=10)
    assert_greater(scores.mean(), 0.90)

    scores = cross_val_score(MultinomialNB(alpha=10), X_3v8, y_3v8, cv=10)
    assert_greater(scores.mean(), 0.95)

    # Bernoulli NB
    scores = cross_val_score(BernoulliNB(alpha=10), X > 4, y, cv=10)
    assert_greater(scores.mean(), 0.85)

    scores = cross_val_score(BernoulliNB(alpha=10), X_3v8 > 4, y_3v8, cv=10)
    assert_greater(scores.mean(), 0.94)

    # Gaussian NB
    scores = cross_val_score(GaussianNB(), X, y, cv=10)
    assert_greater(scores.mean(), 0.81)

    scores = cross_val_score(GaussianNB(), X_3v8, y_3v8, cv=10)
    assert_greater(scores.mean(), 0.86)
