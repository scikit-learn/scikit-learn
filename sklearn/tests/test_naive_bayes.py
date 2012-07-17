import pickle
from io import BytesIO
import numpy as np
import scipy.sparse

from cStringIO import StringIO
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from nose.tools import assert_raises

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
    """
    Multinomial Naive Bayes classification.

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


def test_discretenb_pickle():
    """Test picklability of discrete naive Bayes classifiers"""

    for cls in [BernoulliNB, MultinomialNB, GaussianNB]:
        clf = cls().fit(X2, y2)
        y_pred = clf.predict(X2)

        store = StringIO()
        pickle.dump(clf, store)
        clf = pickle.load(StringIO(store.getvalue()))

        assert_array_equal(y_pred, clf.predict(X2))

    store = BytesIO()
    pickle.dump(clf, store)
    clf = pickle.load(BytesIO(store.getvalue()))


def test_input_check():
    """Test input checks"""
    for cls in [BernoulliNB, MultinomialNB, GaussianNB]:
        clf = cls()
        assert_raises(ValueError, clf.fit, X2, y2[:-1])


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
        clf = cls()
        clf.fit([[0], [0], [1]], [0, 0, 1], class_prior=[0.5, 0.5])
        prior = np.exp(clf.class_log_prior_)
        assert_array_equal(prior, np.array([.5, .5]))


def test_sample_weight():
    clf = MultinomialNB()
    clf.fit([[1, 2], [1, 2], [1, 0]],
            [0, 0, 1],
            sample_weight=[1, 1, 4])
    assert_array_equal(clf.predict([1, 0]), [1])
    positive_prior = np.exp(clf.intercept_)
    assert_array_almost_equal([1 - positive_prior, positive_prior],
                              [1 / 3., 2 / 3.])
