import cPickle as pickle
from cStringIO import StringIO
import numpy as np
import scipy.sparse
from numpy.testing import assert_array_equal, assert_array_almost_equal

from .. import naive_bayes

# Data is just 6 separable points in the plane
X = np.array([[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]])
y = np.array([1, 1, 1, 2, 2, 2])


def test_gnb():
    """
    Gaussian Naive Bayes classification.

    This checks that GaussianNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf = naive_bayes.GaussianNB()
    y_pred = clf.fit(X, y).predict(X)
    assert_array_equal(y_pred, y)

    y_pred_proba = clf.predict_proba(X)
    y_pred_log_proba = clf.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)


# Data is 6 random points in an 100 dimensional space classified to
# three classes.
X2 = np.random.randint(5, size=(6, 100))
y2 = np.array([1, 1, 2, 2, 3, 3])


def test_mnnb():
    """
    Multinomial Naive Bayes classification.

    This checks that MultinomialNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    #
    # Check the ability to predict the learning set.
    #
    clf = naive_bayes.MultinomialNB()
    y_pred = clf.fit(X2, y2).predict(X2)

    assert_array_equal(y_pred, y2)

    #
    # Verify that np.log(clf.predict_proba(X)) gives the same results as
    # clf.predict_log_proba(X)
    #
    y_pred_proba = clf.predict_proba(X2)
    y_pred_log_proba = clf.predict_log_proba(X2)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)


def test_sparse_mnnb():
    """
    Multinomial Naive Bayes classification for sparse data.

    This checks that MultinomialNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    X2S = scipy.sparse.csr_matrix(X2)

    #
    # Check the ability to predict the learning set.
    #
    clf = naive_bayes.MultinomialNB()
    y_pred = clf.fit(X2S, y2).predict(X2S)

    assert_array_equal(y_pred, y2)

    #
    # Verify that np.log(clf.predict_proba(X)) gives the same results as
    # clf.predict_log_proba(X)
    #
    y_pred_proba = clf.predict_proba(X2S)
    y_pred_log_proba = clf.predict_log_proba(X2S)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)


def test_mnnb_pickle():
    '''Test picklability of multinomial NB'''

    clf = naive_bayes.MultinomialNB(alpha=2, fit_prior=False).fit(X, y)
    y_pred = clf.predict(X)

    store = StringIO()
    pickle.dump(clf, store)
    clf = pickle.load(StringIO(store.getvalue()))

    assert_array_equal(y_pred, clf.predict(X))


def test_mnnb_predict_proba():
    '''Test multinomial NB's probability scores'''

    clf = naive_bayes.MultinomialNB().fit([[0,1], [0,1], [1,0]], [0,0,1])
    assert clf.predict([0,1]) == 0
    assert np.sum(clf.predict_proba([0,1])) == 1
    assert np.sum(clf.predict_proba([1,0])) == 1
