import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal

from .. import naive_bayes

# Data is just 6 separable points in the plane
X = np.array( [[-2,-1], [-1, -1], [-1, -2], [1,1], [1,2], [2, 1]])
y = np.array( [1, 1, 1, 2, 2, 2])

def test_gnb():
    """
    Gaussian Naive Bayes classification.

    This checks that GNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    clf =  naive_bayes.GNB()
    y_pred = clf.fit(X, y).predict(X)
    assert_array_equal(y_pred, y)

    y_pred_proba = clf.predict_proba(X)
    y_pred_log_proba = clf.predict_log_proba(X)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)
    
# Data is 6 random points in an 100 dimensional space classified to
# three classes.
X2 = np.random.randint( 5, size=(6, 100) )
y2 = np.array( [1, 1, 2, 2, 3, 3] )

def test_mnnb():
    """
    Multinomial Naive Bayes classification.

    This checks that MNNB implements fit and predict and returns
    correct values for a simple toy dataset.
    """

    #
    # Check the ability to predict the learning set.
    #
    clf =  naive_bayes.MNNB()
    y_pred = clf.fit(X2, y2).predict(X2)

    assert_array_equal(y_pred, y2)
    
    #
    # Verify that np.log(clf.predict_proba(X)) gives the same results as
    # clf.predict_log_proba(X)
    #
    y_pred_proba = clf.predict_proba(X2)
    y_pred_log_proba = clf.predict_log_proba(X2)
    assert_array_almost_equal(np.log(y_pred_proba), y_pred_log_proba, 8)
