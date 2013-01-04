import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import raises
from sklearn.utils.testing import assert_greater

from sklearn.linear_model import logistic
from sklearn import datasets

X = [[-1, 0], [0, 1], [1, 1]]
X_sp = sp.csr_matrix(X)
Y1 = [0, 1, 1]
Y2 = [2, 1, 0]
iris = datasets.load_iris()


def check_predictions(clf, X, y):
    """Check that the model is able to fit the classification data"""
    n_samples = len(y)
    classes = np.unique(y)
    n_classes = classes.shape[0]

    predicted = clf.fit(X, y).predict(X)
    assert_array_equal(clf.classes_, classes)

    assert_equal(predicted.shape, (n_samples,))
    assert_array_equal(predicted, y)

    probabilities = clf.predict_proba(X)
    assert_equal(probabilities.shape, (n_samples, n_classes))
    assert_array_almost_equal(probabilities.sum(axis=1), np.ones(n_samples))
    assert_array_equal(probabilities.argmax(axis=1), y)


def test_predict_2_classes():
    """Simple sanity check on a 2 classes dataset

    Make sure it predicts the correct result on simple datasets.
    """
    check_predictions(logistic.LogisticRegression(), X, Y1)
    check_predictions(logistic.LogisticRegression(), X_sp, Y1)

    check_predictions(logistic.LogisticRegression(C=100), X, Y1)
    check_predictions(logistic.LogisticRegression(C=100), X_sp, Y1)

    check_predictions(logistic.LogisticRegression(fit_intercept=False),
                      X, Y1)
    check_predictions(logistic.LogisticRegression(fit_intercept=False),
                      X_sp, Y1)


def test_error():
    """Test for appropriate exception on errors"""
    assert_raises(ValueError, logistic.LogisticRegression(C=-1).fit, X, Y1)


def test_predict_3_classes():
    check_predictions(logistic.LogisticRegression(C=10), X, Y2)
    check_predictions(logistic.LogisticRegression(C=10), X_sp, Y2)


def test_predict_iris():
    """Test logisic regression with the iris dataset"""
    n_samples, n_features = iris.data.shape

    target = iris.target_names[iris.target]
    clf = logistic.LogisticRegression(C=len(iris.data)).fit(iris.data, target)
    assert_array_equal(np.unique(target), clf.classes_)

    pred = clf.predict(iris.data)
    assert_greater(np.mean(pred == target), .95)

    probabilities = clf.predict_proba(iris.data)
    assert_array_almost_equal(probabilities.sum(axis=1), np.ones(n_samples))

    pred = iris.target_names[probabilities.argmax(axis=1)]
    assert_greater(np.mean(pred == target), .95)


def test_inconsistent_input():
    """Test that an exception is raised on inconsistent input"""
    rng = np.random.RandomState(0)
    X_ = rng.random_sample((5, 10))
    y_ = np.ones(X_.shape[0])
    y_[0] = 0

    clf = logistic.LogisticRegression()

    # Wrong dimensions for training data
    y_wrong = y_[:-1]
    assert_raises(ValueError, clf.fit, X, y_wrong)

    # Wrong dimensions for test data
    assert_raises(ValueError, clf.fit(X_, y_).predict,
                  rng.random_sample((3, 12)))


def test_write_parameters():
    """Test that we can write to coef_ and intercept_"""
    #rng = np.random.RandomState(0)
    #X = rng.random_sample((5, 10))
    #y = np.ones(X.shape[0])
    clf = logistic.LogisticRegression()
    clf.fit(X, Y1)
    clf.coef_[:] = 0
    clf.intercept_[:] = 0
    assert_array_equal(clf.decision_function(X), 0)


@raises(ValueError)
def test_nan():
    """Test proper NaN handling.

    Regression test for Issue #252: fit used to go into an infinite loop.
    """
    Xnan = np.array(X, dtype=np.float64)
    Xnan[0, 1] = np.nan
    logistic.LogisticRegression().fit(Xnan, Y1)


def test_liblinear_random_state():
    X, y = datasets.make_classification(n_samples=20)
    lr1 = logistic.LogisticRegression(random_state=0)
    lr1.fit(X, y)
    lr2 = logistic.LogisticRegression(random_state=0)
    lr2.fit(X, y)
    assert_array_almost_equal(lr1.coef_, lr2.coef_)
