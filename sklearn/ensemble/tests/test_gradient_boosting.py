"""
Testing for the gradient boosting module (sklearn.ensemble.gradient_boosting).
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_equal

from nose.tools import assert_raises

from sklearn.utils import check_random_state
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressor
from sklearn import datasets

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = np.random.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_classification_toy():
    """Check classification on a toy dataset."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    assert_raises(ValueError, clf.predict, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(100, len(clf.estimators_))

    deviance_decrease = (clf.train_deviance[:-1] - clf.train_deviance[1:])
    assert np.any(deviance_decrease >= 0.0), \
           "Train deviance does not monotonically decrease."


def test_parameter_checks():
    """Check input parameter validation."""
    assert_raises(ValueError, GradientBoostingClassifier, n_estimators=0)
    assert_raises(ValueError, GradientBoostingClassifier, n_estimators=-1)

    assert_raises(ValueError, GradientBoostingClassifier, learn_rate=0.0)
    assert_raises(ValueError, GradientBoostingClassifier, learn_rate=-1.0)

    assert_raises(ValueError, GradientBoostingRegressor, loss='foobar')

    assert_raises(ValueError, GradientBoostingClassifier,
                  min_samples_split=0.0)
    assert_raises(ValueError, GradientBoostingClassifier,
                  min_samples_split=-1.0)

    assert_raises(ValueError, GradientBoostingClassifier, subsample=0.0)
    assert_raises(ValueError, GradientBoostingClassifier, subsample=1.1)
    assert_raises(ValueError, GradientBoostingClassifier, subsample=-0.1)

    assert_raises(ValueError, GradientBoostingClassifier, max_depth=-0.1)
    assert_raises(ValueError, GradientBoostingClassifier, max_depth=0)

    assert_raises(ValueError, GradientBoostingClassifier, init={})

    # test fit before feature importance
    assert_raises(ValueError,
                  lambda: GradientBoostingClassifier().feature_importances_)

    # test value error on multi-class
    assert_raises(ValueError,
                  lambda X, y: GradientBoostingClassifier().fit(X, y),
                  X, [0, 0, 1, 1, 2, 2])


def test_classification_synthetic():
    """Test GradientBoostingClassifier on synthetic dataset used by
    Hastie et al. in ESLII Example 12.7. """
    rs = check_random_state(1)
    shape = (12000, 10)
    X = rs.normal(size=shape).reshape(shape)
    y = ((X ** 2.0).sum(axis=1) > 9.34).astype(np.float64)
    y[y == 0.0] = -1.0

    X_train, X_test = X[:2000], X[2000:]
    y_train, y_test = y[:2000], y[2000:]

    gbrt = GradientBoostingClassifier(n_estimators=100, min_samples_split=1,
                                      max_depth=1,
                                      learn_rate=1.0, random_state=0)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test))
    assert error_rate < 0.085, \
           "GB failed with error %.4f" % error_rate

    gbrt = GradientBoostingClassifier(n_estimators=200, min_samples_split=1,
                                      max_depth=1,
                                      learn_rate=1.0, subsample=0.5,
                                      random_state=0)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test))
    assert error_rate < 0.08, \
           "Stochastic GB failed with error %.4f" % error_rate


def test_boston():
    """Check consistency on dataset boston house prices with least squares
    and least absolute deviation. """
    for loss in ("ls", "lad"):
        clf = GradientBoostingRegressor(n_estimators=100, loss=loss,
                                        max_depth=4,
                                        min_samples_split=1, random_state=1)
        assert_raises(ValueError, clf.predict, boston.data)
        clf.fit(boston.data, boston.target)
        y_pred = clf.predict(boston.data)
        mse = np.mean((y_pred - boston.target) ** 2.0)
        assert mse < 6.0, "Failed with loss %s and mse = %.4f" % (loss, mse)


def test_regression_synthetic():
    """Test on synthetic regression datasets used in Leo Breiman,
    `Bagging Predictors?. Machine Learning 24(2): 123-140 (1996). """
    random_state = check_random_state(1)
    regression_params = {'n_estimators': 100, 'max_depth': 4,
                         'min_samples_split': 1, 'learn_rate': 0.1,
                         'loss': 'ls'}

    # Friedman1
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=random_state, noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor()
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    assert mse < 5.0, "Failed on Friedman1 with mse = %.4f" % mse

    # Friedman2
    X, y = datasets.make_friedman2(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    assert mse < 1700.0, "Failed on Friedman2 with mse = %.4f" % mse

    # Friedman3
    X, y = datasets.make_friedman3(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    assert mse < 0.015, "Failed on Friedman3 with mse = %.4f" % mse


def test_feature_importances():
    clf = GradientBoostingRegressor(n_estimators=100, max_depth=4,
                                    min_samples_split=1, random_state=1)
    clf.fit(boston.data, boston.target)
    feature_importances = clf.feature_importances_

    # true feature importance ranking
    true_ranking = np.array([3,  1,  8, 10,  2,  9,  4, 11,  0,  6,  7,  5, 12])

    assert_array_equal(true_ranking, feature_importances.argsort())


def test_probability():
    """Predict probabilities."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)

    assert_raises(ValueError, clf.predict_proba, T)

    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    # check if probabilities are in [0, 1].
    y_proba = clf.predict_proba(T)
    assert np.all(y_proba >= 0.0)
    assert np.all(y_proba <= 1.0)

    # derive predictions from probabilities
    y_pred = clf.classes.take(y_proba.argmax(axis=1), axis=0)
    assert_array_equal(y_pred, true_result)


def test_check_inputs():
    """Test input checks (shape and type of X and y)."""
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    assert_raises(ValueError, clf.fit, X, y + [0, 1])

    from scipy import sparse
    X_sparse = sparse.csr_matrix(X)
    clf = GradientBoostingClassifier(n_estimators=100, random_state=1)
    assert_raises(ValueError, clf.fit, X_sparse, y)

    clf = GradientBoostingClassifier().fit(X, y)
    assert_raises(ValueError, clf.predict, X_sparse)
