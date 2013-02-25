import warnings
import numpy as np

from sklearn.base import clone
from sklearn.externals.six.moves import xrange
from sklearn.utils.testing import (assert_array_equal,
                                   assert_equal,
                                   assert_almost_equal,
                                   assert_raises)

from sklearn.dummy import DummyClassifier, DummyRegressor


def _check_predict_proba(clf, X, y):
    proba = clf.predict_proba(X)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # We know that we can have division by zero
        log_proba = clf.predict_log_proba(X)

    y = np.atleast_1d(y)
    if y.ndim == 1:
        y = np.reshape(y, (-1, 1))

    n_outputs = y.shape[1]
    n_samples = len(X)

    if n_outputs == 1:
        proba = [proba]
        log_proba = [log_proba]

    for k in xrange(n_outputs):
        assert_equal(proba[k].shape[0], n_samples)
        assert_equal(proba[k].shape[1], len(np.unique(y[:, k])))
        assert_array_equal(proba[k].sum(axis=1), np.ones(len(X)))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # We know that we can have division by zero
            assert_array_equal(np.log(proba[k]), log_proba[k])


def _check_behavior_2d(clf):
    # 1d case
    X = np.array([[0], [0], [0], [0]])  # ignored
    y = np.array([1, 2, 1, 1])
    est = clone(clf)
    est.fit(X, y)
    y_pred = est.predict(X)
    assert_equal(y.shape, y_pred.shape)

    # 2d case
    y = np.array([[1, 0],
                  [2, 0],
                  [1, 0],
                  [1, 3]])
    est = clone(clf)
    est.fit(X, y)
    y_pred = est.predict(X)
    assert_equal(y.shape, y_pred.shape)


def test_most_frequent_strategy():
    X = [[0], [0], [0], [0]]  # ignored
    y = [1, 2, 1, 1]

    clf = DummyClassifier(strategy="most_frequent", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), np.ones(len(X)))
    _check_predict_proba(clf, X, y)


def test_most_frequent_strategy_multioutput():
    X = [[0], [0], [0], [0]]  # ignored
    y = np.array([[1, 0],
                  [2, 0],
                  [1, 0],
                  [1, 3]])

    n_samples = len(X)

    clf = DummyClassifier(strategy="most_frequent", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X),
                       np.hstack([np.ones((n_samples, 1)),
                                  np.zeros((n_samples, 1))]))
    _check_predict_proba(clf, X, y)
    _check_behavior_2d(clf)


def test_stratified_strategy():
    X = [[0]] * 5  # ignored
    y = [1, 2, 1, 1, 2]
    clf = DummyClassifier(strategy="stratified", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)
    p = np.bincount(y_pred) / float(len(X))
    assert_almost_equal(p[1], 3. / 5, decimal=1)
    assert_almost_equal(p[2], 2. / 5, decimal=1)
    _check_predict_proba(clf, X, y)


def test_stratified_strategy_multioutput():
    X = [[0]] * 5  # ignored
    y = np.array([[2, 1],
                  [2, 2],
                  [1, 1],
                  [1, 2],
                  [1, 1]])

    clf = DummyClassifier(strategy="stratified", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)

    for k in xrange(y.shape[1]):
        p = np.bincount(y_pred[:, k]) / float(len(X))
        assert_almost_equal(p[1], 3. / 5, decimal=1)
        assert_almost_equal(p[2], 2. / 5, decimal=1)
        _check_predict_proba(clf, X, y)

    _check_behavior_2d(clf)


def test_uniform_strategy():
    X = [[0]] * 4  # ignored
    y = [1, 2, 1, 1]
    clf = DummyClassifier(strategy="uniform", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)
    p = np.bincount(y_pred) / float(len(X))
    assert_almost_equal(p[1], 0.5, decimal=1)
    assert_almost_equal(p[2], 0.5, decimal=1)
    _check_predict_proba(clf, X, y)


def test_uniform_strategy_multioutput():
    X = [[0]] * 4  # ignored
    y = np.array([[2, 1],
                  [2, 2],
                  [1, 2],
                  [1, 1]])
    clf = DummyClassifier(strategy="uniform", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 500
    y_pred = clf.predict(X)

    for k in xrange(y.shape[1]):
        p = np.bincount(y_pred[:, k]) / float(len(X))
        assert_almost_equal(p[1], 0.5, decimal=1)
        assert_almost_equal(p[2], 0.5, decimal=1)
        _check_predict_proba(clf, X, y)

    _check_behavior_2d(clf)


def test_string_labels():
    X = [[0]] * 5
    y = ["paris", "paris", "tokyo", "amsterdam", "berlin"]
    clf = DummyClassifier(strategy="most_frequent")
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), ["paris"] * 5)


def test_classifier_exceptions():
    clf = DummyClassifier(strategy="unknown")
    assert_raises(ValueError, clf.fit, [], [])

    assert_raises(ValueError, clf.predict, [])
    assert_raises(ValueError, clf.predict_proba, [])


def test_regressor():
    X = [[0]] * 4  # ignored
    y = [1, 2, 1, 1]

    reg = DummyRegressor()
    reg.fit(X, y)
    assert_array_equal(reg.predict(X), [5. / 4] * len(X))


def test_multioutput_regressor():

    X_learn = np.random.randn(10, 10)
    y_learn = np.random.randn(10, 5)

    mean = np.mean(y_learn, axis=0).reshape((1, -1))

    X_test = np.random.randn(20, 10)
    y_test = np.random.randn(20, 5)

    # Correctness oracle
    est = DummyRegressor()
    est.fit(X_learn, y_learn)
    y_pred_learn = est.predict(X_learn)
    y_pred_test = est.predict(X_test)

    assert_array_equal(np.tile(mean, (y_learn.shape[0], 1)), y_pred_learn)
    assert_array_equal(np.tile(mean, (y_test.shape[0], 1)), y_pred_test)
    _check_behavior_2d(est)


def test_regressor_exceptions():
    reg = DummyRegressor()
    assert_raises(ValueError, reg.predict, [])
