import numpy as np

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_raises

from sklearn.dummy import DummyClassifier
from sklearn.dummy import DummyRegressor


def _check_predict_proba(clf, X, y):
    out = clf.predict_proba(X)
    assert_equal(out.shape[0], len(X))
    assert_equal(out.shape[1], len(np.unique(y)))
    assert_array_equal(out.sum(axis=1), np.ones(len(X)))


def test_most_frequent_strategy():
    X = [[0], [0], [0], [0]]  # ignored
    y = [1, 2, 1, 1]

    clf = DummyClassifier(strategy="most_frequent", random_state=0)
    clf.fit(X, y)
    assert_array_equal(clf.predict(X), np.ones(len(X)))
    _check_predict_proba(clf, X, y)


def test_stratified_strategy():
    X = [[0]] * 5  # ignored
    y = [1, 2, 1, 1, 2]
    clf = DummyClassifier(strategy="stratified", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 1000
    y_pred = clf.predict(X)
    p = np.bincount(y_pred) / float(len(X))
    assert_almost_equal(p[1], 3. / 5, decimal=1)
    assert_almost_equal(p[2], 2. / 5, decimal=1)
    _check_predict_proba(clf, X, y)


def test_uniform_strategy():
    X = [[0]] * 4  # ignored
    y = [1, 2, 1, 1]
    clf = DummyClassifier(strategy="uniform", random_state=0)
    clf.fit(X, y)

    X = [[0]] * 1000
    y_pred = clf.predict(X)
    p = np.bincount(y_pred) / float(len(X))
    assert_almost_equal(p[1], 0.5, decimal=1)
    assert_almost_equal(p[2], 0.5, decimal=1)
    _check_predict_proba(clf, X, y)


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


def test_regressor_exceptions():
    reg = DummyRegressor()
    assert_raises(ValueError, reg.predict, [])
