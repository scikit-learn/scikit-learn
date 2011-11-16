"""
Testing for the forest module (sklearn.ensemble.forest).
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from nose.tools import assert_raises

from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn import datasets

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
np.random.seed([1])
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = np.random.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_classification_toy():
    """Check classification on a toy dataset."""
    # Random forest
    clf = RandomForestClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    clf = RandomForestClassifier(n_estimators=10, max_features=1, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    # Extra-trees
    clf = ExtraTreesClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    clf = ExtraTreesClassifier(n_estimators=10, max_features=1, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)


def test_iris():
    """Check consistency on dataset iris."""
    for c in ("gini", "entropy"):
        # Random forest
        clf = RandomForestClassifier(n_estimators=10, criterion=c, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

        clf = RandomForestClassifier(n_estimators=10, criterion=c, max_features=2, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.5, "Failed with criterion %s and score = %f" % (c, score)

        # Extra-trees
        clf = ExtraTreesClassifier(n_estimators=10, criterion=c, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

        clf = ExtraTreesClassifier(n_estimators=10, criterion=c, max_features=2, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)


def test_boston():
    """Check consistency on dataset boston house prices."""
    for c in ("mse",):
        # Random forest
        clf = RandomForestRegressor(n_estimators=10, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        clf = RandomForestRegressor(n_estimators=10, criterion=c, max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        # Extra-trees
        clf = ExtraTreesRegressor(n_estimators=10, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        clf = ExtraTreesRegressor(n_estimators=10, criterion=c, max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)


def test_probability():
    """Predict probabilities."""
    # Random forest
    clf = RandomForestClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)
    assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1), np.ones(iris.data.shape[0]))
    assert_array_almost_equal(clf.predict_proba(iris.data), np.exp(clf.predict_log_proba(iris.data)))

    # Extra-trees
    clf = ExtraTreesClassifier(n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)
    assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1), np.ones(iris.data.shape[0]))
    assert_array_almost_equal(clf.predict_proba(iris.data), np.exp(clf.predict_log_proba(iris.data)))


def test_gridsearch():
    """Check that base trees can be grid-searched."""
    # Random forest
    forest = RandomForestClassifier()
    parameters = {'n_estimators': (1,2),
                  'base_estimator__max_depth': (1,2)}
    clf = GridSearchCV(forest, parameters)
    clf.fit(iris.data, iris.target)

    # Extra-trees
    forest = ExtraTreesClassifier()
    parameters = {'n_estimators': (1,2),
                  'base_estimator__max_depth': (1,2)}
    clf = GridSearchCV(forest, parameters)
    clf.fit(iris.data, iris.target)


def test_error():
    """Check that proper errors are triggered."""
    def instantiate(class_name, **params):
        return class_name(**params)

    # Random forest
    assert_raises(ValueError, instantiate, class_name=RandomForestClassifier, n_estimators=-1)
    assert_raises(ValueError, instantiate, class_name=RandomForestRegressor, n_estimators=-1)

    # Extra-trees
    assert_raises(ValueError, instantiate, class_name=ExtraTreesClassifier, n_estimators=-1)
    assert_raises(ValueError, instantiate, class_name=ExtraTreesRegressor, n_estimators=-1)


def test_pickle():
    """Check pickability."""
    import pickle

    # Random forest
    obj = RandomForestClassifier()
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(iris.data, iris.target)
    assert score == score2

    obj = RandomForestRegressor()
    obj.fit(boston.data, boston.target)
    score = obj.score(boston.data, boston.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(boston.data, boston.target)
    assert score == score2

    # Extra-trees
    obj = ExtraTreesClassifier()
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(iris.data, iris.target)
    assert score == score2

    obj = ExtraTreesRegressor()
    obj.fit(boston.data, boston.target)
    score = obj.score(boston.data, boston.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(boston.data, boston.target)
    assert score == score2


if __name__ == "__main__":
    import nose
    nose.runmodule()
