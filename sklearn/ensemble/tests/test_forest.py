"""
Testing for the forest module (sklearn.ensemble.forest).
"""

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_almost_equal
from nose.tools import assert_raises

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

def test_classification_toy_rf():
    """Check classification on a toy dataset (random forest)."""
    clf = RandomForestClassifier(n_trees=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    clf = RandomForestClassifier(n_trees=10, max_features=1, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

def test_classification_toy_et():
    """Check classification on a toy dataset (extra-trees)."""
    clf = ExtraTreesClassifier(n_trees=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    clf = ExtraTreesClassifier(n_trees=10, max_features=1, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

def test_iris_rf():
    """Check consistency on dataset iris (random forest)."""
    for c in ("gini", "entropy"):
        clf = RandomForestClassifier(n_trees=10, criterion=c, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

        clf = RandomForestClassifier(n_trees=10, criterion=c, max_features=2, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.5, "Failed with criterion %s and score = %f" % (c, score)

def test_iris_et():
    """Check consistency on dataset iris (extra-trees)."""
    for c in ("gini", "entropy"):
        clf = ExtraTreesClassifier(n_trees=10, criterion=c, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

        clf = ExtraTreesClassifier(n_trees=10, criterion=c, max_features=2, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

def test_boston_rf():
    """Check consistency on dataset boston house prices (random forest)."""
    for c in ("mse",):
        clf = RandomForestRegressor(n_trees=10, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        clf = RandomForestRegressor(n_trees=10, criterion=c, max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

def test_boston_et():
    """Check consistency on dataset boston house prices (extra-trees)."""
    for c in ("mse",):
        clf = ExtraTreesRegressor(n_trees=10, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        clf = ExtraTreesRegressor(n_trees=10, criterion=c, max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

def test_error_rf():
    """Check that proper errors are triggered (random forest)."""
    assert_raises(ValueError, RandomForestClassifier(n_trees=-1).fit, X, y)
    assert_raises(ValueError, RandomForestRegressor(n_trees=-1).fit, X, y)

def test_error_et():
    """Check that proper errors are triggered (random forest)."""
    assert_raises(ValueError, ExtraTreesClassifier(n_trees=-1).fit, X, y)
    assert_raises(ValueError, ExtraTreesRegressor(n_trees=-1).fit, X, y)
