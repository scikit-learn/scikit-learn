import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, assert_almost_equal

from sklearn import datasets, ensemble

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
    clf = ensemble.RandomForestClassifier(n_trees=10)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    clf = ensemble.RandomForestClassifier(n_trees=10, max_features=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

def test_classification_toy_et():
    """Check classification on a toy dataset (extra-trees)."""
    clf = ensemble.ExtraTreesClassifier(n_trees=10)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

    clf = ensemble.ExtraTreesClassifier(n_trees=10, max_features=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)

def test_iris_rf():
    """Check consistency on dataset iris (random forest)."""
    for c in ("gini", "entropy"):
        clf = ensemble.RandomForestClassifier(n_trees=10, criterion=c)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

        clf = ensemble.RandomForestClassifier(n_trees=10, criterion=c, max_features=2)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.5, "Failed with criterion %s and score = %f" % (c, score)

def test_iris_et():
    """Check consistency on dataset iris (extra-trees)."""
    for c in ("gini", "entropy"):
        clf = ensemble.ExtraTreesClassifier(n_trees=10, criterion=c)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

        clf = ensemble.ExtraTreesClassifier(n_trees=10, criterion=c, max_features=2)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)

def test_boston_rf():
    """Check consistency on dataset boston house prices (random forest)."""
    for c in ("mse",):
        clf = ensemble.RandomForestRegressor(n_trees=10, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        clf = ensemble.RandomForestRegressor(n_trees=10, criterion=c, max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

def test_boston_et():
    """Check consistency on dataset boston house prices (extra-trees)."""
    for c in ("mse",):
        clf = ensemble.ExtraTreesRegressor(n_trees=10, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)

        clf = ensemble.ExtraTreesRegressor(n_trees=10, criterion=c, max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, "Failed with max_features=None, criterion %s and score = %f" % (c, score)
