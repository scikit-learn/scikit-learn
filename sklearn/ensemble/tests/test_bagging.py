"""
Testing for the forest module (sklearn.ensemble.forest).
"""

# Authors: Gilles Louppe, Brian Holt
# License: BSD 3

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal

from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import BaggedClassifier
from sklearn.ensemble import BaggedRegressor
from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
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
    clf = BaggedClassifier(base_estimator=DecisionTreeClassifier(), n_estimators=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))


def test_iris():
    """Check consistency on dataset iris."""
    clf = BaggedClassifier(base_estimator=DecisionTreeClassifier(), n_estimators=10, random_state=1)
    clf.fit(iris.data, iris.target)
    score = clf.score(iris.data, iris.target)
    assert score > 0.9, "Failed with criterion %s and score = %f" % (c, score)


def test_boston():
    """Check consistency on dataset boston house prices."""
    clf = BaggedRegressor(base_estimator=DecisionTreeRegressor(), n_estimators=10, random_state=1)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert score < 3, ("Failed with max_features=None, "
                       "criterion %s and score = %f" % (c, score))


def test_oob_score_classification():
    """Check that oob prediction is as acurate as
    usual prediction on the training set.
    Not really a good test that prediction is independent."""
    clf = BaggedClassifier(base_estimator=DecisionTreeClassifier(), bootstrap=True, oob_score=True, n_estimators=10, random_state=1)
    clf.fit(X, y)
    training_score = clf.score(X, y)
    assert_almost_equal(training_score, clf.oob_score_)


def test_oob_score_regression():
    """Check that oob prediction is pessimistic estimate.
    Not really a good test that prediction is independent."""
    clf = BaggedRegressor(base_estimator=DecisionTreeRegressor(), bootstrap=True, oob_score=True, n_estimators=30, random_state=1)
    n_samples = boston.data.shape[0]
    clf.fit(boston.data[:n_samples / 2, :], boston.target[:n_samples / 2])
    test_score = clf.score(boston.data[n_samples / 2:, :],
            boston.target[n_samples / 2:])
    assert(test_score > clf.oob_score_)
    assert(clf.oob_score_ > .8)
