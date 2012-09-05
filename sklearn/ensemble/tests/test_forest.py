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
from nose.tools import assert_true

from sklearn.utils.testing import assert_less, assert_greater

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
rng = np.random.RandomState(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]


def test_classification_toy():
    """Check classification on a toy dataset."""
    # Random forest
    clf = RandomForestClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))

    clf = RandomForestClassifier(n_estimators=10, max_features=1,
                                 random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))

    # Extra-trees
    clf = ExtraTreesClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))

    clf = ExtraTreesClassifier(n_estimators=10, max_features=1,
                               random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))


def test_iris():
    """Check consistency on dataset iris."""
    for c in ("gini", "entropy"):
        # Random forest
        clf = RandomForestClassifier(n_estimators=10, criterion=c,
                                     random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c,
                                                                         score)

        clf = RandomForestClassifier(n_estimators=10, criterion=c,
                                     max_features=2, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.5, "Failed with criterion %s and score = %f" % (c,
                                                                         score)

        # Extra-trees
        clf = ExtraTreesClassifier(n_estimators=10, criterion=c,
                                   random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c,
                                                                         score)

        clf = ExtraTreesClassifier(n_estimators=10, criterion=c,
                                   max_features=2, random_state=1)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)
        assert score > 0.9, "Failed with criterion %s and score = %f" % (c,
                                                                         score)


def test_boston():
    """Check consistency on dataset boston house prices."""
    for c in ("mse",):
        # Random forest
        clf = RandomForestRegressor(n_estimators=5, criterion=c,
                                    random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, ("Failed with max_features=None, "
                           "criterion %s and score = %f" % (c, score))

        clf = RandomForestRegressor(n_estimators=5, criterion=c,
                                    max_features=6, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, ("Failed with max_features=None, "
                           "criterion %s and score = %f" % (c, score))

        # Extra-trees
        clf = ExtraTreesRegressor(n_estimators=5, criterion=c, random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, ("Failed with max_features=None, "
                           "criterion %s and score = %f" % (c, score))

        clf = ExtraTreesRegressor(n_estimators=5, criterion=c, max_features=6,
                                  random_state=1)
        clf.fit(boston.data, boston.target)
        score = clf.score(boston.data, boston.target)
        assert score < 3, ("Failed with max_features=None, "
                           "criterion %s and score = %f" % (c, score))


def test_probability():
    """Predict probabilities."""
    olderr = np.seterr(divide="ignore")

    # Random forest
    clf = RandomForestClassifier(n_estimators=10, random_state=1,
            max_features=1, max_depth=1)
    clf.fit(iris.data, iris.target)
    assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1),
                              np.ones(iris.data.shape[0]))
    assert_array_almost_equal(clf.predict_proba(iris.data),
                              np.exp(clf.predict_log_proba(iris.data)))

    # Extra-trees
    clf = ExtraTreesClassifier(n_estimators=10, random_state=1, max_features=1,
            max_depth=1)
    clf.fit(iris.data, iris.target)
    assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1),
                              np.ones(iris.data.shape[0]))
    assert_array_almost_equal(clf.predict_proba(iris.data),
                              np.exp(clf.predict_log_proba(iris.data)))

    np.seterr(**olderr)


def test_importances():
    """Check variable importances."""
    X, y = datasets.make_classification(n_samples=1000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=0)

    clf = RandomForestClassifier(n_estimators=10, compute_importances=True)
    clf.fit(X, y)
    importances = clf.feature_importances_
    n_important = sum(importances > 0.1)

    assert_equal(importances.shape[0], 10)
    assert_equal(n_important, 3)

    X_new = clf.transform(X, threshold="mean")
    assert_less(0 < X_new.shape[1], X.shape[1])

    clf = RandomForestClassifier(n_estimators=10)
    clf.fit(X, y)
    assert_true(clf.feature_importances_ is None)


def test_oob_score_classification():
    """Check that oob prediction is as acurate as
    usual prediction on the training set.
    Not really a good test that prediction is independent."""
    clf = RandomForestClassifier(oob_score=True, random_state=rng)
    clf.fit(X, y)
    training_score = clf.score(X, y)
    assert_almost_equal(training_score, clf.oob_score_)


def test_oob_score_regression():
    """Check that oob prediction is pessimistic estimate.
    Not really a good test that prediction is independent."""
    clf = RandomForestRegressor(n_estimators=50, oob_score=True,
            random_state=rng)
    n_samples = boston.data.shape[0]
    clf.fit(boston.data[:n_samples / 2, :], boston.target[:n_samples / 2])
    test_score = clf.score(boston.data[n_samples / 2:, :],
                           boston.target[n_samples / 2:])
    assert_greater(test_score, clf.oob_score_)
    assert_greater(clf.oob_score_, .8)


def test_gridsearch():
    """Check that base trees can be grid-searched."""
    # Random forest
    forest = RandomForestClassifier()
    parameters = {'n_estimators': (1, 2),
                  'max_depth': (1, 2)}
    clf = GridSearchCV(forest, parameters)
    clf.fit(iris.data, iris.target)

    # Extra-trees
    forest = ExtraTreesClassifier()
    parameters = {'n_estimators': (1, 2),
                  'max_depth': (1, 2)}
    clf = GridSearchCV(forest, parameters)
    clf.fit(iris.data, iris.target)


def test_parallel():
    """Check parallel computations."""
    # Classification
    forest = RandomForestClassifier(n_estimators=10, n_jobs=3, random_state=0)

    forest.fit(iris.data, iris.target)
    assert_true(10 == len(forest))

    forest.set_params(n_jobs=1)
    y1 = forest.predict(iris.data)
    forest.set_params(n_jobs=2)
    y2 = forest.predict(iris.data)
    assert_array_equal(y1, y2)

    # Regression
    forest = RandomForestRegressor(n_estimators=10, n_jobs=3, random_state=0)

    forest.fit(boston.data, boston.target)
    assert_true(10 == len(forest))

    forest.set_params(n_jobs=1)
    y1 = forest.predict(boston.data)
    forest.set_params(n_jobs=2)
    y2 = forest.predict(boston.data)
    assert_array_almost_equal(y1, y2, 3)

    # Use all cores on the classification dataset
    forest = RandomForestClassifier(n_jobs=-1)
    forest.fit(iris.data, iris.target)


def test_pickle():
    """Check pickability."""
    import pickle

    # Random forest
    obj = RandomForestClassifier(random_state=0)
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(iris.data, iris.target)
    assert_true(score == score2)

    obj = RandomForestRegressor(random_state=0)
    obj.fit(boston.data, boston.target)
    score = obj.score(boston.data, boston.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(boston.data, boston.target)
    assert_true(score == score2)

    # Extra-trees
    obj = ExtraTreesClassifier(random_state=0)
    obj.fit(iris.data, iris.target)
    score = obj.score(iris.data, iris.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(iris.data, iris.target)
    assert_true(score == score2)

    obj = ExtraTreesRegressor(random_state=0)
    obj.fit(boston.data, boston.target)
    score = obj.score(boston.data, boston.target)
    s = pickle.dumps(obj)

    obj2 = pickle.loads(s)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(boston.data, boston.target)
    assert_true(score == score2)


def test_multioutput():
    """Check estimators on multi-output problems."""
    olderr = np.seterr(divide="ignore")

    X = [[-2, -1],
         [-1, -1],
         [-1, -2],
         [1, 1],
         [1, 2],
         [2, 1],
         [-2, 1],
         [-1, 1],
         [-1, 2],
         [2, -1],
         [1, -1],
         [1, -2]]

    y = [[-1, 0],
         [-1, 0],
         [-1, 0],
         [1, 1],
         [1, 1],
         [1, 1],
         [-1, 2],
         [-1, 2],
         [-1, 2],
         [1, 3],
         [1, 3],
         [1, 3]]

    T = [[-1, -1], [1, 1], [-1, 1], [1, -1]]
    y_true = [[-1, 0], [1, 1], [-1, 2], [1, 3]]

    # toy classification problem
    clf = ExtraTreesClassifier(random_state=0)
    y_hat = clf.fit(X, y).predict(T)
    assert_array_equal(y_hat, y_true)
    assert_equal(y_hat.shape, (4, 2))

    proba = clf.predict_proba(T)
    assert_equal(len(proba), 2)
    assert_equal(proba[0].shape, (4, 2))
    assert_equal(proba[1].shape, (4, 4))

    log_proba = clf.predict_log_proba(T)
    assert_equal(len(log_proba), 2)
    assert_equal(log_proba[0].shape, (4, 2))
    assert_equal(log_proba[1].shape, (4, 4))

    # toy regression problem
    clf = ExtraTreesRegressor(random_state=5)
    y_hat = clf.fit(X, y).predict(T)
    assert_almost_equal(y_hat, y_true)
    assert_equal(y_hat.shape, (4, 2))

    np.seterr(**olderr)


if __name__ == "__main__":
    import nose
    nose.runmodule()
