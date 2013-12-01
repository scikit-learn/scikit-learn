"""
Testing for the forest module (sklearn.ensemble.forest).
"""

# Authors: Gilles Louppe, Brian Holt, Andreas Mueller
# License: BSD 3 clause

from collections import defaultdict

import numpy as np
from numpy.testing import assert_array_equal
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_equal
from numpy.testing import assert_almost_equal
from nose.tools import assert_false, assert_true

from sklearn.utils.testing import assert_less, assert_greater

from sklearn.grid_search import GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomTreesEmbedding
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.svm import LinearSVC
from sklearn.decomposition import TruncatedSVD
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

    # also test apply
    leaf_indices = clf.apply(X)
    assert_equal(leaf_indices.shape, (len(X), clf.n_estimators))

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

    # also test apply
    leaf_indices = clf.apply(X)
    assert_equal(leaf_indices.shape, (len(X), clf.n_estimators))


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


def test_regressor_attributes():
    """Regression models should not have a classes_ attribute."""
    r = RandomForestRegressor()
    assert_false(hasattr(r, "classes_"))
    assert_false(hasattr(r, "n_classes_"))

    r.fit([[1, 2, 3], [4, 5, 6]], [1, 2])
    assert_false(hasattr(r, "classes_"))
    assert_false(hasattr(r, "n_classes_"))


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

    clf = RandomForestClassifier(n_estimators=10)
    clf.fit(X, y)
    importances = clf.feature_importances_
    n_important = sum(importances > 0.1)

    assert_equal(importances.shape[0], 10)
    assert_equal(n_important, 3)

    X_new = clf.transform(X, threshold="mean")
    assert_less(0 < X_new.shape[1], X.shape[1])


def test_oob_score_classification():
    """Check that oob prediction is a good estimation of the generalization
    error."""
    clf = RandomForestClassifier(oob_score=True, random_state=rng)
    n_samples = iris.data.shape[0]
    clf.fit(iris.data[:n_samples / 2, :], iris.target[:n_samples / 2])
    test_score = clf.score(iris.data[n_samples / 2:, :],
                           iris.target[n_samples / 2:])
    assert_less(abs(test_score - clf.oob_score_), 0.1)


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
    clf = ExtraTreesRegressor(random_state=0)
    y_hat = clf.fit(X, y).predict(T)
    assert_almost_equal(y_hat, y_true)
    assert_equal(y_hat.shape, (4, 2))

    np.seterr(**olderr)


def test_classes_shape():
    """Test that n_classes_ and classes_ have proper shape."""
    # Classification, single output
    clf = RandomForestClassifier()
    clf.fit(X, y)

    assert_equal(clf.n_classes_, 2)
    assert_equal(clf.classes_, [-1, 1])

    # Classification, multi-output
    _y = np.vstack((y, np.array(y) * 2)).T
    clf = RandomForestClassifier()
    clf.fit(X, _y)

    assert_equal(len(clf.n_classes_), 2)
    assert_equal(len(clf.classes_), 2)
    assert_equal(clf.n_classes_, [2, 2])
    assert_equal(clf.classes_, [[-1, 1], [-2, 2]])


def test_random_hasher():
    # test random forest hashing on circles dataset
    # make sure that it is linearly separable.
    # even after projected to two SVD dimensions
    # Note: Not all random_states produce perfect results.
    hasher = RandomTreesEmbedding(n_estimators=30, random_state=1)
    X, y = datasets.make_circles(factor=0.5)
    X_transformed = hasher.fit_transform(X)

    # test fit and transform:
    hasher = RandomTreesEmbedding(n_estimators=30, random_state=1)
    assert_array_equal(hasher.fit(X).transform(X).toarray(),
                       X_transformed.toarray())

    # one leaf active per data point per forest
    assert_equal(X_transformed.shape[0], X.shape[0])
    assert_array_equal(X_transformed.sum(axis=1), hasher.n_estimators)
    svd = TruncatedSVD(n_components=2)
    X_reduced = svd.fit_transform(X_transformed)
    linear_clf = LinearSVC()
    linear_clf.fit(X_reduced, y)
    assert_equal(linear_clf.score(X_reduced, y), 1.)


def test_parallel_train():
    rng = np.random.RandomState(12321)
    X = rng.randn(100, 1000)
    y = rng.randint(0, 2, 100)

    clfs = [
        RandomForestClassifier(n_estimators=20,
                               n_jobs=n_jobs,
                               random_state=12345)
        for n_jobs in range(1, 9)
    ]

    for clf in clfs:
        clf.fit(X, y)

    X2 = rng.randn(100, 1000)

    probas = []
    for clf in clfs:
        proba = clf.predict_proba(X2)
        probas.append(proba)

    for proba1, proba2 in zip(probas, probas[1:]):
        assert_true(np.allclose(proba1, proba2))


def test_distribution():
    rng = np.random.RandomState(12321)

    # Single variable with 4 values
    X = rng.randint(0, 4, size=(1000, 1))
    y = rng.rand(1000)
    n_trees = 200

    clf = ExtraTreesRegressor(n_estimators=n_trees, random_state=1).fit(X, y)

    uniques = defaultdict(int)
    for tree in clf.estimators_:
        tree = "".join(("%d,%d/" % (f, int(t)) if f >= 0 else "-")
                       for f, t in zip(tree.tree_.feature,
                                       tree.tree_.threshold))

        uniques[tree] += 1

    uniques = sorted([(1. * count / n_trees, tree)
                      for tree, count in uniques.items()])

    # On a single variable problem where X_0 has 4 equiprobable values, there
    # are 5 ways to build a random tree. The more compact (0,1/0,0/--0,2/--) of
    # them has probability 1/3 while the 4 others have probability 1/6.

    assert_equal(len(uniques), 5)
    assert_greater(0.20, uniques[0][0]) # Rough approximation of 1/6.
    assert_greater(0.20, uniques[1][0])
    assert_greater(0.20, uniques[2][0])
    assert_greater(0.20, uniques[3][0])
    assert_greater(uniques[4][0], 0.3)
    assert_equal(uniques[4][1], "0,1/0,0/--0,2/--")

    # Two variables, one with 2 values, one with 3 values
    X = np.empty((1000, 2))
    X[:, 0] = np.random.randint(0, 2, 1000)
    X[:, 1] = np.random.randint(0, 3, 1000)
    y = rng.rand(1000)

    clf = ExtraTreesRegressor(n_estimators=100,
                              max_features=1,
                              random_state=1).fit(X, y)

    uniques = defaultdict(int)
    for tree in clf.estimators_:
        tree = "".join(("%d,%d/" % (f, int(t)) if f >= 0 else "-")
                       for f, t in zip(tree.tree_.feature,
                                       tree.tree_.threshold))

        uniques[tree] += 1

    uniques = [(count, tree) for tree, count in uniques.items()]
    assert_equal(len(uniques), 8)


if __name__ == "__main__":
    import nose
    nose.runmodule()
