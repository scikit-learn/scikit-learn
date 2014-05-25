"""
Testing for the forest module (sklearn.ensemble.forest).
"""

# Authors: Gilles Louppe,
#          Brian Holt,
#          Andreas Mueller,
#          Arnaud Joly
# License: BSD 3 clause

import pickle
from collections import defaultdict
from itertools import product

import numpy as np
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_false, assert_true
from sklearn.utils.testing import assert_less, assert_greater

from sklearn import datasets
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomTreesEmbedding
from sklearn.grid_search import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.utils.validation import check_random_state

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = check_random_state(0)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

FOREST_CLASSIFIERS = {
    "ExtraTreesClassifier": ExtraTreesClassifier,
    "RandomForestClassifier": RandomForestClassifier,
}

FOREST_REGRESSORS = {
    "ExtraTreesRegressor": ExtraTreesRegressor,
    "RandomForestRegressor": RandomForestRegressor,
}

FOREST_TRANSFORMERS = {
    "RandomTreesEmbedding": RandomTreesEmbedding,
}

FOREST_ESTIMATORS = dict()
FOREST_ESTIMATORS.update(FOREST_CLASSIFIERS)
FOREST_ESTIMATORS.update(FOREST_REGRESSORS)
FOREST_ESTIMATORS.update(FOREST_TRANSFORMERS)


def check_classification_toy(name):
    """Check classification on a toy dataset."""
    ForestClassifier = FOREST_CLASSIFIERS[name]

    clf = ForestClassifier(n_estimators=10, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))

    clf = ForestClassifier(n_estimators=10, max_features=1, random_state=1)
    clf.fit(X, y)
    assert_array_equal(clf.predict(T), true_result)
    assert_equal(10, len(clf))

    # also test apply
    leaf_indices = clf.apply(X)
    assert_equal(leaf_indices.shape, (len(X), clf.n_estimators))


def test_classification_toy():
    for name in FOREST_CLASSIFIERS:
        yield check_classification_toy, name


def check_iris_criterion(name, criterion):
    """Check consistency on dataset iris."""
    ForestClassifier = FOREST_CLASSIFIERS[name]

    clf = ForestClassifier(n_estimators=10, criterion=criterion,
                           random_state=1)
    clf.fit(iris.data, iris.target)
    score = clf.score(iris.data, iris.target)
    assert_greater(score, 0.9, "Failed with criterion %s and score = %f"
                               % (criterion, score))

    clf = ForestClassifier(n_estimators=10, criterion=criterion,
                           max_features=2, random_state=1)
    clf.fit(iris.data, iris.target)
    score = clf.score(iris.data, iris.target)
    assert_greater(score, 0.5, "Failed with criterion %s and score = %f"
                               % (criterion, score))


def test_iris():
    for name, criterion in product(FOREST_CLASSIFIERS, ("gini", "entropy")):
        yield check_iris_criterion, name, criterion


def check_boston_criterion(name, criterion):
    """Check consistency on dataset boston house prices."""
    ForestRegressor = FOREST_REGRESSORS[name]

    clf = ForestRegressor(n_estimators=5, criterion=criterion, random_state=1)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert_greater(score, 0.95, "Failed with max_features=None, criterion %s "
                                "and score = %f" % (criterion, score))

    clf = ForestRegressor(n_estimators=5, criterion=criterion,
                          max_features=6, random_state=1)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert_greater(score, 0.95, "Failed with max_features=6, criterion %s "
                                "and score = %f" % (criterion, score))


def test_boston():
    for name, criterion in product(FOREST_REGRESSORS, ("mse", )):
        yield check_boston_criterion, name, criterion


def check_regressor_attributes(name):
    """Regression models should not have a classes_ attribute."""
    r = FOREST_REGRESSORS[name](random_state=0)
    assert_false(hasattr(r, "classes_"))
    assert_false(hasattr(r, "n_classes_"))

    r.fit([[1, 2, 3], [4, 5, 6]], [1, 2])
    assert_false(hasattr(r, "classes_"))
    assert_false(hasattr(r, "n_classes_"))


def test_regressor_attributes():
    for name in FOREST_REGRESSORS:
        yield check_regressor_attributes, name


def check_probability(name):
    """Predict probabilities."""
    ForestClassifier = FOREST_CLASSIFIERS[name]
    with np.errstate(divide="ignore"):
        clf = ForestClassifier(n_estimators=10, random_state=1, max_features=1,
                               max_depth=1)
        clf.fit(iris.data, iris.target)
        assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1),
                                  np.ones(iris.data.shape[0]))
        assert_array_almost_equal(clf.predict_proba(iris.data),
                                  np.exp(clf.predict_log_proba(iris.data)))


def test_probability():
    for name in FOREST_CLASSIFIERS:
        yield check_probability, name


def check_importance(name, X, y):
    """Check variable importances."""

    ForestClassifier = FOREST_CLASSIFIERS[name]

    clf = ForestClassifier(n_estimators=10)
    clf.fit(X, y)
    importances = clf.feature_importances_
    n_important = np.sum(importances > 0.1)
    assert_equal(importances.shape[0], 10)
    assert_equal(n_important, 3)

    X_new = clf.transform(X, threshold="mean")
    assert_less(0 < X_new.shape[1], X.shape[1])

    # Check with sample weights
    sample_weight = np.ones(y.shape)
    sample_weight[y == 1] *= 100

    clf = ForestClassifier(n_estimators=50, random_state=0)
    clf.fit(X, y, sample_weight=sample_weight)
    importances = clf.feature_importances_
    assert_true(np.all(importances >= 0.0))

    clf = ForestClassifier(n_estimators=50, random_state=0)
    clf.fit(X, y, sample_weight=3 * sample_weight)
    importances_bis = clf.feature_importances_
    assert_almost_equal(importances, importances_bis)


def test_importances():
    X, y = datasets.make_classification(n_samples=1000, n_features=10,
                                        n_informative=3, n_redundant=0,
                                        n_repeated=0, shuffle=False,
                                        random_state=0)

    for name in FOREST_CLASSIFIERS:
        yield check_importance, name, X, y

def check_oob_score(name, X, y, n_estimators=20):
    """Check that oob prediction is a good estimation of the generalization
       error."""
    est = FOREST_ESTIMATORS[name](oob_score=True, random_state=0,
                                  n_estimators=n_estimators)
    n_samples = X.shape[0]
    est.fit(X[:n_samples // 2, :], y[:n_samples // 2])
    test_score = est.score(X[n_samples // 2:, :], y[n_samples // 2:])

    if name in FOREST_CLASSIFIERS:
        assert_less(abs(test_score - est.oob_score_), 0.1)
    else:
        assert_greater(test_score, est.oob_score_)
        assert_greater(est.oob_score_, .8)

def test_oob_score():
    yield check_oob_score, "RandomForestClassifier", iris.data, iris.target
    yield (check_oob_score, "RandomForestRegressor", boston.data,
           boston.target, 50)

    # non-contiguous targets in classification
    yield (check_oob_score, "RandomForestClassifier", iris.data,
           iris.target * 2 + 1)

def check_gridsearch(name):
    forest = FOREST_CLASSIFIERS[name]()
    clf = GridSearchCV(forest, {'n_estimators': (1, 2), 'max_depth': (1, 2)})
    clf.fit(iris.data, iris.target)


def test_gridsearch():
    """Check that base trees can be grid-searched."""
    for name in FOREST_CLASSIFIERS:
        yield check_gridsearch, name


def check_parallel(name, X, y):
    """Check parallel computations in classification"""
    ForestEstimator = FOREST_ESTIMATORS[name]
    forest = ForestEstimator(n_estimators=10, n_jobs=3, random_state=0)

    forest.fit(X, y)
    assert_equal(len(forest), 10)

    forest.set_params(n_jobs=1)
    y1 = forest.predict(X)
    forest.set_params(n_jobs=2)
    y2 = forest.predict(X)
    assert_array_almost_equal(y1, y2, 3)


def test_parallel():
    for name in FOREST_CLASSIFIERS:
        yield check_parallel, name, iris.data, iris.target

    for name in FOREST_REGRESSORS:
        yield check_parallel, name,  boston.data, boston.target


def check_pickle(name, X, y):
    """Check pickability."""

    ForestEstimator = FOREST_ESTIMATORS[name]
    obj = ForestEstimator(random_state=0)
    obj.fit(X, y)
    score = obj.score(X, y)
    pickle_object = pickle.dumps(obj)

    obj2 = pickle.loads(pickle_object)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(X, y)
    assert_equal(score, score2)


def test_pickle():
    for name in FOREST_CLASSIFIERS:
        yield check_pickle, name, iris.data[::2], iris.target[::2]

    for name in FOREST_REGRESSORS:
        yield check_pickle, name,  boston.data[::2], boston.target[::2]


def check_multioutput(name,  X_train, X_test, y_train, y_test):
    """Check estimators on multi-output problems."""
    est = FOREST_ESTIMATORS[name](random_state=0, bootstrap=False)
    y_pred = est.fit(X_train, y_train).predict(X_test)
    assert_array_almost_equal(y_pred, y_test)

    if name in FOREST_CLASSIFIERS:
        with np.errstate(divide="ignore"):
            proba = est.predict_proba(X_test)
            assert_equal(len(proba), 2)
            assert_equal(proba[0].shape, (4, 2))
            assert_equal(proba[1].shape, (4, 4))

            log_proba = est.predict_log_proba(X_test)
            assert_equal(len(log_proba), 2)
            assert_equal(log_proba[0].shape, (4, 2))
            assert_equal(log_proba[1].shape, (4, 4))


def test_multioutput():
    X_train = [[-2, -1],
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

    y_train = [[-1, 0],
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

    X_test = [[-1, -1], [1, 1], [-1, 1], [1, -1]]
    y_test = [[-1, 0], [1, 1], [-1, 2], [1, 3]]

    for name in FOREST_CLASSIFIERS:
        yield check_multioutput, name, X_train, X_test, y_train, y_test

    for name in FOREST_REGRESSORS:
        yield check_multioutput, name, X_train, X_test, y_train, y_test


def check_classes_shape(name):
    """Test that n_classes_ and classes_ have proper shape."""
    ForestClassifier = FOREST_CLASSIFIERS[name]

    # Classification, single output
    clf = ForestClassifier(random_state=0).fit(X, y)

    assert_equal(clf.n_classes_, 2)
    assert_array_equal(clf.classes_, [-1, 1])

    # Classification, multi-output
    _y = np.vstack((y, np.array(y) * 2)).T
    clf = ForestClassifier(random_state=0).fit(X, _y)

    assert_array_equal(clf.n_classes_, [2, 2])
    assert_array_equal(clf.classes_, [[-1, 1], [-2, 2]])


def test_classes_shape():
    for name in FOREST_CLASSIFIERS:
        yield check_classes_shape, name


def test_random_trees_dense_type():
    '''
    Test that the `sparse_output` parameter of RandomTreesEmbedding
    works by returning a dense array.
    '''

    # Create the RTE with sparse=False
    hasher = RandomTreesEmbedding(n_estimators=10, sparse_output=False)
    X, y = datasets.make_circles(factor=0.5)
    X_transformed = hasher.fit_transform(X)

    # Assert that type is ndarray, not scipy.sparse.csr.csr_matrix
    assert_equal(type(X_transformed), np.ndarray)


def test_random_trees_dense_equal():
    '''
    Test that the `sparse_output` parameter of RandomTreesEmbedding
    works by returning the same array for both argument
    values.
    '''

    # Create the RTEs
    hasher_dense = RandomTreesEmbedding(n_estimators=10, sparse_output=False,
                                        random_state=0)
    hasher_sparse = RandomTreesEmbedding(n_estimators=10, sparse_output=True,
                                         random_state=0)
    X, y = datasets.make_circles(factor=0.5)
    X_transformed_dense = hasher_dense.fit_transform(X)
    X_transformed_sparse = hasher_sparse.fit_transform(X)

    # Assert that dense and sparse hashers have same array.
    assert_array_equal(X_transformed_sparse.toarray(), X_transformed_dense)


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
    rng = check_random_state(12321)
    n_samples, n_features = 80, 30
    X_train = rng.randn(n_samples, n_features)
    y_train = rng.randint(0, 2, n_samples)

    clfs = [
        RandomForestClassifier(n_estimators=20, n_jobs=n_jobs,
                               random_state=12345).fit(X_train, y_train)
        for n_jobs in [1, 2, 3, 8, 16, 32]
    ]

    X_test = rng.randn(n_samples, n_features)
    probas = [clf.predict_proba(X_test) for clf in clfs]
    for proba1, proba2 in zip(probas, probas[1:]):
        assert_array_almost_equal(proba1, proba2)


def test_distribution():
    rng = check_random_state(12321)

    # Single variable with 4 values
    X = rng.randint(0, 4, size=(1000, 1))
    y = rng.rand(1000)
    n_trees = 500

    clf = ExtraTreesRegressor(n_estimators=n_trees, random_state=42).fit(X, y)

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
    assert_greater(0.20, uniques[0][0])  # Rough approximation of 1/6.
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

    clf = ExtraTreesRegressor(n_estimators=100, max_features=1,
                              random_state=1).fit(X, y)

    uniques = defaultdict(int)
    for tree in clf.estimators_:
        tree = "".join(("%d,%d/" % (f, int(t)) if f >= 0 else "-")
                       for f, t in zip(tree.tree_.feature,
                                       tree.tree_.threshold))

        uniques[tree] += 1

    uniques = [(count, tree) for tree, count in uniques.items()]
    assert_equal(len(uniques), 8)


def check_max_leaf_nodes_max_depth(name, X, y):
    """Test precedence of max_leaf_nodes over max_depth. """
    ForestEstimator = FOREST_ESTIMATORS[name]
    est = ForestEstimator(max_depth=1, max_leaf_nodes=4,
                          n_estimators=1).fit(X, y)
    assert_greater(est.estimators_[0].tree_.max_depth, 1)

    est = ForestEstimator(max_depth=1, n_estimators=1).fit(X, y)
    assert_equal(est.estimators_[0].tree_.max_depth, 1)


def test_max_leaf_nodes_max_depth():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    for name in FOREST_ESTIMATORS:
        yield check_max_leaf_nodes_max_depth, name, X, y


if __name__ == "__main__":
    import nose
    nose.runmodule()
