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
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import ignore_warnings

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
    for n_jobs in [1, 2]:
        clf = ForestClassifier(n_estimators=10, n_jobs=n_jobs)
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

        clf = ForestClassifier(n_estimators=50, n_jobs=n_jobs, random_state=0)
        clf.fit(X, y, sample_weight=sample_weight)
        importances = clf.feature_importances_
        assert_true(np.all(importances >= 0.0))

        clf = ForestClassifier(n_estimators=50, n_jobs=n_jobs, random_state=0)
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
    # Proper behavior
    est = FOREST_ESTIMATORS[name](oob_score=True, random_state=0,
                                  n_estimators=n_estimators, bootstrap=True)
    n_samples = X.shape[0]
    est.fit(X[:n_samples // 2, :], y[:n_samples // 2])
    test_score = est.score(X[n_samples // 2:, :], y[n_samples // 2:])

    if name in FOREST_CLASSIFIERS:
        assert_less(abs(test_score - est.oob_score_), 0.1)
    else:
        assert_greater(test_score, est.oob_score_)
        assert_greater(est.oob_score_, .8)

    # Check warning if not enough estimators
    with np.errstate(divide="ignore", invalid="ignore"):
        est = FOREST_ESTIMATORS[name](oob_score=True, random_state=0,
                                      n_estimators=1, bootstrap=True)
        assert_warns(UserWarning, est.fit, X, y)


def test_oob_score():
    for name in FOREST_CLASSIFIERS:
        yield check_oob_score, name, iris.data, iris.target

        # non-contiguous targets in classification
        yield check_oob_score, name, iris.data, iris.target * 2 + 1

    for name in FOREST_REGRESSORS:
        yield check_oob_score, name, boston.data, boston.target, 50


def check_oob_score_raise_error(name):
    ForestEstimator = FOREST_ESTIMATORS[name]

    if name in FOREST_TRANSFORMERS:
        for oob_score in [True, False]:
            assert_raises(TypeError, ForestEstimator, oob_score=oob_score)

        assert_raises(NotImplementedError, ForestEstimator()._set_oob_score,
                      X, y)

    else:
        # Unfitted /  no bootstrap / no oob_score
        for oob_score, bootstrap in [(True, False), (False, True),
                                     (False, False)]:
            est = ForestEstimator(oob_score=oob_score, bootstrap=bootstrap,
                                  random_state=0)
            assert_false(hasattr(est, "oob_score_"))

        # No bootstrap
        assert_raises(ValueError, ForestEstimator(oob_score=True,
                                                  bootstrap=False).fit, X, y)


def test_oob_score_raise_error():
    for name in FOREST_ESTIMATORS:
        yield check_oob_score_raise_error, name


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


def check_multioutput(name):
    """Check estimators on multi-output problems."""

    X_train = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1], [-2, 1],
               [-1, 1], [-1, 2], [2, -1], [1, -1], [1, -2]]

    y_train = [[-1, 0], [-1, 0], [-1, 0], [1, 1], [1, 1], [1, 1], [-1, 2],
               [-1, 2], [-1, 2], [1, 3], [1, 3], [1, 3]]
    X_test = [[-1, -1], [1, 1], [-1, 1], [1, -1]]
    y_test = [[-1, 0], [1, 1], [-1, 2], [1, 3]]

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
    for name in FOREST_CLASSIFIERS:
        yield check_multioutput, name

    for name in FOREST_REGRESSORS:
        yield check_multioutput, name


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


def check_min_samples_leaf(name, X, y):
    """Test if leaves contain more than leaf_count training examples"""
    ForestEstimator = FOREST_ESTIMATORS[name]

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes in (None, 1000):
        est = ForestEstimator(min_samples_leaf=5,
                              max_leaf_nodes=max_leaf_nodes,
                              random_state=0)
        est.fit(X, y)
        out = est.estimators_[0].tree_.apply(X)
        node_counts = np.bincount(out)
        # drop inner nodes
        leaf_count = node_counts[node_counts != 0]
        assert_greater(np.min(leaf_count), 4,
                       "Failed with {0}".format(name))


def test_min_samples_leaf():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    X = X.astype(np.float32)
    for name in FOREST_ESTIMATORS:
        yield check_min_samples_leaf, name, X, y


def check_min_weight_fraction_leaf(name, X, y):
    """Test if leaves contain at least min_weight_fraction_leaf of the
    training set"""
    ForestEstimator = FOREST_ESTIMATORS[name]
    rng = np.random.RandomState(0)
    weights = rng.rand(X.shape[0])
    total_weight = np.sum(weights)

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes in (None, 1000):
        for frac in np.linspace(0, 0.5, 6):
            est = ForestEstimator(min_weight_fraction_leaf=frac,
                                  max_leaf_nodes=max_leaf_nodes,
                                  random_state=0)
            if isinstance(est, (RandomForestClassifier,
                                RandomForestRegressor)):
                est.bootstrap = False
            est.fit(X, y, sample_weight=weights)
            out = est.estimators_[0].tree_.apply(X)
            node_weights = np.bincount(out, weights=weights)
            # drop inner nodes
            leaf_weights = node_weights[node_weights != 0]
            assert_greater_equal(
                np.min(leaf_weights),
                total_weight * est.min_weight_fraction_leaf,
                "Failed with {0} "
                "min_weight_fraction_leaf={1}".format(
                    name, est.min_weight_fraction_leaf))


def test_min_weight_fraction_leaf():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    X = X.astype(np.float32)
    for name in FOREST_ESTIMATORS:
        yield check_min_weight_fraction_leaf, name, X, y


def check_memory_layout(name, dtype):
    """Check that it works no matter the memory layout"""

    est = FOREST_ESTIMATORS[name](random_state=0, bootstrap=False)

    # Nothing
    X = np.asarray(iris.data, dtype=dtype)
    y = iris.target
    assert_array_equal(est.fit(X, y).predict(X), y)

    # C-order
    X = np.asarray(iris.data, order="C", dtype=dtype)
    y = iris.target
    assert_array_equal(est.fit(X, y).predict(X), y)

    # F-order
    X = np.asarray(iris.data, order="F", dtype=dtype)
    y = iris.target
    assert_array_equal(est.fit(X, y).predict(X), y)

    # Contiguous
    X = np.ascontiguousarray(iris.data, dtype=dtype)
    y = iris.target
    assert_array_equal(est.fit(X, y).predict(X), y)

    # Strided
    X = np.asarray(iris.data[::3], dtype=dtype)
    y = iris.target[::3]
    assert_array_equal(est.fit(X, y).predict(X), y)


def test_memory_layout():
    for name, dtype in product(FOREST_CLASSIFIERS, [np.float64, np.float32]):
        yield check_memory_layout, name, dtype

    for name, dtype in product(FOREST_REGRESSORS, [np.float64, np.float32]):
        yield check_memory_layout, name, dtype


def check_1d_input(name, X, X_2d, y):
    ForestEstimator = FOREST_ESTIMATORS[name]
    assert_raises(ValueError, ForestEstimator(random_state=0).fit, X, y)

    est = ForestEstimator(random_state=0)
    est.fit(X_2d, y)

    if name in FOREST_CLASSIFIERS or name in FOREST_REGRESSORS:
        assert_raises(ValueError, est.predict, X)


def test_1d_input():
    X = iris.data[:, 0].ravel()
    X_2d = iris.data[:, 0].reshape((-1, 1))
    y = iris.target

    for name in FOREST_ESTIMATORS:
        yield check_1d_input, name, X, X_2d, y


def check_warm_start(name, random_state=42):
    """Test if fitting incrementally with warm start gives a forest of the
    right size and the same results as a normal fit."""
    X, y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
    ForestEstimator = FOREST_ESTIMATORS[name]
    clf_ws = None
    for n_estimators in [5, 10]:
        if clf_ws is None:
            clf_ws = ForestEstimator(n_estimators=n_estimators,
                                     random_state=random_state,
                                     warm_start=True)
        else:
            clf_ws.set_params(n_estimators=n_estimators)
        clf_ws.fit(X, y)
        assert_equal(len(clf_ws), n_estimators)

    clf_no_ws = ForestEstimator(n_estimators=10, random_state=random_state,
                                warm_start=False)
    clf_no_ws.fit(X, y)

    assert_equal(set([tree.random_state for tree in clf_ws]),
                 set([tree.random_state for tree in clf_no_ws]))

    assert_array_equal(clf_ws.apply(X), clf_no_ws.apply(X),
                       err_msg="Failed with {0}".format(name))


def test_warm_start():
    for name in FOREST_ESTIMATORS:
        yield check_warm_start, name


def check_warm_start_clear(name):
    """Test if fit clears state and grows a new forest when warm_start==False.
    """
    X, y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
    ForestEstimator = FOREST_ESTIMATORS[name]
    clf = ForestEstimator(n_estimators=5, max_depth=1, warm_start=False,
                          random_state=1)
    clf.fit(X, y)

    clf_2 = ForestEstimator(n_estimators=5, max_depth=1, warm_start=True,
                            random_state=2)
    clf_2.fit(X, y)  # inits state
    clf_2.set_params(warm_start=False, random_state=1)
    clf_2.fit(X, y)  # clears old state and equals clf

    assert_array_almost_equal(clf_2.apply(X), clf.apply(X))


def test_warm_start_clear():
    for name in FOREST_ESTIMATORS:
        yield check_warm_start_clear, name


def check_warm_start_smaller_n_estimators(name):
    """Test if warm start second fit with smaller n_estimators raises error."""
    X, y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
    ForestEstimator = FOREST_ESTIMATORS[name]
    clf = ForestEstimator(n_estimators=5, max_depth=1, warm_start=True)
    clf.fit(X, y)
    clf.set_params(n_estimators=4)
    assert_raises(ValueError, clf.fit, X, y)


def test_warm_start_smaller_n_estimators():
    for name in FOREST_ESTIMATORS:
        yield check_warm_start_smaller_n_estimators, name


def check_warm_start_equal_n_estimators(name):
    """Test if warm start with equal n_estimators does nothing and returns the
    same forest and raises a warning."""
    X, y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
    ForestEstimator = FOREST_ESTIMATORS[name]
    clf = ForestEstimator(n_estimators=5, max_depth=3, warm_start=True,
                          random_state=1)
    clf.fit(X, y)

    clf_2 = ForestEstimator(n_estimators=5, max_depth=3, warm_start=True,
                            random_state=1)
    clf_2.fit(X, y)
    # Now clf_2 equals clf.

    clf_2.set_params(random_state=2)
    assert_warns(UserWarning, clf_2.fit, X, y)
    # If we had fit the trees again we would have got a different forest as we
    # changed the random state.
    assert_array_equal(clf.apply(X), clf_2.apply(X))


def test_warm_start_equal_n_estimators():
    for name in FOREST_ESTIMATORS:
        yield check_warm_start_equal_n_estimators, name


def check_warm_start_oob(name):
    """Test that the warm start computes oob score when asked."""
    X, y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
    ForestEstimator = FOREST_ESTIMATORS[name]
    # Use 15 estimators to avoid 'some inputs do not have OOB scores' warning.
    clf = ForestEstimator(n_estimators=15, max_depth=3, warm_start=False,
                          random_state=1, bootstrap=True, oob_score=True)
    clf.fit(X, y)

    clf_2 = ForestEstimator(n_estimators=5, max_depth=3, warm_start=False,
                            random_state=1, bootstrap=True, oob_score=False)
    clf_2.fit(X, y)

    clf_2.set_params(warm_start=True, oob_score=True, n_estimators=15)
    clf_2.fit(X, y)

    assert_true(hasattr(clf_2, 'oob_score_'))
    assert_equal(clf.oob_score_, clf_2.oob_score_)

    # Test that oob_score is computed even if we don't need to train
    # additional trees.
    clf_3 = ForestEstimator(n_estimators=15, max_depth=3, warm_start=True,
                            random_state=1, bootstrap=True, oob_score=False)
    clf_3.fit(X, y)
    assert_true(not(hasattr(clf_3, 'oob_score_')))

    clf_3.set_params(oob_score=True)
    ignore_warnings(clf_3.fit)(X, y)

    assert_equal(clf.oob_score_, clf_3.oob_score_)


def test_warm_start_oob():
    for name in FOREST_CLASSIFIERS:
        yield check_warm_start_oob, name
    for name in FOREST_REGRESSORS:
        yield check_warm_start_oob, name


if __name__ == "__main__":
    import nose
    nose.runmodule()
