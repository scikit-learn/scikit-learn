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
from distutils.version import LooseVersion
import itertools
from itertools import combinations
from itertools import product

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix

import pytest

from sklearn.utils._joblib import joblib
from sklearn.utils._joblib import parallel_backend
from sklearn.utils._joblib import register_parallel_backend
from sklearn.utils._joblib import __version__ as __joblib_version__

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_less, assert_greater
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import assert_no_warnings
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import skip_if_no_parallel

from sklearn import datasets
from sklearn.decomposition import TruncatedSVD
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import ExtraTreesRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomTreesEmbedding
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.utils.validation import check_random_state
from sklearn.utils.fixes import comb

from sklearn.tree.tree import SPARSE_SPLITTERS


# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# Larger classification sample used for testing feature importances
X_large, y_large = datasets.make_classification(
    n_samples=500, n_features=10, n_informative=3, n_redundant=0,
    n_repeated=0, shuffle=False, random_state=0)

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

# also make a hastie_10_2 dataset
hastie_X, hastie_y = datasets.make_hastie_10_2(n_samples=20, random_state=1)
hastie_X = hastie_X.astype(np.float32)

# Get the default backend in joblib to test parallelism and interaction with
# different backends
DEFAULT_JOBLIB_BACKEND = joblib.parallel.get_active_backend()[0].__class__

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

FOREST_CLASSIFIERS_REGRESSORS = FOREST_CLASSIFIERS.copy()
FOREST_CLASSIFIERS_REGRESSORS.update(FOREST_REGRESSORS)


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


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_classification_toy(name):
    check_classification_toy(name)


def check_iris_criterion(name, criterion):
    # Check consistency on dataset iris.
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


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
@pytest.mark.parametrize('criterion', ("gini", "entropy"))
def test_iris(name, criterion):
    check_iris_criterion(name, criterion)


def check_boston_criterion(name, criterion):
    # Check consistency on dataset boston house prices.
    ForestRegressor = FOREST_REGRESSORS[name]

    clf = ForestRegressor(n_estimators=5, criterion=criterion,
                          random_state=1)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert_greater(score, 0.94, "Failed with max_features=None, criterion %s "
                                "and score = %f" % (criterion, score))

    clf = ForestRegressor(n_estimators=5, criterion=criterion,
                          max_features=6, random_state=1)
    clf.fit(boston.data, boston.target)
    score = clf.score(boston.data, boston.target)
    assert_greater(score, 0.95, "Failed with max_features=6, criterion %s "
                                "and score = %f" % (criterion, score))


@pytest.mark.parametrize('name', FOREST_REGRESSORS)
@pytest.mark.parametrize('criterion', ("mse", "mae", "friedman_mse"))
def test_boston(name, criterion):
    check_boston_criterion(name, criterion)


def check_regressor_attributes(name):
    # Regression models should not have a classes_ attribute.
    r = FOREST_REGRESSORS[name](random_state=0)
    assert not hasattr(r, "classes_")
    assert not hasattr(r, "n_classes_")

    r.fit([[1, 2, 3], [4, 5, 6]], [1, 2])
    assert not hasattr(r, "classes_")
    assert not hasattr(r, "n_classes_")


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_REGRESSORS)
def test_regressor_attributes(name):
    check_regressor_attributes(name)


def check_probability(name):
    # Predict probabilities.
    ForestClassifier = FOREST_CLASSIFIERS[name]
    with np.errstate(divide="ignore"):
        clf = ForestClassifier(n_estimators=10, random_state=1, max_features=1,
                               max_depth=1)
        clf.fit(iris.data, iris.target)
        assert_array_almost_equal(np.sum(clf.predict_proba(iris.data), axis=1),
                                  np.ones(iris.data.shape[0]))
        assert_array_almost_equal(clf.predict_proba(iris.data),
                                  np.exp(clf.predict_log_proba(iris.data)))


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_probability(name):
    check_probability(name)


def check_importances(name, criterion, dtype, tolerance):
    # cast as dype
    X = X_large.astype(dtype, copy=False)
    y = y_large.astype(dtype, copy=False)

    ForestEstimator = FOREST_ESTIMATORS[name]

    est = ForestEstimator(n_estimators=10, criterion=criterion,
                          random_state=0)
    est.fit(X, y)
    importances = est.feature_importances_

    # The forest estimator can detect that only the first 3 features of the
    # dataset are informative:
    n_important = np.sum(importances > 0.1)
    assert_equal(importances.shape[0], 10)
    assert_equal(n_important, 3)
    assert np.all(importances[:3] > 0.1)

    # Check with parallel
    importances = est.feature_importances_
    est.set_params(n_jobs=2)
    importances_parallel = est.feature_importances_
    assert_array_almost_equal(importances, importances_parallel)

    # Check with sample weights
    sample_weight = check_random_state(0).randint(1, 10, len(X))
    est = ForestEstimator(n_estimators=10, random_state=0, criterion=criterion)
    est.fit(X, y, sample_weight=sample_weight)
    importances = est.feature_importances_
    assert np.all(importances >= 0.0)

    for scale in [0.5, 100]:
        est = ForestEstimator(n_estimators=10, random_state=0,
                              criterion=criterion)
        est.fit(X, y, sample_weight=scale * sample_weight)
        importances_bis = est.feature_importances_
        assert_less(np.abs(importances - importances_bis).mean(), tolerance)


@pytest.mark.parametrize('dtype', (np.float64, np.float32))
@pytest.mark.parametrize(
        'name, criterion',
        itertools.chain(product(FOREST_CLASSIFIERS,
                                ["gini", "entropy"]),
                        product(FOREST_REGRESSORS,
                                ["mse", "friedman_mse", "mae"])))
def test_importances(dtype, name, criterion):
    tolerance = 0.01
    if name in FOREST_REGRESSORS and criterion == "mae":
        tolerance = 0.05
    check_importances(name, criterion, dtype, tolerance)


def test_importances_asymptotic():
    # Check whether variable importances of totally randomized trees
    # converge towards their theoretical values (See Louppe et al,
    # Understanding variable importances in forests of randomized trees, 2013).

    def binomial(k, n):
        return 0 if k < 0 or k > n else comb(int(n), int(k), exact=True)

    def entropy(samples):
        n_samples = len(samples)
        entropy = 0.

        for count in np.bincount(samples):
            p = 1. * count / n_samples
            if p > 0:
                entropy -= p * np.log2(p)

        return entropy

    def mdi_importance(X_m, X, y):
        n_samples, n_features = X.shape

        features = list(range(n_features))
        features.pop(X_m)
        values = [np.unique(X[:, i]) for i in range(n_features)]

        imp = 0.

        for k in range(n_features):
            # Weight of each B of size k
            coef = 1. / (binomial(k, n_features) * (n_features - k))

            # For all B of size k
            for B in combinations(features, k):
                # For all values B=b
                for b in product(*[values[B[j]] for j in range(k)]):
                    mask_b = np.ones(n_samples, dtype=np.bool)

                    for j in range(k):
                        mask_b &= X[:, B[j]] == b[j]

                    X_, y_ = X[mask_b, :], y[mask_b]
                    n_samples_b = len(X_)

                    if n_samples_b > 0:
                        children = []

                        for xi in values[X_m]:
                            mask_xi = X_[:, X_m] == xi
                            children.append(y_[mask_xi])

                        imp += (coef
                                * (1. * n_samples_b / n_samples)  # P(B=b)
                                * (entropy(y_) -
                                   sum([entropy(c) * len(c) / n_samples_b
                                        for c in children])))

        return imp

    data = np.array([[0, 0, 1, 0, 0, 1, 0, 1],
                     [1, 0, 1, 1, 1, 0, 1, 2],
                     [1, 0, 1, 1, 0, 1, 1, 3],
                     [0, 1, 1, 1, 0, 1, 0, 4],
                     [1, 1, 0, 1, 0, 1, 1, 5],
                     [1, 1, 0, 1, 1, 1, 1, 6],
                     [1, 0, 1, 0, 0, 1, 0, 7],
                     [1, 1, 1, 1, 1, 1, 1, 8],
                     [1, 1, 1, 1, 0, 1, 1, 9],
                     [1, 1, 1, 0, 1, 1, 1, 0]])

    X, y = np.array(data[:, :7], dtype=np.bool), data[:, 7]
    n_features = X.shape[1]

    # Compute true importances
    true_importances = np.zeros(n_features)

    for i in range(n_features):
        true_importances[i] = mdi_importance(i, X, y)

    # Estimate importances with totally randomized trees
    clf = ExtraTreesClassifier(n_estimators=500,
                               max_features=1,
                               criterion="entropy",
                               random_state=0).fit(X, y)

    importances = sum(tree.tree_.compute_feature_importances(normalize=False)
                      for tree in clf.estimators_) / clf.n_estimators

    # Check correctness
    assert_almost_equal(entropy(y), sum(importances))
    assert_less(np.abs(true_importances - importances).mean(), 0.01)


def check_unfitted_feature_importances(name):
    assert_raises(ValueError, getattr, FOREST_ESTIMATORS[name](random_state=0),
                  "feature_importances_")


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_unfitted_feature_importances(name):
    check_unfitted_feature_importances(name)


def check_oob_score(name, X, y, n_estimators=20):
    # Check that oob prediction is a good estimation of the generalization
    # error.

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


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_oob_score_classifiers(name):
    check_oob_score(name, iris.data, iris.target)

    # csc matrix
    check_oob_score(name, csc_matrix(iris.data), iris.target)

    # non-contiguous targets in classification
    check_oob_score(name, iris.data, iris.target * 2 + 1)


@pytest.mark.parametrize('name', FOREST_REGRESSORS)
def test_oob_score_regressors(name):
    check_oob_score(name, boston.data, boston.target, 50)

    # csc matrix
    check_oob_score(name, csc_matrix(boston.data), boston.target, 50)


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
            assert not hasattr(est, "oob_score_")

        # No bootstrap
        assert_raises(ValueError, ForestEstimator(oob_score=True,
                                                  bootstrap=False).fit, X, y)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_oob_score_raise_error(name):
    check_oob_score_raise_error(name)


def check_gridsearch(name):
    forest = FOREST_CLASSIFIERS[name]()
    clf = GridSearchCV(forest, {'n_estimators': (1, 2), 'max_depth': (1, 2)})
    clf.fit(iris.data, iris.target)


@pytest.mark.filterwarnings('ignore: The default of the `iid`')  # 0.22
@pytest.mark.filterwarnings('ignore: The default value of cv')  # 0.22
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_gridsearch(name):
    # Check that base trees can be grid-searched.
    check_gridsearch(name)


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


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS_REGRESSORS)
def test_parallel(name):
    if name in FOREST_CLASSIFIERS:
        ds = iris
    elif name in FOREST_REGRESSORS:
        ds = boston

    check_parallel(name, ds.data, ds.target)


def check_pickle(name, X, y):
    # Check pickability.

    ForestEstimator = FOREST_ESTIMATORS[name]
    obj = ForestEstimator(random_state=0)
    obj.fit(X, y)
    score = obj.score(X, y)
    pickle_object = pickle.dumps(obj)

    obj2 = pickle.loads(pickle_object)
    assert_equal(type(obj2), obj.__class__)
    score2 = obj2.score(X, y)
    assert_equal(score, score2)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS_REGRESSORS)
def test_pickle(name):
    if name in FOREST_CLASSIFIERS:
        ds = iris
    elif name in FOREST_REGRESSORS:
        ds = boston

    check_pickle(name, ds.data[::2], ds.target[::2])


def check_multioutput(name):
    # Check estimators on multi-output problems.

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
            assert len(proba) == 2
            assert proba[0].shape == (4, 2)
            assert proba[1].shape == (4, 4)

            log_proba = est.predict_log_proba(X_test)
            assert len(log_proba) == 2
            assert log_proba[0].shape == (4, 2)
            assert log_proba[1].shape == (4, 4)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS_REGRESSORS)
def test_multioutput(name):
    check_multioutput(name)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_multioutput_string(name):
    # Check estimators on multi-output problems with string outputs.

    X_train = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1], [-2, 1],
               [-1, 1], [-1, 2], [2, -1], [1, -1], [1, -2]]
    y_train = [["red", "blue"], ["red", "blue"], ["red", "blue"],
               ["green", "green"], ["green", "green"], ["green", "green"],
               ["red", "purple"], ["red", "purple"], ["red", "purple"],
               ["green", "yellow"], ["green", "yellow"], ["green", "yellow"]]
    X_test = [[-1, -1], [1, 1], [-1, 1], [1, -1]]
    y_test = [["red", "blue"], ["green", "green"],
              ["red", "purple"], ["green", "yellow"]]

    est = FOREST_ESTIMATORS[name](random_state=0, bootstrap=False)
    y_pred = est.fit(X_train, y_train).predict(X_test)
    assert_array_equal(y_pred, y_test)

    with np.errstate(divide="ignore"):
        proba = est.predict_proba(X_test)
        assert len(proba) == 2
        assert proba[0].shape == (4, 2)
        assert proba[1].shape == (4, 4)

        log_proba = est.predict_log_proba(X_test)
        assert len(log_proba) == 2
        assert log_proba[0].shape == (4, 2)
        assert log_proba[1].shape == (4, 4)


def check_classes_shape(name):
    # Test that n_classes_ and classes_ have proper shape.
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


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_classes_shape(name):
    check_classes_shape(name)


def test_random_trees_dense_type():
    # Test that the `sparse_output` parameter of RandomTreesEmbedding
    # works by returning a dense array.

    # Create the RTE with sparse=False
    hasher = RandomTreesEmbedding(n_estimators=10, sparse_output=False)
    X, y = datasets.make_circles(factor=0.5)
    X_transformed = hasher.fit_transform(X)

    # Assert that type is ndarray, not scipy.sparse.csr.csr_matrix
    assert_equal(type(X_transformed), np.ndarray)


def test_random_trees_dense_equal():
    # Test that the `sparse_output` parameter of RandomTreesEmbedding
    # works by returning the same array for both argument values.

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


# Ignore warnings from switching to more power iterations in randomized_svd
@ignore_warnings
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


def test_random_hasher_sparse_data():
    X, y = datasets.make_multilabel_classification(random_state=0)
    hasher = RandomTreesEmbedding(n_estimators=30, random_state=1)
    X_transformed = hasher.fit_transform(X)
    X_transformed_sparse = hasher.fit_transform(csc_matrix(X))
    assert_array_equal(X_transformed_sparse.toarray(), X_transformed.toarray())


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


def check_max_leaf_nodes_max_depth(name):
    X, y = hastie_X, hastie_y

    # Test precedence of max_leaf_nodes over max_depth.
    ForestEstimator = FOREST_ESTIMATORS[name]
    est = ForestEstimator(max_depth=1, max_leaf_nodes=4,
                          n_estimators=1, random_state=0).fit(X, y)
    assert_equal(est.estimators_[0].get_depth(), 1)

    est = ForestEstimator(max_depth=1, n_estimators=1,
                          random_state=0).fit(X, y)
    assert_equal(est.estimators_[0].get_depth(), 1)


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_max_leaf_nodes_max_depth(name):
    check_max_leaf_nodes_max_depth(name)


def check_min_samples_split(name):
    X, y = hastie_X, hastie_y
    ForestEstimator = FOREST_ESTIMATORS[name]

    # test boundary value
    assert_raises(ValueError,
                  ForestEstimator(min_samples_split=-1).fit, X, y)
    assert_raises(ValueError,
                  ForestEstimator(min_samples_split=0).fit, X, y)
    assert_raises(ValueError,
                  ForestEstimator(min_samples_split=1.1).fit, X, y)

    est = ForestEstimator(min_samples_split=10, n_estimators=1, random_state=0)
    est.fit(X, y)
    node_idx = est.estimators_[0].tree_.children_left != -1
    node_samples = est.estimators_[0].tree_.n_node_samples[node_idx]

    assert_greater(np.min(node_samples), len(X) * 0.5 - 1,
                   "Failed with {0}".format(name))

    est = ForestEstimator(min_samples_split=0.5, n_estimators=1,
                          random_state=0)
    est.fit(X, y)
    node_idx = est.estimators_[0].tree_.children_left != -1
    node_samples = est.estimators_[0].tree_.n_node_samples[node_idx]

    assert_greater(np.min(node_samples), len(X) * 0.5 - 1,
                   "Failed with {0}".format(name))


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_min_samples_split(name):
    check_min_samples_split(name)


def check_min_samples_leaf(name):
    X, y = hastie_X, hastie_y

    # Test if leaves contain more than leaf_count training examples
    ForestEstimator = FOREST_ESTIMATORS[name]

    # test boundary value
    assert_raises(ValueError,
                  ForestEstimator(min_samples_leaf=-1).fit, X, y)
    assert_raises(ValueError,
                  ForestEstimator(min_samples_leaf=0).fit, X, y)

    est = ForestEstimator(min_samples_leaf=5, n_estimators=1, random_state=0)
    est.fit(X, y)
    out = est.estimators_[0].tree_.apply(X)
    node_counts = np.bincount(out)
    # drop inner nodes
    leaf_count = node_counts[node_counts != 0]
    assert_greater(np.min(leaf_count), 4,
                   "Failed with {0}".format(name))

    est = ForestEstimator(min_samples_leaf=0.25, n_estimators=1,
                          random_state=0)
    est.fit(X, y)
    out = est.estimators_[0].tree_.apply(X)
    node_counts = np.bincount(out)
    # drop inner nodes
    leaf_count = node_counts[node_counts != 0]
    assert_greater(np.min(leaf_count), len(X) * 0.25 - 1,
                   "Failed with {0}".format(name))


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_min_samples_leaf(name):
    check_min_samples_leaf(name)


def check_min_weight_fraction_leaf(name):
    X, y = hastie_X, hastie_y

    # Test if leaves contain at least min_weight_fraction_leaf of the
    # training set
    ForestEstimator = FOREST_ESTIMATORS[name]
    rng = np.random.RandomState(0)
    weights = rng.rand(X.shape[0])
    total_weight = np.sum(weights)

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for frac in np.linspace(0, 0.5, 6):
        est = ForestEstimator(min_weight_fraction_leaf=frac, n_estimators=1,
                              random_state=0)
        if "RandomForest" in name:
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


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_min_weight_fraction_leaf(name):
    check_min_weight_fraction_leaf(name)


def check_sparse_input(name, X, X_sparse, y):
    ForestEstimator = FOREST_ESTIMATORS[name]

    dense = ForestEstimator(random_state=0, max_depth=2).fit(X, y)
    sparse = ForestEstimator(random_state=0, max_depth=2).fit(X_sparse, y)

    assert_array_almost_equal(sparse.apply(X), dense.apply(X))

    if name in FOREST_CLASSIFIERS or name in FOREST_REGRESSORS:
        assert_array_almost_equal(sparse.predict(X), dense.predict(X))
        assert_array_almost_equal(sparse.feature_importances_,
                                  dense.feature_importances_)

    if name in FOREST_CLASSIFIERS:
        assert_array_almost_equal(sparse.predict_proba(X),
                                  dense.predict_proba(X))
        assert_array_almost_equal(sparse.predict_log_proba(X),
                                  dense.predict_log_proba(X))

    if name in FOREST_TRANSFORMERS:
        assert_array_almost_equal(sparse.transform(X).toarray(),
                                  dense.transform(X).toarray())
        assert_array_almost_equal(sparse.fit_transform(X).toarray(),
                                  dense.fit_transform(X).toarray())


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
@pytest.mark.parametrize('sparse_matrix',
                         (csr_matrix, csc_matrix, coo_matrix))
def test_sparse_input(name, sparse_matrix):
    X, y = datasets.make_multilabel_classification(random_state=0,
                                                   n_samples=50)

    check_sparse_input(name, X, sparse_matrix(X), y)


def check_memory_layout(name, dtype):
    # Check that it works no matter the memory layout

    est = FOREST_ESTIMATORS[name](random_state=0, bootstrap=False)

    # Nothing
    X = np.asarray(iris.data, dtype=dtype)
    y = iris.target
    assert_array_almost_equal(est.fit(X, y).predict(X), y)

    # C-order
    X = np.asarray(iris.data, order="C", dtype=dtype)
    y = iris.target
    assert_array_almost_equal(est.fit(X, y).predict(X), y)

    # F-order
    X = np.asarray(iris.data, order="F", dtype=dtype)
    y = iris.target
    assert_array_almost_equal(est.fit(X, y).predict(X), y)

    # Contiguous
    X = np.ascontiguousarray(iris.data, dtype=dtype)
    y = iris.target
    assert_array_almost_equal(est.fit(X, y).predict(X), y)

    if est.base_estimator.splitter in SPARSE_SPLITTERS:
        # csr matrix
        X = csr_matrix(iris.data, dtype=dtype)
        y = iris.target
        assert_array_almost_equal(est.fit(X, y).predict(X), y)

        # csc_matrix
        X = csc_matrix(iris.data, dtype=dtype)
        y = iris.target
        assert_array_almost_equal(est.fit(X, y).predict(X), y)

        # coo_matrix
        X = coo_matrix(iris.data, dtype=dtype)
        y = iris.target
        assert_array_almost_equal(est.fit(X, y).predict(X), y)

    # Strided
    X = np.asarray(iris.data[::3], dtype=dtype)
    y = iris.target[::3]
    assert_array_almost_equal(est.fit(X, y).predict(X), y)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS_REGRESSORS)
@pytest.mark.parametrize('dtype', (np.float64, np.float32))
def test_memory_layout(name, dtype):
    check_memory_layout(name, dtype)


@ignore_warnings
def check_1d_input(name, X, X_2d, y):
    ForestEstimator = FOREST_ESTIMATORS[name]
    assert_raises(ValueError, ForestEstimator(n_estimators=1,
                                              random_state=0).fit, X, y)

    est = ForestEstimator(random_state=0)
    est.fit(X_2d, y)

    if name in FOREST_CLASSIFIERS or name in FOREST_REGRESSORS:
        assert_raises(ValueError, est.predict, X)


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_1d_input(name):
    X = iris.data[:, 0]
    X_2d = iris.data[:, 0].reshape((-1, 1))
    y = iris.target

    with ignore_warnings():
        check_1d_input(name, X, X_2d, y)


def check_class_weights(name):
    # Check class_weights resemble sample_weights behavior.
    ForestClassifier = FOREST_CLASSIFIERS[name]

    # Iris is balanced, so no effect expected for using 'balanced' weights
    clf1 = ForestClassifier(random_state=0)
    clf1.fit(iris.data, iris.target)
    clf2 = ForestClassifier(class_weight='balanced', random_state=0)
    clf2.fit(iris.data, iris.target)
    assert_almost_equal(clf1.feature_importances_, clf2.feature_importances_)

    # Make a multi-output problem with three copies of Iris
    iris_multi = np.vstack((iris.target, iris.target, iris.target)).T
    # Create user-defined weights that should balance over the outputs
    clf3 = ForestClassifier(class_weight=[{0: 2., 1: 2., 2: 1.},
                                          {0: 2., 1: 1., 2: 2.},
                                          {0: 1., 1: 2., 2: 2.}],
                            random_state=0)
    clf3.fit(iris.data, iris_multi)
    assert_almost_equal(clf2.feature_importances_, clf3.feature_importances_)
    # Check against multi-output "balanced" which should also have no effect
    clf4 = ForestClassifier(class_weight='balanced', random_state=0)
    clf4.fit(iris.data, iris_multi)
    assert_almost_equal(clf3.feature_importances_, clf4.feature_importances_)

    # Inflate importance of class 1, check against user-defined weights
    sample_weight = np.ones(iris.target.shape)
    sample_weight[iris.target == 1] *= 100
    class_weight = {0: 1., 1: 100., 2: 1.}
    clf1 = ForestClassifier(random_state=0)
    clf1.fit(iris.data, iris.target, sample_weight)
    clf2 = ForestClassifier(class_weight=class_weight, random_state=0)
    clf2.fit(iris.data, iris.target)
    assert_almost_equal(clf1.feature_importances_, clf2.feature_importances_)

    # Check that sample_weight and class_weight are multiplicative
    clf1 = ForestClassifier(random_state=0)
    clf1.fit(iris.data, iris.target, sample_weight ** 2)
    clf2 = ForestClassifier(class_weight=class_weight, random_state=0)
    clf2.fit(iris.data, iris.target, sample_weight)
    assert_almost_equal(clf1.feature_importances_, clf2.feature_importances_)

    # Using a Python 2.x list as the sample_weight parameter used to raise
    # an exception. This test makes sure such code will now run correctly.
    clf = ForestClassifier()
    sample_weight = [1.] * len(iris.data)
    clf.fit(iris.data, iris.target, sample_weight=sample_weight)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_class_weights(name):
    check_class_weights(name)


def check_class_weight_balanced_and_bootstrap_multi_output(name):
    # Test class_weight works for multi-output"""
    ForestClassifier = FOREST_CLASSIFIERS[name]
    _y = np.vstack((y, np.array(y) * 2)).T
    clf = ForestClassifier(class_weight='balanced', random_state=0)
    clf.fit(X, _y)
    clf = ForestClassifier(class_weight=[{-1: 0.5, 1: 1.}, {-2: 1., 2: 1.}],
                           random_state=0)
    clf.fit(X, _y)
    # smoke test for balanced subsample
    clf = ForestClassifier(class_weight='balanced_subsample', random_state=0)
    clf.fit(X, _y)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_class_weight_balanced_and_bootstrap_multi_output(name):
    check_class_weight_balanced_and_bootstrap_multi_output(name)


def check_class_weight_errors(name):
    # Test if class_weight raises errors and warnings when expected.
    ForestClassifier = FOREST_CLASSIFIERS[name]
    _y = np.vstack((y, np.array(y) * 2)).T

    # Invalid preset string
    clf = ForestClassifier(class_weight='the larch', random_state=0)
    assert_raises(ValueError, clf.fit, X, y)
    assert_raises(ValueError, clf.fit, X, _y)

    # Warning warm_start with preset
    clf = ForestClassifier(class_weight='balanced', warm_start=True,
                           random_state=0)
    assert_warns(UserWarning, clf.fit, X, y)
    assert_warns(UserWarning, clf.fit, X, _y)

    # Not a list or preset for multi-output
    clf = ForestClassifier(class_weight=1, random_state=0)
    assert_raises(ValueError, clf.fit, X, _y)

    # Incorrect length list for multi-output
    clf = ForestClassifier(class_weight=[{-1: 0.5, 1: 1.}], random_state=0)
    assert_raises(ValueError, clf.fit, X, _y)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
def test_class_weight_errors(name):
    check_class_weight_errors(name)


def check_warm_start(name, random_state=42):
    # Test if fitting incrementally with warm start gives a forest of the
    # right size and the same results as a normal fit.
    X, y = hastie_X, hastie_y
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


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_warm_start(name):
    check_warm_start(name)


def check_warm_start_clear(name):
    # Test if fit clears state and grows a new forest when warm_start==False.
    X, y = hastie_X, hastie_y
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


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_warm_start_clear(name):
    check_warm_start_clear(name)


def check_warm_start_smaller_n_estimators(name):
    # Test if warm start second fit with smaller n_estimators raises error.
    X, y = hastie_X, hastie_y
    ForestEstimator = FOREST_ESTIMATORS[name]
    clf = ForestEstimator(n_estimators=5, max_depth=1, warm_start=True)
    clf.fit(X, y)
    clf.set_params(n_estimators=4)
    assert_raises(ValueError, clf.fit, X, y)


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_warm_start_smaller_n_estimators(name):
    check_warm_start_smaller_n_estimators(name)


def check_warm_start_equal_n_estimators(name):
    # Test if warm start with equal n_estimators does nothing and returns the
    # same forest and raises a warning.
    X, y = hastie_X, hastie_y
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


@pytest.mark.parametrize('name', FOREST_ESTIMATORS)
def test_warm_start_equal_n_estimators(name):
    check_warm_start_equal_n_estimators(name)


def check_warm_start_oob(name):
    # Test that the warm start computes oob score when asked.
    X, y = hastie_X, hastie_y
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

    assert hasattr(clf_2, 'oob_score_')
    assert_equal(clf.oob_score_, clf_2.oob_score_)

    # Test that oob_score is computed even if we don't need to train
    # additional trees.
    clf_3 = ForestEstimator(n_estimators=15, max_depth=3, warm_start=True,
                            random_state=1, bootstrap=True, oob_score=False)
    clf_3.fit(X, y)
    assert not hasattr(clf_3, 'oob_score_')

    clf_3.set_params(oob_score=True)
    ignore_warnings(clf_3.fit)(X, y)

    assert_equal(clf.oob_score_, clf_3.oob_score_)


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS_REGRESSORS)
def test_warm_start_oob(name):
    check_warm_start_oob(name)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
def test_dtype_convert(n_classes=15):
    classifier = RandomForestClassifier(random_state=0, bootstrap=False)

    X = np.eye(n_classes)
    y = [ch for ch in 'ABCDEFGHIJKLMNOPQRSTU'[:n_classes]]

    result = classifier.fit(X, y).predict(X)
    assert_array_equal(classifier.classes_, y)
    assert_array_equal(result, y)


def check_decision_path(name):
    X, y = hastie_X, hastie_y
    n_samples = X.shape[0]
    ForestEstimator = FOREST_ESTIMATORS[name]
    est = ForestEstimator(n_estimators=5, max_depth=1, warm_start=False,
                          random_state=1)
    est.fit(X, y)
    indicator, n_nodes_ptr = est.decision_path(X)

    assert_equal(indicator.shape[1], n_nodes_ptr[-1])
    assert_equal(indicator.shape[0], n_samples)
    assert_array_equal(np.diff(n_nodes_ptr),
                       [e.tree_.node_count for e in est.estimators_])

    # Assert that leaves index are correct
    leaves = est.apply(X)
    for est_id in range(leaves.shape[1]):
        leave_indicator = [indicator[i, n_nodes_ptr[est_id] + j]
                           for i, j in enumerate(leaves[:, est_id])]
        assert_array_almost_equal(leave_indicator, np.ones(shape=n_samples))


@pytest.mark.parametrize('name', FOREST_CLASSIFIERS_REGRESSORS)
def test_decision_path(name):
    check_decision_path(name)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
def test_min_impurity_split():
    # Test if min_impurity_split of base estimators is set
    # Regression test for #8006
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    all_estimators = [RandomForestClassifier, RandomForestRegressor,
                      ExtraTreesClassifier, ExtraTreesRegressor]

    for Estimator in all_estimators:
        est = Estimator(min_impurity_split=0.1)
        est = assert_warns_message(DeprecationWarning, "min_impurity_decrease",
                                   est.fit, X, y)
        for tree in est.estimators_:
            assert_equal(tree.min_impurity_split, 0.1)


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
def test_min_impurity_decrease():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    all_estimators = [RandomForestClassifier, RandomForestRegressor,
                      ExtraTreesClassifier, ExtraTreesRegressor]

    for Estimator in all_estimators:
        est = Estimator(min_impurity_decrease=0.1)
        est.fit(X, y)
        for tree in est.estimators_:
            # Simply check if the parameter is passed on correctly. Tree tests
            # will suffice for the actual working of this param
            assert_equal(tree.min_impurity_decrease, 0.1)


@pytest.mark.parametrize('forest',
                         [RandomForestClassifier, RandomForestRegressor,
                          ExtraTreesClassifier, ExtraTreesRegressor,
                          RandomTreesEmbedding])
def test_nestimators_future_warning(forest):
    # FIXME: to be removed 0.22

    # When n_estimators default value is used
    msg_future = ("The default value of n_estimators will change from "
                  "10 in version 0.20 to 100 in 0.22.")
    est = forest()
    est = assert_warns_message(FutureWarning, msg_future, est.fit, X, y)

    # When n_estimators is a valid value not equal to the default
    est = forest(n_estimators=100)
    est = assert_no_warnings(est.fit, X, y)


class MyBackend(DEFAULT_JOBLIB_BACKEND):
    def __init__(self, *args, **kwargs):
        self.count = 0
        super().__init__(*args, **kwargs)

    def start_call(self):
        self.count += 1
        return super().start_call()


register_parallel_backend('testing', MyBackend)


@pytest.mark.skipif(__joblib_version__ < LooseVersion('0.12'),
                    reason='tests not yet supported in joblib <0.12')
@skip_if_no_parallel
def test_backend_respected():
    clf = RandomForestClassifier(n_estimators=10, n_jobs=2)

    with parallel_backend("testing") as (ba, n_jobs):
        clf.fit(X, y)

    assert ba.count > 0

    # predict_proba requires shared memory. Ensure that's honored.
    with parallel_backend("testing") as (ba, _):
        clf.predict_proba(X)

    assert ba.count == 0


@pytest.mark.filterwarnings('ignore:The default value of n_estimators')
@pytest.mark.parametrize('name', FOREST_CLASSIFIERS)
@pytest.mark.parametrize('oob_score', (True, False))
def test_multi_target(name, oob_score):
    ForestClassifier = FOREST_CLASSIFIERS[name]

    clf = ForestClassifier(bootstrap=True, oob_score=oob_score)

    X = iris.data

    # Make multi column mixed type target.
    y = np.vstack([
        iris.target.astype(float),
        iris.target.astype(int),
        iris.target.astype(str),
    ]).T

    # Try to fit and predict.
    clf.fit(X, y)
    clf.predict(X)
