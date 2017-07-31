"""
Testing for the tree module (sklearn.tree).
"""
import copy
import pickle
from functools import partial
from itertools import product
import struct

import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

from sklearn.random_projection import sparse_random_matrix

from sklearn.metrics import accuracy_score
from sklearn.metrics import mean_squared_error

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_greater_equal
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_less_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_warns
from sklearn.utils.testing import assert_warns_message
from sklearn.utils.testing import raises
from sklearn.utils.testing import ignore_warnings

from sklearn.utils.validation import check_random_state

from sklearn.exceptions import NotFittedError

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import ExtraTreeClassifier
from sklearn.tree import ExtraTreeRegressor

from sklearn import tree
from sklearn.tree._tree import TREE_LEAF
from sklearn.tree.tree import CRITERIA_CLF
from sklearn.tree.tree import CRITERIA_REG
from sklearn import datasets

from sklearn.utils import compute_sample_weight

CLF_CRITERIONS = ("gini", "entropy")
REG_CRITERIONS = ("mse", "mae", "friedman_mse")

CLF_TREES = {
    "DecisionTreeClassifier": DecisionTreeClassifier,
    "Presort-DecisionTreeClassifier": partial(DecisionTreeClassifier,
                                              presort=True),
    "ExtraTreeClassifier": ExtraTreeClassifier,
}

REG_TREES = {
    "DecisionTreeRegressor": DecisionTreeRegressor,
    "Presort-DecisionTreeRegressor": partial(DecisionTreeRegressor,
                                             presort=True),
    "ExtraTreeRegressor": ExtraTreeRegressor,
}

ALL_TREES = dict()
ALL_TREES.update(CLF_TREES)
ALL_TREES.update(REG_TREES)

SPARSE_TREES = ["DecisionTreeClassifier", "DecisionTreeRegressor",
                "ExtraTreeClassifier", "ExtraTreeRegressor"]


X_small = np.array([
    [0, 0, 4, 0, 0, 0, 1, -14, 0, -4, 0, 0, 0, 0, ],
    [0, 0, 5, 3, 0, -4, 0, 0, 1, -5, 0.2, 0, 4, 1, ],
    [-1, -1, 0, 0, -4.5, 0, 0, 2.1, 1, 0, 0, -4.5, 0, 1, ],
    [-1, -1, 0, -1.2, 0, 0, 0, 0, 0, 0, 0.2, 0, 0, 1, ],
    [-1, -1, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1, ],
    [-1, -2, 0, 4, -3, 10, 4, 0, -3.2, 0, 4, 3, -4, 1, ],
    [2.11, 0, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0.5, 0, -3, 1, ],
    [2.11, 0, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0, 0, -2, 1, ],
    [2.11, 8, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0, 0, -2, 1, ],
    [2.11, 8, -6, -0.5, 0, 11, 0, 0, -3.2, 6, 0.5, 0, -1, 0, ],
    [2, 8, 5, 1, 0.5, -4, 10, 0, 1, -5, 3, 0, 2, 0, ],
    [2, 0, 1, 1, 1, -1, 1, 0, 0, -2, 3, 0, 1, 0, ],
    [2, 0, 1, 2, 3, -1, 10, 2, 0, -1, 1, 2, 2, 0, ],
    [1, 1, 0, 2, 2, -1, 1, 2, 0, -5, 1, 2, 3, 0, ],
    [3, 1, 0, 3, 0, -4, 10, 0, 1, -5, 3, 0, 3, 1, ],
    [2.11, 8, -6, -0.5, 0, 1, 0, 0, -3.2, 6, 0.5, 0, -3, 1, ],
    [2.11, 8, -6, -0.5, 0, 1, 0, 0, -3.2, 6, 1.5, 1, -1, -1, ],
    [2.11, 8, -6, -0.5, 0, 10, 0, 0, -3.2, 6, 0.5, 0, -1, -1, ],
    [2, 0, 5, 1, 0.5, -2, 10, 0, 1, -5, 3, 1, 0, -1, ],
    [2, 0, 1, 1, 1, -2, 1, 0, 0, -2, 0, 0, 0, 1, ],
    [2, 1, 1, 1, 2, -1, 10, 2, 0, -1, 0, 2, 1, 1, ],
    [1, 1, 0, 0, 1, -3, 1, 2, 0, -5, 1, 2, 1, 1, ],
    [3, 1, 0, 1, 0, -4, 1, 0, 1, -2, 0, 0, 1, 0, ]])

y_small = [1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0,
           0, 0]
y_small_reg = [1.0, 2.1, 1.2, 0.05, 10, 2.4, 3.1, 1.01, 0.01, 2.98, 3.1, 1.1,
               0.0, 1.2, 2, 11, 0, 0, 4.5, 0.201, 1.06, 0.9, 0]

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
rng = np.random.RandomState(1)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

digits = datasets.load_digits()
perm = rng.permutation(digits.target.size)
digits.data = digits.data[perm]
digits.target = digits.target[perm]

random_state = check_random_state(0)
X_multilabel, y_multilabel = datasets.make_multilabel_classification(
    random_state=0, n_samples=30, n_features=10)

X_sparse_pos = random_state.uniform(size=(20, 5))
X_sparse_pos[X_sparse_pos <= 0.8] = 0.
y_random = random_state.randint(0, 4, size=(20, ))
X_sparse_mix = sparse_random_matrix(20, 10, density=0.25, random_state=0)


DATASETS = {
    "iris": {"X": iris.data, "y": iris.target},
    "boston": {"X": boston.data, "y": boston.target},
    "digits": {"X": digits.data, "y": digits.target},
    "toy": {"X": X, "y": y},
    "clf_small": {"X": X_small, "y": y_small},
    "reg_small": {"X": X_small, "y": y_small_reg},
    "multilabel": {"X": X_multilabel, "y": y_multilabel},
    "sparse-pos": {"X": X_sparse_pos, "y": y_random},
    "sparse-neg": {"X": - X_sparse_pos, "y": y_random},
    "sparse-mix": {"X": X_sparse_mix, "y": y_random},
    "zeros": {"X": np.zeros((20, 3)), "y": y_random}
}

for name in DATASETS:
    DATASETS[name]["X_sparse"] = csc_matrix(DATASETS[name]["X"])


def assert_tree_equal(d, s, message):
    assert_equal(s.node_count, d.node_count,
                 "{0}: inequal number of node ({1} != {2})"
                 "".format(message, s.node_count, d.node_count))

    assert_array_equal(d.children_right, s.children_right,
                       message + ": inequal children_right")
    assert_array_equal(d.children_left, s.children_left,
                       message + ": inequal children_left")

    external = d.children_right == TREE_LEAF
    internal = np.logical_not(external)

    assert_array_equal(d.feature[internal], s.feature[internal],
                       message + ": inequal features")
    assert_array_equal(d.threshold[internal], s.threshold[internal],
                       message + ": inequal threshold")
    assert_array_equal(d.n_node_samples.sum(), s.n_node_samples.sum(),
                       message + ": inequal sum(n_node_samples)")
    assert_array_equal(d.n_node_samples, s.n_node_samples,
                       message + ": inequal n_node_samples")

    assert_almost_equal(d.impurity, s.impurity,
                        err_msg=message + ": inequal impurity")

    assert_array_almost_equal(d.value[external], s.value[external],
                              err_msg=message + ": inequal value")


def test_classification_toy():
    # Check classification on a toy dataset.
    for name, Tree in CLF_TREES.items():
        clf = Tree(random_state=0)
        clf.fit(X, y)
        assert_array_equal(clf.predict(T), true_result,
                           "Failed with {0}".format(name))

        clf = Tree(max_features=1, random_state=1)
        clf.fit(X, y)
        assert_array_equal(clf.predict(T), true_result,
                           "Failed with {0}".format(name))


def test_weighted_classification_toy():
    # Check classification on a weighted toy dataset.
    for name, Tree in CLF_TREES.items():
        clf = Tree(random_state=0)

        clf.fit(X, y, sample_weight=np.ones(len(X)))
        assert_array_equal(clf.predict(T), true_result,
                           "Failed with {0}".format(name))

        clf.fit(X, y, sample_weight=np.ones(len(X)) * 0.5)
        assert_array_equal(clf.predict(T), true_result,
                           "Failed with {0}".format(name))


def test_regression_toy():
    # Check regression on a toy dataset.
    for name, Tree in REG_TREES.items():
        reg = Tree(random_state=1)
        reg.fit(X, y)
        assert_almost_equal(reg.predict(T), true_result,
                            err_msg="Failed with {0}".format(name))

        clf = Tree(max_features=1, random_state=1)
        clf.fit(X, y)
        assert_almost_equal(reg.predict(T), true_result,
                            err_msg="Failed with {0}".format(name))


def test_xor():
    # Check on a XOR problem
    y = np.zeros((10, 10))
    y[:5, :5] = 1
    y[5:, 5:] = 1

    gridx, gridy = np.indices(y.shape)

    X = np.vstack([gridx.ravel(), gridy.ravel()]).T
    y = y.ravel()

    for name, Tree in CLF_TREES.items():
        clf = Tree(random_state=0)
        clf.fit(X, y)
        assert_equal(clf.score(X, y), 1.0,
                     "Failed with {0}".format(name))

        clf = Tree(random_state=0, max_features=1)
        clf.fit(X, y)
        assert_equal(clf.score(X, y), 1.0,
                     "Failed with {0}".format(name))


def test_iris():
    # Check consistency on dataset iris.
    for (name, Tree), criterion in product(CLF_TREES.items(), CLF_CRITERIONS):
        clf = Tree(criterion=criterion, random_state=0)
        clf.fit(iris.data, iris.target)
        score = accuracy_score(clf.predict(iris.data), iris.target)
        assert_greater(score, 0.9,
                       "Failed with {0}, criterion = {1} and score = {2}"
                       "".format(name, criterion, score))

        clf = Tree(criterion=criterion, max_features=2, random_state=0)
        clf.fit(iris.data, iris.target)
        score = accuracy_score(clf.predict(iris.data), iris.target)
        assert_greater(score, 0.5,
                       "Failed with {0}, criterion = {1} and score = {2}"
                       "".format(name, criterion, score))


def test_boston():
    # Check consistency on dataset boston house prices.

    for (name, Tree), criterion in product(REG_TREES.items(), REG_CRITERIONS):
        reg = Tree(criterion=criterion, random_state=0)
        reg.fit(boston.data, boston.target)
        score = mean_squared_error(boston.target, reg.predict(boston.data))
        assert_less(score, 1,
                    "Failed with {0}, criterion = {1} and score = {2}"
                    "".format(name, criterion, score))

        # using fewer features reduces the learning ability of this tree,
        # but reduces training time.
        reg = Tree(criterion=criterion, max_features=6, random_state=0)
        reg.fit(boston.data, boston.target)
        score = mean_squared_error(boston.target, reg.predict(boston.data))
        assert_less(score, 2,
                    "Failed with {0}, criterion = {1} and score = {2}"
                    "".format(name, criterion, score))


def test_probability():
    # Predict probabilities using DecisionTreeClassifier.

    for name, Tree in CLF_TREES.items():
        clf = Tree(max_depth=1, max_features=1, random_state=42)
        clf.fit(iris.data, iris.target)

        prob_predict = clf.predict_proba(iris.data)
        assert_array_almost_equal(np.sum(prob_predict, 1),
                                  np.ones(iris.data.shape[0]),
                                  err_msg="Failed with {0}".format(name))
        assert_array_equal(np.argmax(prob_predict, 1),
                           clf.predict(iris.data),
                           err_msg="Failed with {0}".format(name))
        assert_almost_equal(clf.predict_proba(iris.data),
                            np.exp(clf.predict_log_proba(iris.data)), 8,
                            err_msg="Failed with {0}".format(name))


def test_arrayrepr():
    # Check the array representation.
    # Check resize
    X = np.arange(10000)[:, np.newaxis]
    y = np.arange(10000)

    for name, Tree in REG_TREES.items():
        reg = Tree(max_depth=None, random_state=0)
        reg.fit(X, y)


def test_pure_set():
    # Check when y is pure.
    X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
    y = [1, 1, 1, 1, 1, 1]

    for name, TreeClassifier in CLF_TREES.items():
        clf = TreeClassifier(random_state=0)
        clf.fit(X, y)
        assert_array_equal(clf.predict(X), y,
                           err_msg="Failed with {0}".format(name))

    for name, TreeRegressor in REG_TREES.items():
        reg = TreeRegressor(random_state=0)
        reg.fit(X, y)
        assert_almost_equal(reg.predict(X), y,
                            err_msg="Failed with {0}".format(name))


def test_numerical_stability():
    # Check numerical stability.
    X = np.array([
        [152.08097839, 140.40744019, 129.75102234, 159.90493774],
        [142.50700378, 135.81935120, 117.82884979, 162.75781250],
        [127.28772736, 140.40744019, 129.75102234, 159.90493774],
        [132.37025452, 143.71923828, 138.35694885, 157.84558105],
        [103.10237122, 143.71928406, 138.35696411, 157.84559631],
        [127.71276855, 143.71923828, 138.35694885, 157.84558105],
        [120.91514587, 140.40744019, 129.75102234, 159.90493774]])

    y = np.array(
        [1., 0.70209277, 0.53896582, 0., 0.90914464, 0.48026916, 0.49622521])

    with np.errstate(all="raise"):
        for name, Tree in REG_TREES.items():
            reg = Tree(random_state=0)
            reg.fit(X, y)
            reg.fit(X, -y)
            reg.fit(-X, y)
            reg.fit(-X, -y)


def test_importances():
    # Check variable importances.
    X, y = datasets.make_classification(n_samples=5000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=0)

    for name, Tree in CLF_TREES.items():
        clf = Tree(random_state=0)

        clf.fit(X, y)
        importances = clf.feature_importances_
        n_important = np.sum(importances > 0.1)

        assert_equal(importances.shape[0], 10, "Failed with {0}".format(name))
        assert_equal(n_important, 3, "Failed with {0}".format(name))

    # Check on iris that importances are the same for all builders
    clf = DecisionTreeClassifier(random_state=0)
    clf.fit(iris.data, iris.target)
    clf2 = DecisionTreeClassifier(random_state=0,
                                  max_leaf_nodes=len(iris.data))
    clf2.fit(iris.data, iris.target)

    assert_array_equal(clf.feature_importances_,
                       clf2.feature_importances_)


@raises(ValueError)
def test_importances_raises():
    # Check if variable importance before fit raises ValueError.
    clf = DecisionTreeClassifier()
    clf.feature_importances_


def test_importances_gini_equal_mse():
    # Check that gini is equivalent to mse for binary output variable

    X, y = datasets.make_classification(n_samples=2000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=0)

    # The gini index and the mean square error (variance) might differ due
    # to numerical instability. Since those instabilities mainly occurs at
    # high tree depth, we restrict this maximal depth.
    clf = DecisionTreeClassifier(criterion="gini", max_depth=5,
                                 random_state=0).fit(X, y)
    reg = DecisionTreeRegressor(criterion="mse", max_depth=5,
                                random_state=0).fit(X, y)

    assert_almost_equal(clf.feature_importances_, reg.feature_importances_)
    assert_array_equal(clf.tree_.feature, reg.tree_.feature)
    assert_array_equal(clf.tree_.children_left, reg.tree_.children_left)
    assert_array_equal(clf.tree_.children_right, reg.tree_.children_right)
    assert_array_equal(clf.tree_.n_node_samples, reg.tree_.n_node_samples)


def test_max_features():
    # Check max_features.
    for name, TreeRegressor in REG_TREES.items():
        reg = TreeRegressor(max_features="auto")
        reg.fit(boston.data, boston.target)
        assert_equal(reg.max_features_, boston.data.shape[1])

    for name, TreeClassifier in CLF_TREES.items():
        clf = TreeClassifier(max_features="auto")
        clf.fit(iris.data, iris.target)
        assert_equal(clf.max_features_, 2)

    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(max_features="sqrt")
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_,
                     int(np.sqrt(iris.data.shape[1])))

        est = TreeEstimator(max_features="log2")
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_,
                     int(np.log2(iris.data.shape[1])))

        est = TreeEstimator(max_features=1)
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_, 1)

        est = TreeEstimator(max_features=3)
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_, 3)

        est = TreeEstimator(max_features=0.01)
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_, 1)

        est = TreeEstimator(max_features=0.5)
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_,
                     int(0.5 * iris.data.shape[1]))

        est = TreeEstimator(max_features=1.0)
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_, iris.data.shape[1])

        est = TreeEstimator(max_features=None)
        est.fit(iris.data, iris.target)
        assert_equal(est.max_features_, iris.data.shape[1])

        # use values of max_features that are invalid
        est = TreeEstimator(max_features=10)
        assert_raises(ValueError, est.fit, X, y)

        est = TreeEstimator(max_features=-1)
        assert_raises(ValueError, est.fit, X, y)

        est = TreeEstimator(max_features=0.0)
        assert_raises(ValueError, est.fit, X, y)

        est = TreeEstimator(max_features=1.5)
        assert_raises(ValueError, est.fit, X, y)

        est = TreeEstimator(max_features="foobar")
        assert_raises(ValueError, est.fit, X, y)


def test_error():
    # Test that it gives proper exception on deficient input.
    for name, TreeEstimator in CLF_TREES.items():
        # predict before fit
        est = TreeEstimator()
        assert_raises(NotFittedError, est.predict_proba, X)

        est.fit(X, y)
        X2 = [[-2, -1, 1]]  # wrong feature shape for sample
        assert_raises(ValueError, est.predict_proba, X2)

    for name, TreeEstimator in ALL_TREES.items():
        # Invalid values for parameters
        assert_raises(ValueError, TreeEstimator(min_samples_leaf=-1).fit, X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_leaf=.6).fit, X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_leaf=0.).fit, X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_leaf=3.).fit, X, y)
        assert_raises(ValueError,
                      TreeEstimator(min_weight_fraction_leaf=-1).fit,
                      X, y)
        assert_raises(ValueError,
                      TreeEstimator(min_weight_fraction_leaf=0.51).fit,
                      X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_split=-1).fit,
                      X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_split=0.0).fit,
                      X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_split=1.1).fit,
                      X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_split=2.5).fit,
                      X, y)
        assert_raises(ValueError, TreeEstimator(max_depth=-1).fit, X, y)
        assert_raises(ValueError, TreeEstimator(max_features=42).fit, X, y)
        assert_raises(ValueError, TreeEstimator(min_impurity_split=-1.0).fit,
                      X, y)
        assert_raises(ValueError,
                      TreeEstimator(min_impurity_decrease=-1.0).fit, X, y)

        # Wrong dimensions
        est = TreeEstimator()
        y2 = y[:-1]
        assert_raises(ValueError, est.fit, X, y2)

        # Test with arrays that are non-contiguous.
        Xf = np.asfortranarray(X)
        est = TreeEstimator()
        est.fit(Xf, y)
        assert_almost_equal(est.predict(T), true_result)

        # predict before fitting
        est = TreeEstimator()
        assert_raises(NotFittedError, est.predict, T)

        # predict on vector with different dims
        est.fit(X, y)
        t = np.asarray(T)
        assert_raises(ValueError, est.predict, t[:, 1:])

        # wrong sample shape
        Xt = np.array(X).T

        est = TreeEstimator()
        est.fit(np.dot(X, Xt), y)
        assert_raises(ValueError, est.predict, X)
        assert_raises(ValueError, est.apply, X)

        clf = TreeEstimator()
        clf.fit(X, y)
        assert_raises(ValueError, clf.predict, Xt)
        assert_raises(ValueError, clf.apply, Xt)

        # apply before fitting
        est = TreeEstimator()
        assert_raises(NotFittedError, est.apply, T)


def test_min_samples_split():
    """Test min_samples_split parameter"""
    X = np.asfortranarray(iris.data.astype(tree._tree.DTYPE))
    y = iris.target

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes, name in product((None, 1000), ALL_TREES.keys()):
        TreeEstimator = ALL_TREES[name]

        # test for integer parameter
        est = TreeEstimator(min_samples_split=10,
                            max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        est.fit(X, y)
        # count samples on nodes, -1 means it is a leaf
        node_samples = est.tree_.n_node_samples[est.tree_.children_left != -1]

        assert_greater(np.min(node_samples), 9,
                       "Failed with {0}".format(name))

        # test for float parameter
        est = TreeEstimator(min_samples_split=0.2,
                            max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        est.fit(X, y)
        # count samples on nodes, -1 means it is a leaf
        node_samples = est.tree_.n_node_samples[est.tree_.children_left != -1]

        assert_greater(np.min(node_samples), 9,
                       "Failed with {0}".format(name))


def test_min_samples_leaf():
    # Test if leaves contain more than leaf_count training examples
    X = np.asfortranarray(iris.data.astype(tree._tree.DTYPE))
    y = iris.target

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes, name in product((None, 1000), ALL_TREES.keys()):
        TreeEstimator = ALL_TREES[name]

        # test integer parameter
        est = TreeEstimator(min_samples_leaf=5,
                            max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        est.fit(X, y)
        out = est.tree_.apply(X)
        node_counts = np.bincount(out)
        # drop inner nodes
        leaf_count = node_counts[node_counts != 0]
        assert_greater(np.min(leaf_count), 4,
                       "Failed with {0}".format(name))

        # test float parameter
        est = TreeEstimator(min_samples_leaf=0.1,
                            max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        est.fit(X, y)
        out = est.tree_.apply(X)
        node_counts = np.bincount(out)
        # drop inner nodes
        leaf_count = node_counts[node_counts != 0]
        assert_greater(np.min(leaf_count), 4,
                       "Failed with {0}".format(name))


def check_min_weight_fraction_leaf(name, datasets, sparse=False):
    """Test if leaves contain at least min_weight_fraction_leaf of the
    training set"""
    if sparse:
        X = DATASETS[datasets]["X_sparse"].astype(np.float32)
    else:
        X = DATASETS[datasets]["X"].astype(np.float32)
    y = DATASETS[datasets]["y"]

    weights = rng.rand(X.shape[0])
    total_weight = np.sum(weights)

    TreeEstimator = ALL_TREES[name]

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes, frac in product((None, 1000), np.linspace(0, 0.5, 6)):
        est = TreeEstimator(min_weight_fraction_leaf=frac,
                            max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        est.fit(X, y, sample_weight=weights)

        if sparse:
            out = est.tree_.apply(X.tocsr())

        else:
            out = est.tree_.apply(X)

        node_weights = np.bincount(out, weights=weights)
        # drop inner nodes
        leaf_weights = node_weights[node_weights != 0]
        assert_greater_equal(
            np.min(leaf_weights),
            total_weight * est.min_weight_fraction_leaf,
            "Failed with {0} "
            "min_weight_fraction_leaf={1}".format(
                name, est.min_weight_fraction_leaf))

    # test case with no weights passed in
    total_weight = X.shape[0]

    for max_leaf_nodes, frac in product((None, 1000), np.linspace(0, 0.5, 6)):
        est = TreeEstimator(min_weight_fraction_leaf=frac,
                            max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        est.fit(X, y)

        if sparse:
            out = est.tree_.apply(X.tocsr())
        else:
            out = est.tree_.apply(X)

        node_weights = np.bincount(out)
        # drop inner nodes
        leaf_weights = node_weights[node_weights != 0]
        assert_greater_equal(
            np.min(leaf_weights),
            total_weight * est.min_weight_fraction_leaf,
            "Failed with {0} "
            "min_weight_fraction_leaf={1}".format(
                name, est.min_weight_fraction_leaf))


def test_min_weight_fraction_leaf():
    # Check on dense input
    for name in ALL_TREES:
        yield check_min_weight_fraction_leaf, name, "iris"

    # Check on sparse input
    for name in SPARSE_TREES:
        yield check_min_weight_fraction_leaf, name, "multilabel", True


def check_min_weight_fraction_leaf_with_min_samples_leaf(name, datasets,
                                                         sparse=False):
    """Test the interaction between min_weight_fraction_leaf and min_samples_leaf
    when sample_weights is not provided in fit."""
    if sparse:
        X = DATASETS[datasets]["X_sparse"].astype(np.float32)
    else:
        X = DATASETS[datasets]["X"].astype(np.float32)
    y = DATASETS[datasets]["y"]

    total_weight = X.shape[0]
    TreeEstimator = ALL_TREES[name]
    for max_leaf_nodes, frac in product((None, 1000), np.linspace(0, 0.5, 3)):
        # test integer min_samples_leaf
        est = TreeEstimator(min_weight_fraction_leaf=frac,
                            max_leaf_nodes=max_leaf_nodes,
                            min_samples_leaf=5,
                            random_state=0)
        est.fit(X, y)

        if sparse:
            out = est.tree_.apply(X.tocsr())
        else:
            out = est.tree_.apply(X)

        node_weights = np.bincount(out)
        # drop inner nodes
        leaf_weights = node_weights[node_weights != 0]
        assert_greater_equal(
            np.min(leaf_weights),
            max((total_weight *
                 est.min_weight_fraction_leaf), 5),
            "Failed with {0} "
            "min_weight_fraction_leaf={1}, "
            "min_samples_leaf={2}".format(name,
                                          est.min_weight_fraction_leaf,
                                          est.min_samples_leaf))
    for max_leaf_nodes, frac in product((None, 1000), np.linspace(0, 0.5, 3)):
        # test float min_samples_leaf
        est = TreeEstimator(min_weight_fraction_leaf=frac,
                            max_leaf_nodes=max_leaf_nodes,
                            min_samples_leaf=.1,
                            random_state=0)
        est.fit(X, y)

        if sparse:
            out = est.tree_.apply(X.tocsr())
        else:
            out = est.tree_.apply(X)

        node_weights = np.bincount(out)
        # drop inner nodes
        leaf_weights = node_weights[node_weights != 0]
        assert_greater_equal(
            np.min(leaf_weights),
            max((total_weight * est.min_weight_fraction_leaf),
                (total_weight * est.min_samples_leaf)),
            "Failed with {0} "
            "min_weight_fraction_leaf={1}, "
            "min_samples_leaf={2}".format(name,
                                          est.min_weight_fraction_leaf,
                                          est.min_samples_leaf))


def test_min_weight_fraction_leaf_with_min_samples_leaf():
    # Check on dense input
    for name in ALL_TREES:
        yield (check_min_weight_fraction_leaf_with_min_samples_leaf,
               name, "iris")

    # Check on sparse input
    for name in SPARSE_TREES:
        yield (check_min_weight_fraction_leaf_with_min_samples_leaf,
               name, "multilabel", True)


def test_min_impurity_split():
    # test if min_impurity_split creates leaves with impurity
    # [0, min_impurity_split) when min_samples_leaf = 1 and
    # min_samples_split = 2.
    X = np.asfortranarray(iris.data.astype(tree._tree.DTYPE))
    y = iris.target

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes, name in product((None, 1000), ALL_TREES.keys()):
        TreeEstimator = ALL_TREES[name]
        min_impurity_split = .5

        # verify leaf nodes without min_impurity_split less than
        # impurity 1e-7
        est = TreeEstimator(max_leaf_nodes=max_leaf_nodes,
                            random_state=0)
        assert_true(est.min_impurity_split is None,
                    "Failed, min_impurity_split = {0} > 1e-7".format(
                        est.min_impurity_split))
        try:
            assert_warns(DeprecationWarning, est.fit, X, y)
        except AssertionError:
            pass
        for node in range(est.tree_.node_count):
            if (est.tree_.children_left[node] == TREE_LEAF or
                    est.tree_.children_right[node] == TREE_LEAF):
                assert_equal(est.tree_.impurity[node], 0.,
                             "Failed with {0} "
                             "min_impurity_split={1}".format(
                                 est.tree_.impurity[node],
                                 est.min_impurity_split))

        # verify leaf nodes have impurity [0,min_impurity_split] when using
        # min_impurity_split
        est = TreeEstimator(max_leaf_nodes=max_leaf_nodes,
                            min_impurity_split=min_impurity_split,
                            random_state=0)
        assert_warns_message(DeprecationWarning,
                             "Use the min_impurity_decrease",
                             est.fit, X, y)
        for node in range(est.tree_.node_count):
            if (est.tree_.children_left[node] == TREE_LEAF or
                    est.tree_.children_right[node] == TREE_LEAF):
                assert_greater_equal(est.tree_.impurity[node], 0,
                                     "Failed with {0}, "
                                     "min_impurity_split={1}".format(
                                         est.tree_.impurity[node],
                                         est.min_impurity_split))
                assert_less_equal(est.tree_.impurity[node], min_impurity_split,
                                  "Failed with {0}, "
                                  "min_impurity_split={1}".format(
                                      est.tree_.impurity[node],
                                      est.min_impurity_split))


def test_min_impurity_decrease():
    # test if min_impurity_decrease ensure that a split is made only if
    # if the impurity decrease is atleast that value
    X, y = datasets.make_classification(n_samples=10000, random_state=42)

    # test both DepthFirstTreeBuilder and BestFirstTreeBuilder
    # by setting max_leaf_nodes
    for max_leaf_nodes, name in product((None, 1000), ALL_TREES.keys()):
        TreeEstimator = ALL_TREES[name]

        # Check default value of min_impurity_decrease, 1e-7
        est1 = TreeEstimator(max_leaf_nodes=max_leaf_nodes, random_state=0)
        # Check with explicit value of 0.05
        est2 = TreeEstimator(max_leaf_nodes=max_leaf_nodes,
                             min_impurity_decrease=0.05, random_state=0)
        # Check with a much lower value of 0.0001
        est3 = TreeEstimator(max_leaf_nodes=max_leaf_nodes,
                             min_impurity_decrease=0.0001, random_state=0)
        # Check with a much lower value of 0.1
        est4 = TreeEstimator(max_leaf_nodes=max_leaf_nodes,
                             min_impurity_decrease=0.1, random_state=0)

        for est, expected_decrease in ((est1, 1e-7), (est2, 0.05),
                                       (est3, 0.0001), (est4, 0.1)):
            assert_less_equal(est.min_impurity_decrease, expected_decrease,
                              "Failed, min_impurity_decrease = {0} > {1}"
                              .format(est.min_impurity_decrease,
                                      expected_decrease))
            est.fit(X, y)
            for node in range(est.tree_.node_count):
                # If current node is a not leaf node, check if the split was
                # justified w.r.t the min_impurity_decrease
                if est.tree_.children_left[node] != TREE_LEAF:
                    imp_parent = est.tree_.impurity[node]
                    wtd_n_node = est.tree_.weighted_n_node_samples[node]

                    left = est.tree_.children_left[node]
                    wtd_n_left = est.tree_.weighted_n_node_samples[left]
                    imp_left = est.tree_.impurity[left]
                    wtd_imp_left = wtd_n_left * imp_left

                    right = est.tree_.children_right[node]
                    wtd_n_right = est.tree_.weighted_n_node_samples[right]
                    imp_right = est.tree_.impurity[right]
                    wtd_imp_right = wtd_n_right * imp_right

                    wtd_avg_left_right_imp = wtd_imp_right + wtd_imp_left
                    wtd_avg_left_right_imp /= wtd_n_node

                    fractional_node_weight = (
                        est.tree_.weighted_n_node_samples[node] / X.shape[0])

                    actual_decrease = fractional_node_weight * (
                        imp_parent - wtd_avg_left_right_imp)

                    assert_greater_equal(actual_decrease, expected_decrease,
                                         "Failed with {0} "
                                         "expected min_impurity_decrease={1}"
                                         .format(actual_decrease,
                                                 expected_decrease))

    for name, TreeEstimator in ALL_TREES.items():
        if "Classifier" in name:
            X, y = iris.data, iris.target
        else:
            X, y = boston.data, boston.target

        est = TreeEstimator(random_state=0)
        est.fit(X, y)
        score = est.score(X, y)
        fitted_attribute = dict()
        for attribute in ["max_depth", "node_count", "capacity"]:
            fitted_attribute[attribute] = getattr(est.tree_, attribute)

        serialized_object = pickle.dumps(est)
        est2 = pickle.loads(serialized_object)
        assert_equal(type(est2), est.__class__)
        score2 = est2.score(X, y)
        assert_equal(score, score2,
                     "Failed to generate same score  after pickling "
                     "with {0}".format(name))

        for attribute in fitted_attribute:
            assert_equal(getattr(est2.tree_, attribute),
                         fitted_attribute[attribute],
                         "Failed to generate same attribute {0} after "
                         "pickling with {1}".format(attribute, name))


def test_multioutput():
    # Check estimators on multi-output problems.
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
    for name, TreeClassifier in CLF_TREES.items():
        clf = TreeClassifier(random_state=0)
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
    for name, TreeRegressor in REG_TREES.items():
        reg = TreeRegressor(random_state=0)
        y_hat = reg.fit(X, y).predict(T)
        assert_almost_equal(y_hat, y_true)
        assert_equal(y_hat.shape, (4, 2))


def test_classes_shape():
    # Test that n_classes_ and classes_ have proper shape.
    for name, TreeClassifier in CLF_TREES.items():
        # Classification, single output
        clf = TreeClassifier(random_state=0)
        clf.fit(X, y)

        assert_equal(clf.n_classes_, 2)
        assert_array_equal(clf.classes_, [-1, 1])

        # Classification, multi-output
        _y = np.vstack((y, np.array(y) * 2)).T
        clf = TreeClassifier(random_state=0)
        clf.fit(X, _y)
        assert_equal(len(clf.n_classes_), 2)
        assert_equal(len(clf.classes_), 2)
        assert_array_equal(clf.n_classes_, [2, 2])
        assert_array_equal(clf.classes_, [[-1, 1], [-2, 2]])


def test_unbalanced_iris():
    # Check class rebalancing.
    unbalanced_X = iris.data[:125]
    unbalanced_y = iris.target[:125]
    sample_weight = compute_sample_weight("balanced", unbalanced_y)

    for name, TreeClassifier in CLF_TREES.items():
        clf = TreeClassifier(random_state=0)
        clf.fit(unbalanced_X, unbalanced_y, sample_weight=sample_weight)
        assert_almost_equal(clf.predict(unbalanced_X), unbalanced_y)


def test_memory_layout():
    # Check that it works no matter the memory layout
    for (name, TreeEstimator), dtype in product(ALL_TREES.items(),
                                                [np.float64, np.float32]):
        est = TreeEstimator(random_state=0)

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

        if not est.presort:
            # csr matrix
            X = csr_matrix(iris.data, dtype=dtype)
            y = iris.target
            assert_array_equal(est.fit(X, y).predict(X), y)

            # csc_matrix
            X = csc_matrix(iris.data, dtype=dtype)
            y = iris.target
            assert_array_equal(est.fit(X, y).predict(X), y)

        # Strided
        X = np.asarray(iris.data[::3], dtype=dtype)
        y = iris.target[::3]
        assert_array_equal(est.fit(X, y).predict(X), y)


def test_sample_weight():
    # Check sample weighting.
    # Test that zero-weighted samples are not taken into account
    X = np.arange(100)[:, np.newaxis]
    y = np.ones(100)
    y[:50] = 0.0

    sample_weight = np.ones(100)
    sample_weight[y == 0] = 0.0

    clf = DecisionTreeClassifier(random_state=0)
    clf.fit(X, y, sample_weight=sample_weight)
    assert_array_equal(clf.predict(X), np.ones(100))

    # Test that low weighted samples are not taken into account at low depth
    X = np.arange(200)[:, np.newaxis]
    y = np.zeros(200)
    y[50:100] = 1
    y[100:200] = 2
    X[100:200, 0] = 200

    sample_weight = np.ones(200)

    sample_weight[y == 2] = .51  # Samples of class '2' are still weightier
    clf = DecisionTreeClassifier(max_depth=1, random_state=0)
    clf.fit(X, y, sample_weight=sample_weight)
    assert_equal(clf.tree_.threshold[0], 149.5)

    sample_weight[y == 2] = .5  # Samples of class '2' are no longer weightier
    clf = DecisionTreeClassifier(max_depth=1, random_state=0)
    clf.fit(X, y, sample_weight=sample_weight)
    assert_equal(clf.tree_.threshold[0], 49.5)  # Threshold should have moved

    # Test that sample weighting is the same as having duplicates
    X = iris.data
    y = iris.target

    duplicates = rng.randint(0, X.shape[0], 100)

    clf = DecisionTreeClassifier(random_state=1)
    clf.fit(X[duplicates], y[duplicates])

    sample_weight = np.bincount(duplicates, minlength=X.shape[0])
    clf2 = DecisionTreeClassifier(random_state=1)
    clf2.fit(X, y, sample_weight=sample_weight)

    internal = clf.tree_.children_left != tree._tree.TREE_LEAF
    assert_array_almost_equal(clf.tree_.threshold[internal],
                              clf2.tree_.threshold[internal])


def test_sample_weight_invalid():
    # Check sample weighting raises errors.
    X = np.arange(100)[:, np.newaxis]
    y = np.ones(100)
    y[:50] = 0.0

    clf = DecisionTreeClassifier(random_state=0)

    sample_weight = np.random.rand(100, 1)
    assert_raises(ValueError, clf.fit, X, y, sample_weight=sample_weight)

    sample_weight = np.array(0)
    assert_raises(ValueError, clf.fit, X, y, sample_weight=sample_weight)

    sample_weight = np.ones(101)
    assert_raises(ValueError, clf.fit, X, y, sample_weight=sample_weight)

    sample_weight = np.ones(99)
    assert_raises(ValueError, clf.fit, X, y, sample_weight=sample_weight)


def check_class_weights(name):
    """Check class_weights resemble sample_weights behavior."""
    TreeClassifier = CLF_TREES[name]

    # Iris is balanced, so no effect expected for using 'balanced' weights
    clf1 = TreeClassifier(random_state=0)
    clf1.fit(iris.data, iris.target)
    clf2 = TreeClassifier(class_weight='balanced', random_state=0)
    clf2.fit(iris.data, iris.target)
    assert_almost_equal(clf1.feature_importances_, clf2.feature_importances_)

    # Make a multi-output problem with three copies of Iris
    iris_multi = np.vstack((iris.target, iris.target, iris.target)).T
    # Create user-defined weights that should balance over the outputs
    clf3 = TreeClassifier(class_weight=[{0: 2., 1: 2., 2: 1.},
                                        {0: 2., 1: 1., 2: 2.},
                                        {0: 1., 1: 2., 2: 2.}],
                          random_state=0)
    clf3.fit(iris.data, iris_multi)
    assert_almost_equal(clf2.feature_importances_, clf3.feature_importances_)
    # Check against multi-output "auto" which should also have no effect
    clf4 = TreeClassifier(class_weight='balanced', random_state=0)
    clf4.fit(iris.data, iris_multi)
    assert_almost_equal(clf3.feature_importances_, clf4.feature_importances_)

    # Inflate importance of class 1, check against user-defined weights
    sample_weight = np.ones(iris.target.shape)
    sample_weight[iris.target == 1] *= 100
    class_weight = {0: 1., 1: 100., 2: 1.}
    clf1 = TreeClassifier(random_state=0)
    clf1.fit(iris.data, iris.target, sample_weight)
    clf2 = TreeClassifier(class_weight=class_weight, random_state=0)
    clf2.fit(iris.data, iris.target)
    assert_almost_equal(clf1.feature_importances_, clf2.feature_importances_)

    # Check that sample_weight and class_weight are multiplicative
    clf1 = TreeClassifier(random_state=0)
    clf1.fit(iris.data, iris.target, sample_weight ** 2)
    clf2 = TreeClassifier(class_weight=class_weight, random_state=0)
    clf2.fit(iris.data, iris.target, sample_weight)
    assert_almost_equal(clf1.feature_importances_, clf2.feature_importances_)


def test_class_weights():
    for name in CLF_TREES:
        yield check_class_weights, name


def check_class_weight_errors(name):
    # Test if class_weight raises errors and warnings when expected.
    TreeClassifier = CLF_TREES[name]
    _y = np.vstack((y, np.array(y) * 2)).T

    # Invalid preset string
    clf = TreeClassifier(class_weight='the larch', random_state=0)
    assert_raises(ValueError, clf.fit, X, y)
    assert_raises(ValueError, clf.fit, X, _y)

    # Not a list or preset for multi-output
    clf = TreeClassifier(class_weight=1, random_state=0)
    assert_raises(ValueError, clf.fit, X, _y)

    # Incorrect length list for multi-output
    clf = TreeClassifier(class_weight=[{-1: 0.5, 1: 1.}], random_state=0)
    assert_raises(ValueError, clf.fit, X, _y)


def test_class_weight_errors():
    for name in CLF_TREES:
        yield check_class_weight_errors, name


def test_max_leaf_nodes():
    # Test greedy trees with max_depth + 1 leafs.
    from sklearn.tree._tree import TREE_LEAF
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    k = 4
    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(max_depth=None, max_leaf_nodes=k + 1).fit(X, y)
        tree = est.tree_
        assert_equal((tree.children_left == TREE_LEAF).sum(), k + 1)

        # max_leaf_nodes in (0, 1) should raise ValueError
        est = TreeEstimator(max_depth=None, max_leaf_nodes=0)
        assert_raises(ValueError, est.fit, X, y)
        est = TreeEstimator(max_depth=None, max_leaf_nodes=1)
        assert_raises(ValueError, est.fit, X, y)
        est = TreeEstimator(max_depth=None, max_leaf_nodes=0.1)
        assert_raises(ValueError, est.fit, X, y)


def test_max_leaf_nodes_max_depth():
    # Test precedence of max_leaf_nodes over max_depth.
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    k = 4
    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(max_depth=1, max_leaf_nodes=k).fit(X, y)
        tree = est.tree_
        assert_greater(tree.max_depth, 1)


def test_arrays_persist():
    # Ensure property arrays' memory stays alive when tree disappears
    # non-regression for #2726
    for attr in ['n_classes', 'value', 'children_left', 'children_right',
                 'threshold', 'impurity', 'feature', 'n_node_samples']:
        value = getattr(DecisionTreeClassifier().fit([[0], [1]], [0, 1]).tree_, attr)
        # if pointing to freed memory, contents may be arbitrary
        assert_true(-3 <= value.flat[0] < 3,
                    'Array points to arbitrary memory')


def test_only_constant_features():
    random_state = check_random_state(0)
    X = np.zeros((10, 20))
    y = random_state.randint(0, 2, (10, ))
    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(random_state=0)
        est.fit(X, y)
        assert_equal(est.tree_.max_depth, 0)


def test_behaviour_constant_feature_after_splits():
    X = np.transpose(np.vstack(([[0, 0, 0, 0, 0, 1, 2, 4, 5, 6, 7]],
                               np.zeros((4, 11)))))
    y = [0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3]
    for name, TreeEstimator in ALL_TREES.items():
        # do not check extra random trees
        if "ExtraTree" not in name:
            est = TreeEstimator(random_state=0, max_features=1)
            est.fit(X, y)
            assert_equal(est.tree_.max_depth, 2)
            assert_equal(est.tree_.node_count, 5)


def test_with_only_one_non_constant_features():
    X = np.hstack([np.array([[1.], [1.], [0.], [0.]]),
                   np.zeros((4, 1000))])

    y = np.array([0., 1., 0., 1.0])
    for name, TreeEstimator in CLF_TREES.items():
        est = TreeEstimator(random_state=0, max_features=1)
        est.fit(X, y)
        assert_equal(est.tree_.max_depth, 1)
        assert_array_equal(est.predict_proba(X), 0.5 * np.ones((4, 2)))

    for name, TreeEstimator in REG_TREES.items():
        est = TreeEstimator(random_state=0, max_features=1)
        est.fit(X, y)
        assert_equal(est.tree_.max_depth, 1)
        assert_array_equal(est.predict(X), 0.5 * np.ones((4, )))


def test_big_input():
    # Test if the warning for too large inputs is appropriate.
    X = np.repeat(10 ** 40., 4).astype(np.float64).reshape(-1, 1)
    clf = DecisionTreeClassifier()
    try:
        clf.fit(X, [0, 1, 0, 1])
    except ValueError as e:
        assert_in("float32", str(e))


def test_realloc():
    from sklearn.tree._utils import _realloc_test
    assert_raises(MemoryError, _realloc_test)


def test_huge_allocations():
    n_bits = 8 * struct.calcsize("P")

    X = np.random.randn(10, 2)
    y = np.random.randint(0, 2, 10)

    # Sanity check: we cannot request more memory than the size of the address
    # space. Currently raises OverflowError.
    huge = 2 ** (n_bits + 1)
    clf = DecisionTreeClassifier(splitter='best', max_leaf_nodes=huge)
    assert_raises(Exception, clf.fit, X, y)

    # Non-regression test: MemoryError used to be dropped by Cython
    # because of missing "except *".
    huge = 2 ** (n_bits - 1) - 1
    clf = DecisionTreeClassifier(splitter='best', max_leaf_nodes=huge)
    assert_raises(MemoryError, clf.fit, X, y)


def check_sparse_input(tree, dataset, max_depth=None):
    TreeEstimator = ALL_TREES[tree]
    X = DATASETS[dataset]["X"]
    X_sparse = DATASETS[dataset]["X_sparse"]
    y = DATASETS[dataset]["y"]

    # Gain testing time
    if dataset in ["digits", "boston"]:
        n_samples = X.shape[0] // 5
        X = X[:n_samples]
        X_sparse = X_sparse[:n_samples]
        y = y[:n_samples]

    for sparse_format in (csr_matrix, csc_matrix, coo_matrix):
        X_sparse = sparse_format(X_sparse)

        # Check the default (depth first search)
        d = TreeEstimator(random_state=0, max_depth=max_depth).fit(X, y)
        s = TreeEstimator(random_state=0, max_depth=max_depth).fit(X_sparse, y)

        assert_tree_equal(d.tree_, s.tree_,
                          "{0} with dense and sparse format gave different "
                          "trees".format(tree))

        y_pred = d.predict(X)
        if tree in CLF_TREES:
            y_proba = d.predict_proba(X)
            y_log_proba = d.predict_log_proba(X)

        for sparse_matrix in (csr_matrix, csc_matrix, coo_matrix):
            X_sparse_test = sparse_matrix(X_sparse, dtype=np.float32)

            assert_array_almost_equal(s.predict(X_sparse_test), y_pred)

            if tree in CLF_TREES:
                assert_array_almost_equal(s.predict_proba(X_sparse_test),
                                          y_proba)
                assert_array_almost_equal(s.predict_log_proba(X_sparse_test),
                                          y_log_proba)


def test_sparse_input():
    for tree_type, dataset in product(SPARSE_TREES, ("clf_small", "toy",
                                                     "digits", "multilabel",
                                                     "sparse-pos",
                                                     "sparse-neg",
                                                     "sparse-mix", "zeros")):
        max_depth = 3 if dataset == "digits" else None
        yield (check_sparse_input, tree_type, dataset, max_depth)

    # Due to numerical instability of MSE and too strict test, we limit the
    # maximal depth
    for tree_type, dataset in product(SPARSE_TREES, ["boston", "reg_small"]):
        if tree_type in REG_TREES:
            yield (check_sparse_input, tree_type, dataset, 2)


def check_sparse_parameters(tree, dataset):
    TreeEstimator = ALL_TREES[tree]
    X = DATASETS[dataset]["X"]
    X_sparse = DATASETS[dataset]["X_sparse"]
    y = DATASETS[dataset]["y"]

    # Check max_features
    d = TreeEstimator(random_state=0, max_features=1, max_depth=2).fit(X, y)
    s = TreeEstimator(random_state=0, max_features=1,
                      max_depth=2).fit(X_sparse, y)
    assert_tree_equal(d.tree_, s.tree_,
                      "{0} with dense and sparse format gave different "
                      "trees".format(tree))
    assert_array_almost_equal(s.predict(X), d.predict(X))

    # Check min_samples_split
    d = TreeEstimator(random_state=0, max_features=1,
                      min_samples_split=10).fit(X, y)
    s = TreeEstimator(random_state=0, max_features=1,
                      min_samples_split=10).fit(X_sparse, y)
    assert_tree_equal(d.tree_, s.tree_,
                      "{0} with dense and sparse format gave different "
                      "trees".format(tree))
    assert_array_almost_equal(s.predict(X), d.predict(X))

    # Check min_samples_leaf
    d = TreeEstimator(random_state=0,
                      min_samples_leaf=X_sparse.shape[0] // 2).fit(X, y)
    s = TreeEstimator(random_state=0,
                      min_samples_leaf=X_sparse.shape[0] // 2).fit(X_sparse, y)
    assert_tree_equal(d.tree_, s.tree_,
                      "{0} with dense and sparse format gave different "
                      "trees".format(tree))
    assert_array_almost_equal(s.predict(X), d.predict(X))

    # Check best-first search
    d = TreeEstimator(random_state=0, max_leaf_nodes=3).fit(X, y)
    s = TreeEstimator(random_state=0, max_leaf_nodes=3).fit(X_sparse, y)
    assert_tree_equal(d.tree_, s.tree_,
                      "{0} with dense and sparse format gave different "
                      "trees".format(tree))
    assert_array_almost_equal(s.predict(X), d.predict(X))


def test_sparse_parameters():
    for tree_type, dataset in product(SPARSE_TREES, ["sparse-pos",
                                                     "sparse-neg",
                                                     "sparse-mix", "zeros"]):
        yield (check_sparse_parameters, tree_type, dataset)


def check_sparse_criterion(tree, dataset):
    TreeEstimator = ALL_TREES[tree]
    X = DATASETS[dataset]["X"]
    X_sparse = DATASETS[dataset]["X_sparse"]
    y = DATASETS[dataset]["y"]

    # Check various criterion
    CRITERIONS = REG_CRITERIONS if tree in REG_TREES else CLF_CRITERIONS
    for criterion in CRITERIONS:
        d = TreeEstimator(random_state=0, max_depth=3,
                          criterion=criterion).fit(X, y)
        s = TreeEstimator(random_state=0, max_depth=3,
                          criterion=criterion).fit(X_sparse, y)

        assert_tree_equal(d.tree_, s.tree_,
                          "{0} with dense and sparse format gave different "
                          "trees".format(tree))
        assert_array_almost_equal(s.predict(X), d.predict(X))


def test_sparse_criterion():
    for tree_type, dataset in product(SPARSE_TREES, ["sparse-pos",
                                                     "sparse-neg",
                                                     "sparse-mix", "zeros"]):
        yield (check_sparse_criterion, tree_type, dataset)


def check_explicit_sparse_zeros(tree, max_depth=3,
                                n_features=10):
    TreeEstimator = ALL_TREES[tree]

    # n_samples set n_feature to ease construction of a simultaneous
    # construction of a csr and csc matrix
    n_samples = n_features
    samples = np.arange(n_samples)

    # Generate X, y
    random_state = check_random_state(0)
    indices = []
    data = []
    offset = 0
    indptr = [offset]
    for i in range(n_features):
        n_nonzero_i = random_state.binomial(n_samples, 0.5)
        indices_i = random_state.permutation(samples)[:n_nonzero_i]
        indices.append(indices_i)
        data_i = random_state.binomial(3, 0.5, size=(n_nonzero_i, )) - 1
        data.append(data_i)
        offset += n_nonzero_i
        indptr.append(offset)

    indices = np.concatenate(indices)
    data = np.array(np.concatenate(data), dtype=np.float32)
    X_sparse = csc_matrix((data, indices, indptr),
                          shape=(n_samples, n_features))
    X = X_sparse.toarray()
    X_sparse_test = csr_matrix((data, indices, indptr),
                               shape=(n_samples, n_features))
    X_test = X_sparse_test.toarray()
    y = random_state.randint(0, 3, size=(n_samples, ))

    # Ensure that X_sparse_test owns its data, indices and indptr array
    X_sparse_test = X_sparse_test.copy()

    # Ensure that we have explicit zeros
    assert_greater((X_sparse.data == 0.).sum(), 0)
    assert_greater((X_sparse_test.data == 0.).sum(), 0)

    # Perform the comparison
    d = TreeEstimator(random_state=0, max_depth=max_depth).fit(X, y)
    s = TreeEstimator(random_state=0, max_depth=max_depth).fit(X_sparse, y)

    assert_tree_equal(d.tree_, s.tree_,
                      "{0} with dense and sparse format gave different "
                      "trees".format(tree))

    Xs = (X_test, X_sparse_test)
    for X1, X2 in product(Xs, Xs):
        assert_array_almost_equal(s.tree_.apply(X1), d.tree_.apply(X2))
        assert_array_almost_equal(s.apply(X1), d.apply(X2))
        assert_array_almost_equal(s.apply(X1), s.tree_.apply(X1))

        assert_array_almost_equal(s.tree_.decision_path(X1).toarray(),
                                  d.tree_.decision_path(X2).toarray())
        assert_array_almost_equal(s.decision_path(X1).toarray(),
                                  d.decision_path(X2).toarray())
        assert_array_almost_equal(s.decision_path(X1).toarray(),
                                  s.tree_.decision_path(X1).toarray())

        assert_array_almost_equal(s.predict(X1), d.predict(X2))

        if tree in CLF_TREES:
            assert_array_almost_equal(s.predict_proba(X1),
                                      d.predict_proba(X2))


def test_explicit_sparse_zeros():
    for tree_type in SPARSE_TREES:
        yield (check_explicit_sparse_zeros, tree_type)


@ignore_warnings
def check_raise_error_on_1d_input(name):
    TreeEstimator = ALL_TREES[name]

    X = iris.data[:, 0].ravel()
    X_2d = iris.data[:, 0].reshape((-1, 1))
    y = iris.target

    assert_raises(ValueError, TreeEstimator(random_state=0).fit, X, y)

    est = TreeEstimator(random_state=0)
    est.fit(X_2d, y)
    assert_raises(ValueError, est.predict, [X])


@ignore_warnings
def test_1d_input():
    for name in ALL_TREES:
        yield check_raise_error_on_1d_input, name


def _check_min_weight_leaf_split_level(TreeEstimator, X, y, sample_weight):
    # Private function to keep pretty printing in nose yielded tests
    est = TreeEstimator(random_state=0)
    est.fit(X, y, sample_weight=sample_weight)
    assert_equal(est.tree_.max_depth, 1)

    est = TreeEstimator(random_state=0, min_weight_fraction_leaf=0.4)
    est.fit(X, y, sample_weight=sample_weight)
    assert_equal(est.tree_.max_depth, 0)


def check_min_weight_leaf_split_level(name):
    TreeEstimator = ALL_TREES[name]

    X = np.array([[0], [0], [0], [0], [1]])
    y = [0, 0, 0, 0, 1]
    sample_weight = [0.2, 0.2, 0.2, 0.2, 0.2]
    _check_min_weight_leaf_split_level(TreeEstimator, X, y, sample_weight)

    if not TreeEstimator().presort:
        _check_min_weight_leaf_split_level(TreeEstimator, csc_matrix(X), y,
                                           sample_weight)


def test_min_weight_leaf_split_level():
    for name in ALL_TREES:
        yield check_min_weight_leaf_split_level, name


def check_public_apply(name):
    X_small32 = X_small.astype(tree._tree.DTYPE)

    est = ALL_TREES[name]()
    est.fit(X_small, y_small)
    assert_array_equal(est.apply(X_small),
                       est.tree_.apply(X_small32))


def check_public_apply_sparse(name):
    X_small32 = csr_matrix(X_small.astype(tree._tree.DTYPE))

    est = ALL_TREES[name]()
    est.fit(X_small, y_small)
    assert_array_equal(est.apply(X_small),
                       est.tree_.apply(X_small32))


def test_public_apply():
    for name in ALL_TREES:
        yield (check_public_apply, name)

    for name in SPARSE_TREES:
        yield (check_public_apply_sparse, name)


def check_presort_sparse(est, X, y):
    assert_raises(ValueError, est.fit, X, y)


def test_presort_sparse():
    ests = (DecisionTreeClassifier(presort=True),
            DecisionTreeRegressor(presort=True))
    sparse_matrices = (csr_matrix, csc_matrix, coo_matrix)

    y, X = datasets.make_multilabel_classification(random_state=0,
                                                   n_samples=50,
                                                   n_features=1,
                                                   n_classes=20)
    y = y[:, 0]

    for est, sparse_matrix in product(ests, sparse_matrices):
        yield check_presort_sparse, est, sparse_matrix(X), y


def test_decision_path_hardcoded():
    X = iris.data
    y = iris.target
    est = DecisionTreeClassifier(random_state=0, max_depth=1).fit(X, y)
    node_indicator = est.decision_path(X[:2]).toarray()
    assert_array_equal(node_indicator, [[1, 1, 0], [1, 0, 1]])


def check_decision_path(name):
    X = iris.data
    y = iris.target
    n_samples = X.shape[0]

    TreeEstimator = ALL_TREES[name]
    est = TreeEstimator(random_state=0, max_depth=2)
    est.fit(X, y)

    node_indicator_csr = est.decision_path(X)
    node_indicator = node_indicator_csr.toarray()
    assert_equal(node_indicator.shape, (n_samples, est.tree_.node_count))

    # Assert that leaves index are correct
    leaves = est.apply(X)
    leave_indicator = [node_indicator[i, j] for i, j in enumerate(leaves)]
    assert_array_almost_equal(leave_indicator, np.ones(shape=n_samples))

    # Ensure only one leave node per sample
    all_leaves = est.tree_.children_left == TREE_LEAF
    assert_array_almost_equal(np.dot(node_indicator, all_leaves),
                              np.ones(shape=n_samples))

    # Ensure max depth is consistent with sum of indicator
    max_depth = node_indicator.sum(axis=1).max()
    assert_less_equal(est.tree_.max_depth, max_depth)


def test_decision_path():
    for name in ALL_TREES:
        yield (check_decision_path, name)


def check_no_sparse_y_support(name):
    X, y = X_multilabel, csr_matrix(y_multilabel)
    TreeEstimator = ALL_TREES[name]
    assert_raises(TypeError, TreeEstimator(random_state=0).fit, X, y)


def test_no_sparse_y_support():
    # Currently we don't support sparse y
    for name in ALL_TREES:
        yield (check_no_sparse_y_support, name)


def test_mae():
    # check MAE criterion produces correct results
    # on small toy dataset
    dt_mae = DecisionTreeRegressor(random_state=0, criterion="mae",
                                   max_leaf_nodes=2)
    dt_mae.fit([[3], [5], [3], [8], [5]], [6, 7, 3, 4, 3])
    assert_array_equal(dt_mae.tree_.impurity, [1.4, 1.5, 4.0/3.0])
    assert_array_equal(dt_mae.tree_.value.flat, [4, 4.5, 4.0])

    dt_mae.fit([[3], [5], [3], [8], [5]], [6, 7, 3, 4, 3],
               [0.6, 0.3, 0.1, 1.0, 0.3])
    assert_array_equal(dt_mae.tree_.impurity, [7.0/2.3, 3.0/0.7, 4.0/1.6])
    assert_array_equal(dt_mae.tree_.value.flat, [4.0, 6.0, 4.0])


def test_criterion_copy():
    # Let's check whether copy of our criterion has the same type
    # and properties as original
    n_outputs = 3
    n_classes = np.arange(3, dtype=np.intp)
    n_samples = 100

    def _pickle_copy(obj):
        return pickle.loads(pickle.dumps(obj))
    for copy_func in [copy.copy, copy.deepcopy, _pickle_copy]:
        for _, typename in CRITERIA_CLF.items():
            criteria = typename(n_outputs, n_classes)
            result = copy_func(criteria).__reduce__()
            typename_, (n_outputs_, n_classes_), _ = result
            assert_equal(typename, typename_)
            assert_equal(n_outputs, n_outputs_)
            assert_array_equal(n_classes, n_classes_)

        for _, typename in CRITERIA_REG.items():
            criteria = typename(n_outputs, n_samples)
            result = copy_func(criteria).__reduce__()
            typename_, (n_outputs_, n_samples_), _ = result
            assert_equal(typename, typename_)
            assert_equal(n_outputs, n_outputs_)
            assert_equal(n_samples, n_samples_)
