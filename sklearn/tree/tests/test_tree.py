"""
Testing for the tree module (sklearn.tree).
"""

import pickle
import sys
sys.path.insert(0, "/Users/fares19/GitTools2/scikit-learn/")
import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, rand
from functools import partial
from itertools import product

from sklearn.metrics import accuracy_score
from sklearn.metrics import mean_squared_error

from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import raises
from sklearn.utils.validation import check_random_state

from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree import ExtraTreeClassifier
from sklearn.tree import ExtraTreeRegressor

from sklearn import tree
from sklearn.tree.tree import SPARSE_SPLITTER
from sklearn import datasets
from sklearn.ensemble import RandomForestClassifier
import time

from sklearn.preprocessing._weights import _balance_weights


CLF_CRITERIONS = ("gini", "entropy")
REG_CRITERIONS = ("mse", )

CLF_TREES = {
    "DecisionTreeClassifier": DecisionTreeClassifier,
    "Presort-DecisionTreeClassifier": partial(DecisionTreeClassifier,
                                              splitter="presort-best"),
    "ExtraTreeClassifier": ExtraTreeClassifier,
}

REG_TREES = {
    "DecisionTreeRegressor": DecisionTreeRegressor,
    "Presort-DecisionTreeRegressor": partial(DecisionTreeRegressor,
                                             splitter="presort-best"),
    "ExtraTreeRegressor": ExtraTreeRegressor,
}

ALL_TREES = dict()
ALL_TREES.update(CLF_TREES)
ALL_TREES.update(REG_TREES)


X_small = [[0, 0, 4, 0, 0, 0, 1, -14, 0, -4, 0, 0, 0, 0, ],
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
          [3, 1, 0, 1, 0, -4, 1, 0, 1, -2, 0, 0, 1, 0, ]]

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


def assert_tree_equal(d, s, message):
    from sklearn.tree._tree import TREE_LEAF
    assert_equal(s.node_count, d.node_count,
                 "{0}: inequal number of node ({1} != {2})"
                 "".format(message, s.node_count, d.node_count))

    assert_array_equal(d.children_right, s.children_right,
                       message + ": inequal n_node_samples")
    assert_array_equal(d.children_left, s.children_left,
                       message + ": inequal children_left")

    external = d.children_right == TREE_LEAF
    internal = np.logical_not(external)

    assert_array_equal(d.feature[internal], s.feature[internal],
                       message + ": inequal features")
    assert_array_equal(d.threshold[internal], s.threshold[internal],
                       message + ": inequal threshold")
    assert_array_equal(d.n_node_samples, s.n_node_samples,
                       message + ": inequal n_node_samples")

    assert_almost_equal(d.impurity, s.impurity,
                        err_msg=message + ": inequal impurity")

    assert_array_almost_equal(d.value[external], s.value[external],
                              err_msg=message + ": inequal value")


def test_classification_toy():
    """Check classification on a toy dataset."""
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
    """Check classification on a weighted toy dataset."""
    for name, Tree in CLF_TREES.items():
        clf = Tree(random_state=0)

        clf.fit(X, y, sample_weight=np.ones(len(X)))
        assert_array_equal(clf.predict(T), true_result,
                           "Failed with {0}".format(name))

        clf.fit(X, y, sample_weight=np.ones(len(X)) * 0.5)
        assert_array_equal(clf.predict(T), true_result,
                           "Failed with {0}".format(name))


def test_regression_toy():
    """Check regression on a toy dataset."""
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
    """Check on a XOR problem"""
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
    """Check consistency on dataset iris."""
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
    """Check consistency on dataset boston house prices."""

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
    """Predict probabilities using DecisionTreeClassifier."""

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
    """Check the array representation."""
    # Check resize
    X = np.arange(10000)[:, np.newaxis]
    y = np.arange(10000)

    for name, Tree in REG_TREES.items():
        reg = Tree(max_depth=None, random_state=0)
        reg.fit(X, y)


def test_pure_set():
    """Check when y is pure."""
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
        assert_almost_equal(clf.predict(X), y,
                            err_msg="Failed with {0}".format(name))


def test_numerical_stability():
    """Check numerical stability."""
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
    """Check variable importances."""
    X, y = datasets.make_classification(n_samples=2000,
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

        X_new = clf.transform(X, threshold="mean")
        assert_less(0, X_new.shape[1], "Failed with {0}".format(name))
        assert_less(X_new.shape[1], X.shape[1], "Failed with {0}".format(name))


@raises(ValueError)
def test_importances_raises():
    """Check if variable importance before fit raises ValueError. """
    clf = DecisionTreeClassifier()
    clf.feature_importances_


def test_importances_gini_equal_mse():
    """Check that gini is equivalent to mse for binary output variable"""

    X, y = datasets.make_classification(n_samples=2000,
                                        n_features=10,
                                        n_informative=3,
                                        n_redundant=0,
                                        n_repeated=0,
                                        shuffle=False,
                                        random_state=0)

    clf = DecisionTreeClassifier(criterion="gini", random_state=0).fit(X, y)
    reg = DecisionTreeRegressor(criterion="mse", random_state=0).fit(X, y)

    assert_almost_equal(clf.feature_importances_, reg.feature_importances_)
    assert_array_equal(clf.tree_.feature, reg.tree_.feature)
    assert_array_equal(clf.tree_.children_left, reg.tree_.children_left)
    assert_array_equal(clf.tree_.children_right, reg.tree_.children_right)
    assert_array_equal(clf.tree_.n_node_samples, reg.tree_.n_node_samples)


def test_max_features():
    """Check max_features."""
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
    """Test that it gives proper exception on deficient input."""
    for name, TreeEstimator in CLF_TREES.items():
        # predict before fit
        est = TreeEstimator()
        assert_raises(Exception, est.predict_proba, X)

        est.fit(X, y)
        X2 = [-2, -1, 1]  # wrong feature shape for sample
        assert_raises(ValueError, est.predict_proba, X2)

    for name, TreeEstimator in ALL_TREES.items():
        # Invalid values for parameters
        assert_raises(ValueError, TreeEstimator(min_samples_leaf=-1).fit, X, y)
        assert_raises(ValueError, TreeEstimator(min_samples_split=-1).fit,
                      X, y)
        assert_raises(ValueError, TreeEstimator(max_depth=-1).fit, X, y)
        assert_raises(ValueError, TreeEstimator(max_features=42).fit, X, y)

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
        assert_raises(Exception, est.predict, T)

        # predict on vector with different dims
        est.fit(X, y)
        t = np.asarray(T)
        assert_raises(ValueError, est.predict, t[:, 1:])

        # wrong sample shape
        Xt = np.array(X).T

        est = TreeEstimator()
        est.fit(np.dot(X, Xt), y)
        assert_raises(ValueError, est.predict, X)

        clf = TreeEstimator()
        clf.fit(X, y)
        assert_raises(ValueError, clf.predict, Xt)


def test_min_samples_leaf():
    """Test if leaves contain more than leaf_count training examples"""
    X = np.asfortranarray(iris.data.astype(tree._tree.DTYPE))
    y = iris.target

    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(min_samples_leaf=5, random_state=0)
        est.fit(X, y)
        out = est.tree_.apply(X)
        node_counts = np.bincount(out)
        leaf_count = node_counts[node_counts != 0]  # drop inner nodes
        assert_greater(np.min(leaf_count), 4,
                       "Failed with {0}".format(name))


def test_pickle():
    """Check that tree estimator are pickable """
    for name, TreeClassifier in CLF_TREES.items():
        clf = TreeClassifier(random_state=0)
        clf.fit(iris.data, iris.target)
        score = clf.score(iris.data, iris.target)

        serialized_object = pickle.dumps(clf)
        clf2 = pickle.loads(serialized_object)
        assert_equal(type(clf2), clf.__class__)
        score2 = clf2.score(iris.data, iris.target)
        assert_equal(score, score2, "Failed to generate same score "
                                    "after pickling (classification) "
                                    "with {0}".format(name))

    for name, TreeRegressor in REG_TREES.items():
        reg = TreeRegressor(random_state=0)
        reg.fit(boston.data, boston.target)
        score = reg.score(boston.data, boston.target)

        serialized_object = pickle.dumps(reg)
        reg2 = pickle.loads(serialized_object)
        assert_equal(type(reg2), reg.__class__)
        score2 = reg2.score(boston.data, boston.target)
        assert_equal(score, score2, "Failed to generate same score "
                                    "after pickling (regression) "
                                    "with {0}".format(name))


def test_multioutput():
    """Check estimators on multi-output problems."""
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
    """Test that n_classes_ and classes_ have proper shape."""
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
    """Check class rebalancing."""
    unbalanced_X = iris.data[:125]
    unbalanced_y = iris.target[:125]
    sample_weight = _balance_weights(unbalanced_y)

    for name, TreeClassifier in CLF_TREES.items():
        clf = TreeClassifier(random_state=0)
        clf.fit(unbalanced_X, unbalanced_y, sample_weight=sample_weight)
        assert_almost_equal(clf.predict(unbalanced_X), unbalanced_y)


def test_memory_layout():
    """Check that it works no matter the memory layout"""
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

        if est.splitter in SPARSE_SPLITTER:
            # csr matrix
            X = csr_matrix(iris.data)
            y = iris.target
            assert_array_equal(est.fit(X, y).predict(X), y)

            # csc_matrix
            X = csc_matrix(iris.data)
            y = iris.target
            assert_array_equal(est.fit(X, y).predict(X), y)

        # Strided
        X = np.asarray(iris.data[::3], dtype=dtype)
        y = iris.target[::3]
        assert_array_equal(est.fit(X, y).predict(X), y)


def test_sample_weight():
    """Check sample weighting."""
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

    sample_weight[y == 2] = .50  # Samples of class '2' are no longer weightier
    clf = DecisionTreeClassifier(max_depth=1, random_state=0)
    clf.fit(X, y, sample_weight=sample_weight)
    assert_equal(clf.tree_.threshold[0], 49.5)  # Threshold should have moved

    # Test that sample weighting is the same as having duplicates
    X = iris.data
    y = iris.target

    duplicates = rng.randint(0, X.shape[0], 200)

    clf = DecisionTreeClassifier(random_state=1)
    clf.fit(X[duplicates], y[duplicates])

    sample_weight = np.bincount(duplicates, minlength=X.shape[0])
    clf2 = DecisionTreeClassifier(random_state=1)
    clf2.fit(X, y, sample_weight=sample_weight)

    internal = clf.tree_.children_left != tree._tree.TREE_LEAF
    assert_array_almost_equal(clf.tree_.threshold[internal],
                              clf2.tree_.threshold[internal])


def test_sample_weight_invalid():
    """Check sample weighting raises errors."""
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


def test_max_leaf_nodes():
    """Test greedy trees with max_depth + 1 leafs. """
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
    """Test preceedence of max_leaf_nodes over max_depth. """
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=1)
    k = 4
    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(max_depth=1, max_leaf_nodes=k).fit(X, y)
        tree = est.tree_
        assert_greater(tree.max_depth, 1)


def test_arrays_persist():
    """Ensure property arrays' memory stays alive when tree disappears

    non-regression for #2726
    """
    for attr in ['n_classes', 'value', 'children_left', 'children_right',
                 'threshold', 'impurity', 'feature', 'n_node_samples']:
        value = getattr(DecisionTreeClassifier().fit([[0]], [0]).tree_, attr)
        # if pointing to freed memory, contents may be arbitrary
        assert_true(-2 <= value.flat[0] < 2,
                    'Array points to arbitrary memory')


def test_only_constant_features():
    random_state = check_random_state(0)
    X = np.zeros((10, 20))
    y = random_state.randint(0, 2, (10, ))
    for name, TreeEstimator in ALL_TREES.items():
        est = TreeEstimator(random_state=0)
        est.fit(X, y)
        assert_equal(est.tree_.max_depth, 0)


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
    """Test if the warning for too large inputs is appropriate."""
    X = np.repeat(10 ** 40., 4).astype(np.float64).reshape(-1, 1)
    clf = DecisionTreeClassifier()
    try:
        clf.fit(X, [0, 1, 0, 1])
    except ValueError as e:
        assert_in("float32", str(e))


# def test_memoryerror():
#     from sklearn.tree._tree import _realloc_test
#     assert_raises(MemoryError, _realloc_test)


def test_sparse_input_on_x_small():

    for max_features in [1, None]:
        s = DecisionTreeClassifier(random_state=0, max_features=max_features)
        s.fit(csc_matrix(X_small), y_small)

        d = DecisionTreeClassifier(random_state=0, max_features=max_features)
        d.fit(X_small, y_small)

        assert_tree_equal(d.tree_, s.tree_,
                          "dense and sparse format gave different trees")
        assert_array_almost_equal(s.predict(X_small), d.predict(X_small))

    for max_features in [1, None]:
        s = DecisionTreeRegressor(random_state=0, max_features=max_features)
        s.fit(csc_matrix(X_small), y_small)

        d = DecisionTreeRegressor(random_state=0, max_features=max_features)
        d.fit(X_small, y_small)

        assert_tree_equal(d.tree_, s.tree_,
                          "dense and sparse format gave different trees")
        assert_array_almost_equal(s.predict(X_small), d.predict(X_small))


def test_sparse_input_on_random_data():
    rng = np.random.RandomState(0)
    X = rng.rand(40, 10)
    X[X < .8] = 0
    X = csr_matrix(X)
    y = (4 * rng.rand(40)).astype(np.int)

    for name, TreeEstimator in ALL_TREES.items():

        d = TreeEstimator(random_state=0).fit(X.toarray(), y)
        if d.splitter not in SPARSE_SPLITTER:
            continue

        s = TreeEstimator(random_state=0).fit(X, y)

        assert_tree_equal(d.tree_, s.tree_,
                          "{0} with dense and sparse format gave different "
                          "trees".format(name))
        assert_array_almost_equal(s.predict(X), d.predict(X))


def test_sparse_input_boston():
    for max_leaf_nodes in [None, 10, 20]:
        d = DecisionTreeClassifier(random_state=0,
                                   max_leaf_nodes=max_leaf_nodes)
        d.fit(csc_matrix(boston.data), boston.target)

        # csc data
        csc = DecisionTreeClassifier(random_state=0,
                                     max_leaf_nodes=max_leaf_nodes)
        csc.fit(csc_matrix(boston.data), boston.target)
        assert_tree_equal(d.tree_, csc.tree_,
                          "dense and sparse format gave different trees")
        assert_array_almost_equal(csc.predict(boston.data),
                                  d.predict(boston.data))

        # csr data
        csr = DecisionTreeClassifier(random_state=0,
                                     max_leaf_nodes=max_leaf_nodes)
        csr.fit(csc_matrix(boston.data), boston.target)

        assert_tree_equal(d.tree_, csr.tree_,
                          "dense and sparse format gave different trees")
        assert_array_almost_equal(csr.predict(boston.data),
                                  d.predict(boston.data))


def test_sparse_input_digits():
    for max_leaf_nodes in [None, 10, 20]:
        d = DecisionTreeClassifier(random_state=0,
                                   max_leaf_nodes=max_leaf_nodes)
        d.fit(digits.data, digits.target)

        # csc data
        csc = DecisionTreeClassifier(random_state=0,
                                     max_leaf_nodes=max_leaf_nodes)
        csc.fit(csc_matrix(digits.data), digits.target)
        assert_tree_equal(d.tree_, csc.tree_,
                          "dense and sparse format gave different trees")
        assert_array_almost_equal(csc.predict(csc_matrix(digits.data)),
                                  d.predict(digits.data))

        # csr data
        csr = DecisionTreeClassifier(random_state=0,
                                     max_leaf_nodes=max_leaf_nodes)
        csr.fit(csr_matrix(digits.data), digits.target)

        assert_tree_equal(d.tree_, csr.tree_,
                          "dense and sparse format gave different trees")
        assert_array_almost_equal(csr.predict(csr_matrix(digits.data)),
                                  d.predict(digits.data))


def test_sparse_with_various_criterion():
    B = datasets.load_digits()
    X_ = B.data
    y_ = B.target
    d = DecisionTreeClassifier(random_state=0,
                               criterion="entropy").fit(X_, y_).tree_
    s = DecisionTreeClassifier(random_state=0,
                               criterion="entropy").fit(csc_matrix(X_),
                                                        y_).tree_
    message = 'Sparse & Dense Trees are not the same, fitting the digits data'
    assert_tree_equal(d, s, message)


def test_equality_of_sparse_and_dense_tree_with_digits_max_depth():

    B = datasets.load_digits()
    X_ = B.data
    y_ = B.target
    d = DecisionTreeClassifier(random_state=0, max_depth=1).fit(X_, y_).tree_
    s = DecisionTreeClassifier(random_state=0,
                               max_depth=1).fit(csc_matrix(X_), y_).tree_
    message = 'Sparse & Dense Trees are not the same, fitting the digits data'
    assert_tree_equal(d, s, message)


def test_equality_of_sparse_and_dense_tree_with_digits_min_samples_split():

    B = datasets.load_digits()
    X_ = B.data
    y_ = B.target
    d = DecisionTreeClassifier(random_state=0,
                               min_samples_split=4).fit(X_, y_).tree_
    s = DecisionTreeClassifier(random_state=0,
                               min_samples_split=4).fit(csc_matrix(X_),
                                                        y_).tree_
    message = 'Sparse & Dense Trees are not the same, fitting the digits data'
    assert_tree_equal(d, s, message)


def test_equality_of_sparse_and_dense_tree_with_digits_min_samples_leaf():

    B = datasets.load_digits()
    X_ = B.data
    y_ = B.target
    d = DecisionTreeClassifier(random_state=0,
                               min_samples_leaf=10).fit(X_, y_).tree_
    s = DecisionTreeClassifier(random_state=0,
                               min_samples_leaf=10).fit(csc_matrix(X_),
                                                        y_).tree_
    message = 'Sparse & Dense Trees are not the same, fitting the digits data'
    assert_tree_equal(d, s, message)


def test_random_sparse_matrix_best_first_search():
    n_samples = 40
    n_features = 20

    n_test = 20
    for dens in [0.0, 0.3, 1.0]:
        X_ = rand(n_samples, n_features, density=dens, format='csc')
        y_ = np.random.randint(2, size=n_samples)
        X_test = rand(n_test, n_features, density=dens, format='csr')

        s = DecisionTreeClassifier(random_state=0,
                                   max_leaf_nodes=5).fit(X_, y_)

        d = DecisionTreeClassifier(random_state=0,
                                   max_leaf_nodes=5).fit(X_.toarray(), y_)

        assert_array_equal(s.predict(X_test),
                           d.predict(X_test))

        assert_array_equal(s.predict_proba(X_test.toarray()),
                           d.predict_proba(X_test.toarray()))

        assert_tree_equal(d.tree_, s.tree_,
                          "Sparse and Dense Trees are not the same,"
                          "fitting the random data")


def test_random_sparse_matrix_depth_first_search():
    n_samples = 20
    n_features = 20
    n_test = 100
    for density in [0.01, 0.1, 0.5]:
        X_ = rand(n_samples, n_features, density=density, format='csc')
        y_ = np.random.randint(2, size=n_samples)
        X_test = rand(n_test, n_features, density=density, format='csr')

        s = DecisionTreeClassifier(random_state=0,
                                   max_depth=100).fit(X_, y_)
        d = DecisionTreeClassifier(random_state=0,
                                   max_depth=100).fit(X_.toarray(), y_)

        assert_tree_equal(d.tree_, s.tree_,
                          "Sparse and Dense Trees are not the same,"
                          "fitting the random data")

        assert_array_equal(s.predict(X_test),
                           d.predict(X_test))

        assert_array_equal(s.predict_proba(X_test.toarray()),
                           d.predict_proba(X_test.toarray()))


def test_random_sparse_matrix_best_first_search_reg():

    n_samples = 40
    n_features = 20

    n_test = 20
    for dens in [0.0, 0.3, 1.0]:
        X_ = rand(n_samples, n_features, density=dens, format='csc')
        y_ = np.random.randint(2, size=n_samples)
        X_test = rand(n_test, n_features, density=dens, format='csr')

        s = DecisionTreeRegressor(random_state=0,
                                  max_leaf_nodes=5).fit(X_, y_)

        d = DecisionTreeRegressor(random_state=0,
                                  max_leaf_nodes=5).fit(X_.toarray(), y_)

        assert_array_equal(s.predict(X_test),
                           d.predict(X_test))

        assert_tree_equal(d.tree_, s.tree_,
                          "Sparse and Dense Trees are not the same, "
                          "fitting the random data")


def test_random_sparse_matrix_depth_first_search_reg():

    n_samples = 40
    n_features = 20

    n_test = 20
    for dens in [0.0, 0.3, 1.0]:
        X_ = rand(n_samples, n_features, density=dens, format='csc')
        y_ = np.random.randint(2, size=n_samples)
        X_test = rand(n_test, n_features, density=dens, format='csr')

        s = DecisionTreeRegressor(random_state=0,
                                  max_depth=100).fit(X_, y_)

        d = DecisionTreeRegressor(random_state=0,
                                  max_depth=100).fit(X_.toarray(), y_)

        assert_array_equal(s.predict(X_test),
                           d.predict(X_test))

        assert_tree_equal(d.tree_, s.tree_,
                          "Sparse and Dense Trees are not the same,"
                          "fitting the random data")


def test_random_sparse_matrix_depth_first_search_negative_input():
    n_samples = 40
    n_features = 20

    n_test = 20
    for dens in [0.0, 0.3, 1.0]:
        X_ = -1*rand(n_samples, n_features, density=dens, format='csc')
        y_ = np.random.randint(2, size=n_samples)
        X_test = -1*rand(n_test, n_features, density=dens, format='csr')

        s = DecisionTreeClassifier(random_state=0,
                                   max_depth=100).fit(X_, y_)

        d = DecisionTreeClassifier(random_state=0,
                                   max_depth=100).fit(X_.toarray(), y_)

        assert_array_equal(s.predict(X_test),
                           d.predict(X_test))

        assert_array_equal(s.predict_proba(X_test.toarray()),
                           d.predict_proba(X_test.toarray()))

        assert_tree_equal(d.tree_, s.tree_,
                          "Sparse and Dense Trees are not the same,"
                          "fitting the random data")


def test_random_sparse_matrix_depth_first_search_mixed_input():
    n_samples = 40
    n_features = 20

    n_test = 20
    X_ = rand(n_samples, n_features, density=0.3, format='csc')
    X_.data = np.array([(np.random.uniform(0, 2)*2-1)*x for x in X_.data])

    X_test = rand(n_test, n_features, density=0.3, format='csr')
    X_test.data = np.array([(np.random.uniform(0, 2)*2-1)*x
                            for x in X_test.data])

    y_ = np.random.randint(2, size=n_samples)
    s = DecisionTreeClassifier(random_state=0,
                               max_depth=100).fit(X_, y_)

    d = DecisionTreeClassifier(random_state=0,
                               max_depth=100).fit(X_.toarray(), y_)

    assert_array_equal(s.predict(X_test),
                       d.predict(X_test))

    assert_array_equal(s.predict_proba(X_test.toarray()),
                       d.predict_proba(X_test.toarray()))

    assert_tree_equal(d.tree_, s.tree_,
                      "Sparse and Dense Trees are not the same,"
                      "fitting the random data")


def test_multi_label_classification():
    n_samples = 40
    n_features = 20

    n_test = 20
    X_ = rand(n_samples, n_features, density=0.3, format='csc')
    X_.data = np.array([(np.random.uniform(0, 2)*2-1)*x for x in X_.data])

    X_test = rand(n_test, n_features, density=0.3, format='csr')
    X_test.data = np.array([(np.random.uniform(0, 2)*2-1)*x
                            for x in X_test.data])

    y_ = np.random.randint(2, size=(n_samples, 3))
    s = DecisionTreeClassifier(random_state=0,
                               max_depth=100).fit(X_, y_)

    d = DecisionTreeClassifier(random_state=0,
                               max_depth=100).fit(X_.toarray(), y_)

    assert_array_equal(s.predict(X_test),
                       d.predict(X_test))

    assert_array_equal(s.predict_proba(X_test.toarray()),
                       d.predict_proba(X_test.toarray()))

    assert_tree_equal(d.tree_, s.tree_,
                      "Sparse and Dense Trees are not the same,"
                      "fitting the random data")


def test_random_sparse_matrix_of_random_forest():
    n_samples = 40
    n_features = 20

    n_test = 20
    for dens in [0.1, 0.3, 1.0]:
        X_ = rand(n_samples, n_features, density=dens, format='csc')
        y_ = np.random.randint(2, size=n_samples)
        X_test = rand(n_test, n_features, density=dens, format='csr')

        s = RandomForestClassifier(n_estimators=2, random_state=0,
                                   max_depth=100).fit(X_, y_)

        d = RandomForestClassifier(n_estimators=2, random_state=0,
                                   max_depth=100).fit(X_.toarray(), y_)

        Ss = s.estimators_
        Ds = d.estimators_
        for i in range(len(Ds)):
            assert_tree_equal(Ds[i].tree_, Ss[i].tree_,
                              "Sparse and Dense Random Forests not the same")

        assert_array_equal(s.predict(X_test.toarray()),
                           d.predict(X_test.toarray()))

        assert_array_almost_equal(s.predict_proba(X_test.toarray()),
                                  d.predict_proba(X_test.toarray()))