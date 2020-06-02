import numpy as np
import pytest
from sklearn import datasets
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree.tests.test_tree import REG_TREES, CLF_TREES


def test_montonic_constraints():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=0)
    train = np.arange(90)
    test = np.arange(90, 100)
    X_train = X[train]
    y_train = y[train]
    X_test_0 = np.copy(X[test])
    X_test_1 = np.copy(X_test_0)
    X_test_1[:, 0] += 10
    X_test_2 = np.copy(X_test_0)
    X_test_2[:, 1] += 10
    monotonic_cst = np.zeros(X.shape[1])
    monotonic_cst[0] = -1
    monotonic_cst[1] = 1

    for name, TreeRegressor in REG_TREES.items():
        est = TreeRegressor(max_depth=None, monotonic_cst=monotonic_cst)
        if hasattr(est, "random_state"):
            est.set_params(**{"random_state": 0})
        est.fit(X_train, y_train)

        y0 = est.predict(X_test_0)
        # decreasing constraint
        y1 = est.predict(X_test_1)
        # y1 should always be lower than y0
        assert(np.max(y1 - y0) <= 0)

        # increasing constraint
        y2 = est.predict(X_test_2)
        # y2 should always be greater than y0
        assert(np.min(y2 - y0) >= 0)

    for name, TreeClassifier in CLF_TREES.items():
        est = TreeClassifier(max_depth=None, monotonic_cst=monotonic_cst)
        if hasattr(est, "random_state"):
            est.set_params(**{"random_state": 0})
        est.fit(X_train, y_train)

        y0 = est.predict_proba(X_test_0)[:, 0]

        # decreasing constraint
        y1 = est.predict_proba(X_test_1)[:, 0]
        # y1 should always be lower than y0
        assert(np.max(y1 - y0) <= 0)

        # increasing constraint
        y2 = est.predict_proba(X_test_2)[:, 0]
        # y2 should always be greater than y0
        assert(np.min(y2 - y0) >= 0)


def test_multiclass_raises():
    X, y = datasets.make_hastie_10_2(n_samples=100, random_state=0)
    y[0] = 0
    monotonic_cst = np.zeros(X.shape[1])
    monotonic_cst[0] = -1
    monotonic_cst[1] = 1
    for name, TreeClassifier in CLF_TREES.items():
        est = TreeClassifier(max_depth=None, monotonic_cst=monotonic_cst,
                             random_state=0)
        if hasattr(est, "random_state"):
            est.set_params(**{"random_state": 0})

        with pytest.raises(ValueError):
            est.fit(X, y)


def is_monotonic(a, cst):
    return (cst * np.diff(a) >= 0.0).all()


def assert_children_values_monotonic_bounded(tree_, monotonic_cst):
    # Flip values so that only need to check for increasing constraint
    values = monotonic_cst * tree_.value

    for i in range(tree_.node_count):
        if tree_.children_left[i] > i and tree_.children_right[i] > i:
            # Check monotonicity
            i_left = tree_.children_left[i]
            i_right = tree_.children_right[i]
            assert(float(values[i_left]) <= float(values[i_right]))
            val_middle = float(values[i])
            # Check bounds
            if tree_.feature[i_left] >= 0:
                i_left_right = tree_.children_right[i_left]
                assert(float(values[i_left_right]) <= val_middle)
            if tree_.feature[i_right] >= 0:
                i_right_left = tree_.children_left[i_right]
                assert(val_middle <= float(values[i_right_left]))


def assert_tree_monotonic(clf, monotonic_cst):
    X_grid = np.arange(0, 1, 0.01).reshape(-1, 1)
    y_pred_grid = clf.predict(X_grid)
    assert is_monotonic(y_pred_grid, monotonic_cst)


@pytest.mark.parametrize('monotonic_cst', (-1, 1))
@pytest.mark.parametrize('splitter', ("best", "random"))
@pytest.mark.parametrize('depth_first', (True, False))
@pytest.mark.parametrize('seed', range(4))
def test_nodes_values(monotonic_cst, splitter, depth_first, seed):
    # Adaptation from test_nodes_values in test_montonic_constraints.py
    # Build a single tree with only one feature, and make sure the nodes
    # values respect the monotonic constraints.

    # Considering the following tree with a monotonic POS constraint, we
    # should have:
    #
    #       root
    #      /    \
    #     5      10
    #    / \    /  \
    #   a   b  c    d
    #
    # a <= b <= root <= c <= d (assert_children_values_monotonic_bounded)

    rng = np.random.RandomState(seed)
    n_samples = 1000
    n_features = 1
    X = rng.rand(n_samples, n_features)
    y = rng.rand(n_samples)

    if depth_first:
        # No max_leaf_nodes, default depth first tree builder
        clf = DecisionTreeRegressor(splitter=splitter,
                                    monotonic_cst=[monotonic_cst],
                                    random_state=seed)
    else:
        # max_leaf_nodes triggers depth first tree builder
        clf = DecisionTreeRegressor(splitter=splitter,
                                    monotonic_cst=[monotonic_cst],
                                    max_leaf_nodes=n_samples,
                                    random_state=seed)
    clf.fit(X, y)

    assert_children_values_monotonic_bounded(clf.tree_, monotonic_cst)
    assert_tree_monotonic(clf, monotonic_cst)
