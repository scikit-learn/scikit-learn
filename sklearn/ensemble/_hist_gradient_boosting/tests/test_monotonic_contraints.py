import numpy as np
import pytest

from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower
from sklearn.ensemble._hist_gradient_boosting.common import G_H_DTYPE
from sklearn.experimental import enable_hist_gradient_boosting  # noqa
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.ensemble._hist_gradient_boosting.common import MonotonicConstraint


def is_increasing(a):
    return (np.diff(a) >= 0.0).all()


def is_decreasing(a):
    return (np.diff(a) <= 0.0).all()


def assert_leaves_values_monotonic(predictor, monotonic_cst):
    # make sure leaves values (from left to right) are either all increasing
    # or all decreasing (or neither) depending on the monotonic constraint.
    nodes = predictor.nodes

    def get_leaves_values():
        """get leaves values from left to right"""
        values = []

        def dfs(node_idx):
            node = nodes[node_idx]
            if node['is_leaf']:
                values.append(node['value'])
                return
            dfs(node['left'])
            dfs(node['right'])

        dfs(0)  # start at root (0)
        return values

    values = get_leaves_values()

    if monotonic_cst == MonotonicConstraint.NO_CST:
        # some increasing, some decreasing
        assert not is_increasing(values) and not is_decreasing(values)
    elif monotonic_cst == MonotonicConstraint.INC:
        # all increasing
        assert is_increasing(values)
    else:  # DEC
        # all decreasing
        assert is_decreasing(values)


def assert_children_values_monotonic(predictor, monotonic_cst):
    # Make sure siblings values respect the monotonic constraints. Left should
    # be lower (resp greater) than right child if constraint is INC (resp.
    # DEC).
    # Note that this property alone isn't enough to ensure full monotonicity,
    # since we also need to guanrantee that all the descendents of the left
    # child won't be greater (resp. lower) than the right child, or its
    # descendents. That's why we need to bound the predicted values (this is
    # tested in assert_children_values_bounded)
    nodes = predictor.nodes

    left_lower = []
    left_greater = []

    def dfs(node_idx):
        node = nodes[node_idx]
        if node['is_leaf']:
            return

        left_idx = node['left']
        right_idx = node['right']

        if nodes[left_idx]['value'] < nodes[right_idx]['value']:
            left_lower.append(node)
        elif nodes[left_idx]['value'] > nodes[right_idx]['value']:
            left_greater.append(node)
        dfs(left_idx)
        dfs(right_idx)

    dfs(0)  # start at root (0)

    if monotonic_cst == MonotonicConstraint.NO_CST:
        assert left_lower and left_greater
    elif monotonic_cst == MonotonicConstraint.INC:
        assert left_lower and not left_greater
    else:  # DEC
        assert not left_lower and left_greater


def assert_children_values_bounded(grower, monotonic_cst):
    # Make sure that the values of the children of a node are bounded by the
    # middle value between that node and its sibling (if there is a monotonic
    # constraint).
    # As a bonus, we also check that the siblings values are properly ordered
    # which is slightly redundant with assert_children_values_monotonic (but
    # this check is done on the grower nodes whereas
    # assert_children_values_monotonic is done on the predictor nodes)

    if monotonic_cst == MonotonicConstraint.NO_CST:
        return

    def dfs(node):
        if node.is_leaf:
            return
        if node is not grower.root and node is node.parent.left_child:
            sibling = node.sibling  # on the right
            middle = (node.value + sibling.value) / 2
            if monotonic_cst == MonotonicConstraint.INC:
                assert (node.left_child.value <=
                        node.right_child.value <=
                        middle)
                if not sibling.is_leaf:
                    assert (middle <=
                            sibling.left_child.value <=
                            sibling.right_child.value)
            else:  # DEC
                assert (node.left_child.value >=
                        node.right_child.value >=
                        middle)
                if not sibling.is_leaf:
                    assert (middle >=
                            sibling.left_child.value >=
                            sibling.right_child.value)

        dfs(node.left_child)
        dfs(node.right_child)

    dfs(grower.root)


@pytest.mark.parametrize('seed', range(3))
@pytest.mark.parametrize('monotonic_cst', (
    MonotonicConstraint.NO_CST,
    MonotonicConstraint.INC,
    MonotonicConstraint.DEC,
))
def test_nodes_values(monotonic_cst, seed):
    # Build a single tree with only one feature, and make sure the nodes
    # values respect the monotonic constraints.

    # Considering the following tree with a monotonic INC constraint, we
    # should have:
    #
    #       root
    #      /    \
    #     5     10    # middle = 7.5
    #    / \   / \
    #   a  b  c  d
    #
    # a < b and c < d  (assert_children_values_monotonic)
    # a, b < middle < c, d (assert_children_values_bounded)
    # a < b < c < d (assert_leaves_values_monotonic)
    #
    # The last one is a consequence of the others, but can't hurt to check

    rng = np.random.RandomState(seed)
    n_samples = 1000
    n_features = 1
    X_binned = rng.randint(0, 256, size=(n_samples, n_features),
                           dtype=np.uint8)
    X_binned = np.asfortranarray(X_binned)

    gradients = rng.normal(size=n_samples).astype(G_H_DTYPE)
    hessians = np.ones(shape=1, dtype=G_H_DTYPE)

    grower = TreeGrower(X_binned, gradients, hessians,
                        monotonic_cst=[monotonic_cst],
                        shrinkage=.1)
    grower.grow()

    # grow() will shrink the leaves values at the very end. For our comparison
    # tests, we need to revert the shrinkage of the leaves, else we would
    # compare the value of a leaf (shrunk) with a node (not shrunk) and the
    # test would not be correct.
    for leave in grower.finalized_leaves:
        leave.value /= grower.shrinkage

    predictor = grower.make_predictor()
    assert_children_values_monotonic(predictor, monotonic_cst)
    assert_children_values_bounded(grower, monotonic_cst)
    assert_leaves_values_monotonic(predictor, monotonic_cst)


@pytest.mark.parametrize('seed', range(3))
def test_predictions(seed):
    # Train a model with an INC constraint on the first feature and a DEC
    # constraint on the second feature, and make sure the constraints are
    # respected by checking the predictions.
    # test adapted from lightgbm's test_monotone_constraint(), itself inspired
    # by https://xgboost.readthedocs.io/en/latest//tutorials/monotonic.html

    rng = np.random.RandomState(seed)

    n_samples = 1000
    f_0 = rng.random(size=n_samples)  # positive correlation with y
    f_1 = rng.random(size=n_samples)  # negative correslation with y
    X = np.c_[f_0, f_1]
    noise = rng.normal(loc=0.0, scale=0.01, size=n_samples)
    y = (5 * f_0 + np.sin(10 * np.pi * f_0) -
         5 * f_1 - np.cos(10 * np.pi * f_1) +
         noise)

    gbdt = HistGradientBoostingRegressor()
    gbdt._monotonic_cst = [1, -1]  # INC, DEC
    gbdt.fit(X, y)

    linspace = np.linspace(0, 1, 100)
    sin = np.sin(linspace)
    constant = np.full_like(linspace, fill_value=.5)

    # We now assert the predictions properly respect the constraints, on each
    # feature. When testing for a feature we need to set the other one to a
    # constant, because the monotonic constraints are only a "all else being
    # equal" type of constraints:
    # a constraint on the first feature only means that
    # x0 < x0' => f(x0, x1) < f(x0', x1)
    # while x1 stays constant.
    # The constraint does not guanrantee that
    # x0 < x0' => f(x0, x1) < f(x0', x1')

    # First feature
    # assert pred is all increasing when f_0 is all increasing
    X = np.c_[linspace, constant]
    pred = gbdt.predict(X)
    assert is_increasing(pred)
    # assert pred actually follows the variations of f_0
    X = np.c_[sin, constant]
    pred = gbdt.predict(X)
    assert np.all((np.diff(pred) >= 0) == (np.diff(sin) >= 0))

    # Second feature
    # assert pred is all decreasing when f_1 is all decreasing
    X = np.c_[constant, linspace]
    pred = gbdt.predict(X)
    assert is_decreasing(pred)
    # assert pred actually follows the inverse variations of f_1
    X = np.c_[constant, sin]
    pred = gbdt.predict(X)
    assert ((np.diff(pred) <= 0) == (np.diff(sin) >= 0)).all()
