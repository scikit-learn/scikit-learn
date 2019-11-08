import numpy as np
import pytest
from pytest import approx

from sklearn.ensemble._hist_gradient_boosting.grower import TreeGrower
from sklearn.ensemble._hist_gradient_boosting.common import G_H_DTYPE


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


    if monotonic_cst == 0:  # NO_CST
        # some increasing, some decreasing
        assert (any(v1 < v2 for (v1, v2) in zip(values, values[1:])) and
                any(v1 > v2 for (v1, v2) in zip(values, values[1:])))
    elif monotonic_cst == 1:  # INC
        # all increasing
        assert all(v1 < v2 for (v1, v2) in zip(values, values[1:]))
    else:  # DEC
        # all decreasing
        assert all(v1 > v2 for (v1, v2) in zip(values, values[1:]))


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
        else:
            left_greater.append(node)
        dfs(left_idx)
        dfs(right_idx)
    dfs(0)  # start at root (0)

    if monotonic_cst == 0:  # NO_CST
        assert left_lower and left_greater
    elif monotonic_cst == 1:  # INC
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

    if monotonic_cst == 0:  # NO_CST
        return

    def dfs(node):
        if node.is_leaf:
            return
        if node is not grower.root and node is node.parent.left_child:
            sibling = node.sibling  # on the right
            middle = (node.value + sibling.value) / 2
            if monotonic_cst == 1:  # INC
                assert node.left_child.value < node.right_child.value < middle
                if not sibling.is_leaf:
                    assert middle < sibling.left_child.value < sibling.right_child.value
            else:  # DEC
                assert node.left_child.value > node.right_child.value > middle
                if not sibling.is_leaf:
                    assert middle > sibling.left_child.value > sibling.right_child.value

        dfs(node.left_child)
        dfs(node.right_child)
    dfs(grower.root)


@pytest.mark.parametrize('seed', range(3))
@pytest.mark.parametrize('monotonic_cst', (
    0, # NO_CST
    1, # INC
    2, # DEC
))
def test_grower(monotonic_cst, seed):
    # Build a single tree with only one feature, and make sure the predictor
    # respects the monotonic constraints.

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
                        monotonic_cst=[monotonic_cst])
    grower.grow()

    predictor = grower.make_predictor()
    assert_children_values_monotonic(predictor, monotonic_cst)
    assert_children_values_bounded(grower, monotonic_cst)
    assert_leaves_values_monotonic(predictor, monotonic_cst)

def is_increasing(y):
    return (np.diff(y) >= 0.0).all()

def is_decreasing(y):
    return (np.diff(y) <= 0.0).all()

# @pytest.mark.parametrize('seed', range(1))
# def test_light(seed):

#     def is_correctly_constrained(learner):
#         n = 10
#         variable_x = np.linspace(0, 1, n).reshape((n, 1))
#         fixed_xs_values = np.linspace(0, 1, n)
#         for i in range(n):
#             fixed_x = fixed_xs_values[i] * np.ones((n, 1))
#             monotonically_increasing_x = np.column_stack((variable_x, fixed_x))
#             monotonically_increasing_y = learner.predict(monotonically_increasing_x)
#             print(monotonically_increasing_y)
#             # monotonically_decreasing_x = np.column_stack((fixed_x, variable_x))
#             # monotonically_decreasing_y = learner.predict(monotonically_decreasing_x)
#             if not (is_increasing(monotonically_increasing_y)): # and is_decreasing(monotonically_decreasing_y)):
#                 return False
#             break
#         return True
#     np.random.seed(seed)

#     number_of_dpoints = 3000
#     x1_positively_correlated_with_y = np.random.random(size=number_of_dpoints)
#     x2_negatively_correlated_with_y = np.random.random(size=number_of_dpoints)
#     x = np.column_stack((x1_positively_correlated_with_y, x2_negatively_correlated_with_y))
#     zs = np.random.normal(loc=0.0, scale=0.01, size=number_of_dpoints)
#     y = (5 * x1_positively_correlated_with_y
#             + np.sin(10 * np.pi * x1_positively_correlated_with_y)
#             - 5 * x2_negatively_correlated_with_y
#             - np.cos(10 * np.pi * x2_negatively_correlated_with_y)
#             + zs)
#     from sklearn.experimental import enable_hist_gradient_boosting  # noqa
#     from sklearn.ensemble import HistGradientBoostingRegressor
#     # trainset = lgb.Dataset(x, label=y)
#     # params = {
#     #     'min_data': 20,
#     #     'num_leaves': 20,
#     #     'monotone_constraints': '1,-1'
#     # }
#     # constrained_model = lgb.train(params, trainset)
#     from sklearn.ensemble._hist_gradient_boosting.utils import (
#         get_equivalent_estimator)
#     gbdt = HistGradientBoostingRegressor(max_iter=42, max_leaf_nodes=20)
#     gbdt.monotonic_cst = [1, 0]
#     gbdt.fit(x, y)
#     l = get_equivalent_estimator(gbdt)
#     l.set_params(monotone_constraints='1,0')
#     l.fit(x, y)

#     assert is_correctly_constrained(l)
#     assert is_correctly_constrained(gbdt)


@pytest.mark.parametrize('seed', range(1))
def test_zob(seed):
    print()

    from sklearn.experimental import enable_hist_gradient_boosting  # noqa
    from sklearn.ensemble import HistGradientBoostingRegressor
    from sklearn.datasets import make_regression

    rng = np.random.RandomState(seed)

    n_samples = 1000
    n_features = 1
    X = rng.normal(size=(n_samples, n_features))
    y = X[:, 0]

    gbdt = HistGradientBoostingRegressor(max_iter=100, min_samples_leaf=20)
    gbdt.monotone_constraints = [1]
    # gbdt.monotone_constraints = [1, 0]
    gbdt.fit(X, y)

    # increasing X_0 (X_1 must be constant)
    n_samples = 10
    X = np.c_[np.linspace(-2, 2, n_samples)]
    # X = np.c_[np.linspace(-2, 2, n_samples), np.zeros(n_samples)]
    for i, predictor in enumerate(gbdt._predictors):
        predictor = predictor[0]
        pred = predictor.predict(X)
        if not is_increasing(pred):
            print(i)
            print(np.diff(pred))
            print(pred)
            assert 0

    # pred = gbdt.predict(X)
    # assert is_increasing(pred)
