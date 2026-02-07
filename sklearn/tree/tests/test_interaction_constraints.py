# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest
from scipy.sparse import csc_matrix

from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.utils._testing import assert_allclose


def _assert_tree_respects_interaction_constraints(tree_, interaction_cst):
    constraints = [set(group) for group in interaction_cst]

    stack = [(0, set())]
    while stack:
        node_id, path_features = stack.pop()
        feature_idx = tree_.feature[node_id]
        if feature_idx < 0:
            continue

        new_path_features = path_features | {feature_idx}
        assert any(new_path_features.issubset(group) for group in constraints)

        stack.append((tree_.children_left[node_id], new_path_features))
        stack.append((tree_.children_right[node_id], new_path_features))


@pytest.mark.parametrize(
    "interaction_cst, n_features, result",
    [
        (None, 12, None),
        ([{0, 1}], 2, [{0, 1}]),
        ("pairwise", 2, [{0, 1}]),
        ("pairwise", 4, [{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}]),
        ("no_interactions", 2, [{0}, {1}]),
        ("no_interactions", 4, [{0}, {1}, {2}, {3}]),
        ([(1, 0), [5, 1]], 6, [{0, 1}, {1, 5}, {2, 3, 4}]),
    ],
)
def test_check_interaction_cst(interaction_cst, n_features, result):
    est = DecisionTreeRegressor()
    est.set_params(interaction_cst=interaction_cst)
    assert est._check_interaction_cst(n_features) == result


@pytest.mark.parametrize(
    "interaction_cst, expected_error",
    [
        ([0, 1], "Interaction constraints must be a sequence"),
        ([{0, 9999}], "Interaction constraints must consist of integer indices"),
        ([{-1, 0}], "Interaction constraints must consist of integer indices"),
        ([{0.5}], "Interaction constraints must consist of integer indices"),
    ],
)
def test_invalid_interaction_cst_raises(interaction_cst, expected_error):
    rng = np.random.RandomState(0)
    X = rng.randn(40, 2)
    y = rng.randn(40)

    est = DecisionTreeRegressor(interaction_cst=interaction_cst, random_state=0)
    with pytest.raises(ValueError, match=expected_error):
        est.fit(X, y)


def test_interaction_cst_random_splitter_raises():
    rng = np.random.RandomState(0)
    X = rng.randn(40, 2)
    y = rng.randn(40)

    est = DecisionTreeRegressor(
        splitter="random", interaction_cst=[{0, 1}], random_state=0
    )
    msg = "Interaction constraints are only supported for splitter='best'."
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)


def test_interaction_cst_best_first_raises():
    rng = np.random.RandomState(0)
    X = rng.randn(40, 2)
    y = rng.randn(40)

    est = DecisionTreeRegressor(
        max_leaf_nodes=8, interaction_cst=[{0, 1}], random_state=0
    )
    msg = "Interaction constraints are only supported for depth-first tree building."
    with pytest.raises(ValueError, match=msg):
        est.fit(X, y)


@pytest.mark.parametrize("Tree", [DecisionTreeRegressor, DecisionTreeClassifier])
@pytest.mark.parametrize("sparse_input", [False, True])
def test_interaction_cst_tree_structure(Tree, sparse_input):
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(400, 6))
    y = (
        X[:, 0] * X[:, 1]
        + X[:, 1] * X[:, 2]
        + X[:, 3] * X[:, 4]
        + 0.1 * X[:, 5]
        + 0.01 * rng.randn(X.shape[0])
    )
    if Tree is DecisionTreeClassifier:
        y = y > np.median(y)

    interaction_cst = [{0, 1}, {1, 2}, {3, 4, 5}]
    est = Tree(max_depth=5, random_state=0, interaction_cst=interaction_cst)
    X_fit = csc_matrix(X) if sparse_input else X
    est.fit(X_fit, y)

    assert est.tree_.node_count > 1
    _assert_tree_respects_interaction_constraints(
        est.tree_, est._check_interaction_cst(X.shape[1])
    )


def test_interaction_cst_overlapping_groups_tree_structure():
    rng = np.random.RandomState(0)
    X = rng.uniform(size=(60, 4))
    y = (
        10 * (X[:, 2] > 0.5)
        + 5 * (X[:, 0] > 0.6)
        + 3 * (X[:, 1] > 0.7)
        + 2 * (X[:, 3] > 0.4)
        + 4 * X[:, 0] * X[:, 2]
        + 2 * X[:, 1] * X[:, 3]
    )

    interaction_cst = [{0, 1, 3}, {1, 2, 3}, {0, 2}]
    est = DecisionTreeRegressor(
        max_depth=4, max_features=4, random_state=0, interaction_cst=interaction_cst
    )
    est.fit(X, y)

    _assert_tree_respects_interaction_constraints(
        est.tree_, est._check_interaction_cst(X.shape[1])
    )


def test_no_interactions_numerically():
    rng = np.random.RandomState(42)
    X = rng.uniform(size=(1000, 2))
    # Build a target with a strong interaction term between feature 0 and 1:
    # y = x0 + x1 + 5 * x0 * x1.
    y = X[:, 0] + X[:, 1] + 5 * X[:, 0] * X[:, 1]

    est = DecisionTreeRegressor(
        max_depth=8,
        interaction_cst="no_interactions",
        random_state=0,
    )
    est.fit(X, y)

    delta = 0.25
    # Keep evaluation points inside the training domain to avoid extrapolation.
    X_test = X[(X[:, 0] < 1 - delta) & (X[:, 1] < 1 - delta)]
    X_delta_d_0 = X_test + [delta, 0]
    X_delta_0_d = X_test + [0, delta]
    X_delta_d_d = X_test + [delta, delta]

    # For unconstrained models, this finite-difference expression captures the
    # second-order mixed term and should be positive with the interaction above.
    # With "no_interactions", features cannot interact along a tree path, so the
    # mixed term must vanish.
    second_order_effect = (
        est.predict(X_delta_d_d)
        + est.predict(X_test)
        - est.predict(X_delta_d_0)
        - est.predict(X_delta_0_d)
    )
    assert_allclose(second_order_effect, 0, atol=1e-12)
