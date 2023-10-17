# License: BSD 3 clause

import numpy as np

from sklearn.callback import ComputationTree

levels = [
    {"descr": "level0", "max_iter": 3},
    {"descr": "level1", "max_iter": 5},
    {"descr": "level2", "max_iter": 7},
    {"descr": "level3", "max_iter": None},
]


def test_computation_tree():
    """Check the construction of the computation tree"""
    computation_tree = ComputationTree(estimator_name="estimator", levels=levels)
    assert computation_tree.estimator_name == "estimator"

    root = computation_tree.root
    assert root.parent is None
    assert root.idx == 0

    assert len(root.children) == root.max_iter == 3
    assert [node.idx for node in root.children] == list(range(3))

    for node1 in root.children:
        assert len(node1.children) == 5
        assert [n.idx for n in node1.children] == list(range(5))

        for node2 in node1.children:
            assert len(node2.children) == 7
            assert [n.idx for n in node2.children] == list(range(7))

            for node3 in node2.children:
                assert not node3.children


def test_n_nodes():
    """Check that the number of node in a comutation tree corresponds to what we expect
    from the level descriptions
    """
    computation_tree = ComputationTree(estimator_name="", levels=levels)

    max_iter_per_level = [level["max_iter"] for level in levels[:-1]]
    expected_n_nodes = 1 + np.sum(np.cumprod(max_iter_per_level))

    assert computation_tree.n_nodes == expected_n_nodes
    assert len(computation_tree.iterate(include_leaves=True)) == expected_n_nodes
    assert computation_tree._tree_status.shape == (expected_n_nodes,)


def test_tree_status_idx():
    """Check that each node has a unique index in the _tree_status array and that their
    order corresponds to the order given by a depth first search.
    """
    computation_tree = ComputationTree(estimator_name="", levels=levels)

    indexes = [
        node.tree_status_idx for node in computation_tree.iterate(include_leaves=True)
    ]
    assert indexes == list(range(computation_tree.n_nodes))


def test_get_ancestors():
    """Check the ancestor search and its propagation to parent trees"""
    parent_levels = [
        {"descr": "parent_level0", "max_iter": 2},
        {"descr": "parent_level1", "max_iter": 4},
        {"descr": "parent_level2", "max_iter": None},
    ]

    parent_computation_tree = ComputationTree(
        estimator_name="parent_estimator", levels=parent_levels
    )
    parent_node = parent_computation_tree.root.children[0].children[2]
    # indices of each node (in its parent children) in this chain are 0, 0, 2.
    # (root is always 0).
    expected_parent_indices = [2, 0, 0]

    computation_tree = ComputationTree(
        estimator_name="estimator", levels=levels, parent_node=parent_node
    )
    node = computation_tree.root.children[1].children[3].children[5]
    expected_node_indices = [5, 3, 1, 0]

    ancestors = node.get_ancestors(include_ancestor_trees=False)
    assert ancestors == [
        node,
        node.parent,
        node.parent.parent,
        node.parent.parent.parent,
    ]
    assert [n.idx for n in ancestors] == expected_node_indices
    assert computation_tree.root in ancestors

    ancestors = node.get_ancestors(include_ancestor_trees=True)
    assert ancestors == [
        node,
        node.parent,
        node.parent.parent,
        node.parent.parent.parent,
        parent_node,
        parent_node.parent,
        parent_node.parent.parent,
    ]
    assert [n.idx for n in ancestors] == expected_node_indices + expected_parent_indices
