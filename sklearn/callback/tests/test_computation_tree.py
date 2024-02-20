# License: BSD 3 clause
# Authors: the scikit-learn developers

import numpy as np

from sklearn.callback import build_computation_tree

TREE_STRUCTURE = [
    {"stage": "stage0", "n_children": 3},
    {"stage": "stage1", "n_children": 5},
    {"stage": "stage2", "n_children": 7},
    {"stage": "stage3", "n_children": None},
]


def test_computation_tree():
    """Check the construction of the computation tree."""
    computation_tree = build_computation_tree(
        estimator_name="estimator", tree_structure=TREE_STRUCTURE
    )
    assert computation_tree.estimator_name == ("estimator",)
    assert computation_tree.parent is None
    assert computation_tree.idx == 0

    assert len(computation_tree.children) == computation_tree.n_children == 3
    assert [node.idx for node in computation_tree.children] == list(range(3))

    for node1 in computation_tree.children:
        assert len(node1.children) == 5
        assert [n.idx for n in node1.children] == list(range(5))

        for node2 in node1.children:
            assert len(node2.children) == 7
            assert [n.idx for n in node2.children] == list(range(7))

            for node3 in node2.children:
                assert not node3.children


def test_n_nodes():
    """Check that the number of node in a computation tree corresponds to what we expect
    from the level descriptions.
    """
    computation_tree = build_computation_tree(
        estimator_name="", tree_structure=TREE_STRUCTURE
    )

    n_children_per_level = [stage["n_children"] for stage in TREE_STRUCTURE[:-1]]
    expected_n_nodes = 1 + np.sum(np.cumprod(n_children_per_level))

    actual_n_nodes = sum(1 for _ in computation_tree)

    assert actual_n_nodes == expected_n_nodes


def test_path():
    """Check that the path from the root to a node is correct."""
    computation_tree = build_computation_tree(
        estimator_name="", tree_structure=TREE_STRUCTURE
    )

    assert computation_tree.path == [computation_tree]

    node = computation_tree.children[1].children[2].children[3]
    expected_path = [
        computation_tree,
        computation_tree.children[1],
        computation_tree.children[1].children[2],
        node,
    ]
    assert node.path == expected_path

    assert all(node.depth == i for i, node in enumerate(expected_path))
