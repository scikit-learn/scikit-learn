# License: BSD 3 clause
# Authors: the scikit-learn developers

import numpy as np

from sklearn.callback import build_computation_tree

LEVELS = [
    {"descr": "level0", "max_iter": 3},
    {"descr": "level1", "max_iter": 5},
    {"descr": "level2", "max_iter": 7},
    {"descr": "level3", "max_iter": None},
]


def test_computation_tree():
    """Check the construction of the computation tree."""
    computation_tree = build_computation_tree(estimator_name="estimator", levels=LEVELS)
    assert computation_tree.estimator_name == ("estimator",)
    assert computation_tree.parent is None
    assert computation_tree.idx == 0

    assert len(computation_tree.children) == computation_tree.max_iter == 3
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
    computation_tree = build_computation_tree(estimator_name="", levels=LEVELS)

    max_iter_per_level = [level["max_iter"] for level in LEVELS[:-1]]
    expected_n_nodes = 1 + np.sum(np.cumprod(max_iter_per_level))

    actual_n_nodes = sum(1 for _ in computation_tree)

    assert actual_n_nodes == expected_n_nodes


def test_path():
    """Check that the path from the root to a node is correct."""
    computation_tree = build_computation_tree(estimator_name="", levels=LEVELS)

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
