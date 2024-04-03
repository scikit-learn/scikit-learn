# License: BSD 3 clause
# Authors: the scikit-learn developers

import numpy as np

from sklearn.callback import TaskNode


def _make_task_tree(n_children, n_grandchildren):
    root = TaskNode(
        estimator_name="estimator", name="root task", max_subtasks=n_children
    )

    for i in range(n_children):
        child = root._add_child(name="child task", max_subtasks=n_grandchildren, idx=i)

        for j in range(n_grandchildren):
            child._add_child(name="grandchild task", max_subtasks=0, idx=j)

    return root


def test_task_tree():
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.max_subtasks == 3
    assert root.idx is None
    assert root.parent is None

    assert len(root.children) == 3
    assert all(len(child.children) == 5 for child in root.children.values())

    # 1 root, 3 children, 3 * 5 grandchildren
    expected_n_nodes = np.sum(np.cumprod([1, 3, 5]))
    actual_n_nodes = sum(1 for _ in root)
    assert actual_n_nodes == expected_n_nodes


def test_path():
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.path == [root]

    # pick a node
    node = root.children[1].children[2]
    expected_path = [root, root.children[1], node]
    assert node.path == expected_path


def test_depth():
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.depth == 0

    assert all(child.depth == 1 for child in root.children.values())
    assert all(
        grandchild.depth == 2
        for child in root.children.values()
        for grandchild in child.children.values()
    )
