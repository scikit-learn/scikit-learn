# License: BSD 3 clause
# Authors: the scikit-learn developers

import numpy as np

from sklearn.callback import TaskNode


def _make_task_tree(n_children, n_grandchildren):
    root = TaskNode(
        task_name="root task", task_id=0, max_tasks=1, estimator_name="estimator"
    )

    for i in range(n_children):
        child = TaskNode(
            task_name="child task",
            task_id=i,
            max_tasks=n_children,
            estimator_name="estimator",
        )
        root._add_child(child)

        for j in range(n_grandchildren):
            grandchild = TaskNode(
                task_name="grandchild task",
                task_id=j,
                max_tasks=n_grandchildren,
                estimator_name="estimator",
            )
            child._add_child(grandchild)

    return root


def test_task_tree():
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.parent is None
    assert root.depth == 0
    assert len(root.children_map) == 3

    for child in root.children_map.values():
        assert child.parent is root
        assert child.depth == 1
        assert len(child.children_map) == 5
        assert root.max_subtasks == child.max_tasks

        for grandchild in child.children_map.values():
            assert grandchild.parent is child
            assert grandchild.depth == 2
            assert len(grandchild.children_map) == 0
            assert child.max_subtasks == grandchild.max_tasks

    # 1 root + 1 * 3 children + 1 * 3 * 5 grandchildren
    expected_n_nodes = np.sum(np.cumprod([1, 3, 5]))
    actual_n_nodes = sum(1 for _ in root)
    assert actual_n_nodes == expected_n_nodes

    # None of the nodes should have been merged with another node
    assert all(node.prev_estimator_name is None for node in root)
    assert all(node.prev_task_name is None for node in root)


def test_path():
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.path == [root]

    # pick an arbitrary node
    node = root.children_map[1].children_map[2]

    expected_path = [root, root.children_map[1], node]
    assert node.path == expected_path
