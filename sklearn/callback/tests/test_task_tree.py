# License: BSD 3 clause
# Authors: the scikit-learn developers

import numpy as np
import pytest

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
    """Check that the task tree is correctly built."""
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
    """Sanity check for the path property."""
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.path == [root]

    # pick an arbitrary node
    node = root.children_map[1].children_map[2]

    expected_path = [root, root.children_map[1], node]
    assert node.path == expected_path


def test_add_task():
    """Check that informative error messages are raised when adding tasks."""
    root = TaskNode(task_name="root task", task_id=0, max_tasks=1, estimator_name="est")

    # Before adding new task, it's considered a leaf
    assert root.max_subtasks == 0

    root._add_child(
        TaskNode(task_name="child task", task_id=0, max_tasks=2, estimator_name="est")
    )
    assert root.max_subtasks == 2
    assert len(root.children_map) == 1

    # root already has a child with id 0
    with pytest.raises(
        ValueError, match=r"Task node .* already has a child with task_id=0"
    ):
        root._add_child(
            TaskNode(
                task_name="child task", task_id=0, max_tasks=2, estimator_name="est"
            )
        )

    root._add_child(
        TaskNode(task_name="child task", task_id=1, max_tasks=2, estimator_name="est")
    )
    assert len(root.children_map) == 2

    # root can have at most 2 children
    with pytest.raises(ValueError, match=r"Cannot add child to task node"):
        root._add_child(
            TaskNode(
                task_name="child task", task_id=2, max_tasks=2, estimator_name="est"
            )
        )


def test_merge_with():
    outer_root = TaskNode(
        task_name="root", task_id=0, max_tasks=1, estimator_name="outer"
    )

    # Add a child task within the same estimator
    outer_child = TaskNode(
        task_name="child", task_id="id", max_tasks=2, estimator_name="outer"
    )
    outer_root._add_child(outer_child)

    # The root task of the inner estimator is merged with (and effectively replaces)
    # a leaf of the outer estimator because they correspond to the same formal task.
    inner_root = TaskNode(
        task_name="root", task_id=0, max_tasks=1, estimator_name="inner"
    )
    inner_root._merge_with(outer_child)

    assert inner_root.parent is outer_root
    assert inner_root.task_id == outer_child.task_id
    assert outer_child not in outer_root.children_map.values()
    assert inner_root in outer_root.children_map.values()

    # The name and estimator name of the tasks it was merged with are stored
    assert inner_root.prev_task_name == outer_child.task_name
    assert inner_root.prev_estimator_name == outer_child.estimator_name
