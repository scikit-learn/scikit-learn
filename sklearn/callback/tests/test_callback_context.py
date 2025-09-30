# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest

from sklearn.callback._callback_context import CallbackContext, get_context_path
from sklearn.callback.tests._utils import (
    Estimator,
    MetaEstimator,
    TestingAutoPropagatedCallback,
    TestingCallback,
)


def test_propagate_callbacks():
    """Sanity check for the `propagate_callbacks` method."""
    not_propagated_callback = TestingCallback()
    propagated_callback = TestingAutoPropagatedCallback()

    estimator = Estimator()
    metaestimator = MetaEstimator(estimator)
    metaestimator.set_callbacks([not_propagated_callback, propagated_callback])

    callback_ctx = metaestimator.__skl_init_callback_context__()
    callback_ctx.propagate_callbacks(estimator)

    assert hasattr(estimator, "_parent_callback_ctx")
    assert not_propagated_callback not in estimator._skl_callbacks
    assert propagated_callback in estimator._skl_callbacks


def test_propagate_callback_no_callback():
    """Check that no callback is propagated if there's no callback."""
    estimator = Estimator()
    metaestimator = MetaEstimator(estimator)

    callback_ctx = metaestimator.__skl_init_callback_context__()
    assert len(callback_ctx._callbacks) == 0

    callback_ctx.propagate_callbacks(estimator)

    assert not hasattr(metaestimator, "_skl_callbacks")
    assert not hasattr(estimator, "_skl_callbacks")


def test_auto_propagated_callbacks():
    """Check that it's not possible to set an auto-propagated callback on the
    sub-estimator of a meta-estimator.
    """
    estimator = Estimator()
    estimator.set_callbacks(TestingAutoPropagatedCallback())
    meta_estimator = MetaEstimator(estimator=estimator)

    match = (
        r"sub-estimator .*of a meta-estimator .*can't have auto-propagated callbacks"
    )
    with pytest.raises(TypeError, match=match):
        meta_estimator.fit(X=None, y=None)


def _make_task_tree(n_children, n_grandchildren):
    """Helper function to create a tree of tasks with a context for each task."""
    estimator = Estimator()
    root = CallbackContext._from_estimator(
        estimator,
        task_name="root task",
        task_id=0,
        max_subtasks=n_children,
    )

    for i in range(n_children):
        child = CallbackContext._from_estimator(
            estimator,
            task_name="child task",
            task_id=i,
            max_subtasks=n_grandchildren,
        )
        root._add_child(child)

        for j in range(n_grandchildren):
            grandchild = CallbackContext._from_estimator(
                estimator,
                task_name="grandchild task",
                task_id=j,
            )
            child._add_child(grandchild)

    return root


def test_task_tree():
    """Check that the task tree is correctly built."""
    root = _make_task_tree(n_children=3, n_grandchildren=5)

    assert root.parent is None
    assert len(get_context_path(root)) == 1
    assert len(root._children_map) == 3

    for child in root._children_map.values():
        assert child.parent is root
        assert len(get_context_path(child)) == 2
        assert len(child._children_map) == 5
        assert root.max_subtasks == 3

        for grandchild in child._children_map.values():
            assert grandchild.parent is child
            assert len(get_context_path(grandchild)) == 3
            assert len(grandchild._children_map) == 0
            assert child.max_subtasks == 5

    # 1 root + 1 * 3 children + 1 * 3 * 5 grandchildren
    expected_n_nodes = np.sum(np.cumprod([1, 3, 5]))
    actual_n_nodes = sum(1 for _ in root)
    assert actual_n_nodes == expected_n_nodes

    # None of the nodes should have been merged with another node
    assert all(node.prev_estimator_name is None for node in root)
    assert all(node.prev_task_name is None for node in root)


def test_add_child():
    """Sanity check for the `_add_child` method."""
    estimator = Estimator()
    root = CallbackContext._from_estimator(
        estimator, task_name="root task", task_id=0, max_subtasks=2
    )

    root._add_child(
        CallbackContext._from_estimator(estimator, task_name="child task", task_id=0)
    )
    assert root.max_subtasks == 2
    assert len(root._children_map) == 1

    # root already has a child with id 0
    with pytest.raises(
        ValueError, match=r"Callback context .* already has a child with task_id=0"
    ):
        root._add_child(
            CallbackContext._from_estimator(
                estimator, task_name="child task", task_id=0
            )
        )

    root._add_child(
        CallbackContext._from_estimator(estimator, task_name="child task", task_id=1)
    )
    assert len(root._children_map) == 2

    # root can have at most 2 children
    with pytest.raises(ValueError, match=r"Cannot add child to callback context"):
        root._add_child(
            CallbackContext._from_estimator(
                estimator, task_name="child task", task_id=2
            )
        )


def test_merge_with():
    """Sanity check for the `_merge_with` method."""
    estimator = Estimator()
    meta_estimator = MetaEstimator(estimator)
    outer_root = CallbackContext._from_estimator(
        meta_estimator, task_name="root", task_id=0, max_subtasks=2
    )

    # Add a child task within the same estimator
    outer_child = CallbackContext._from_estimator(
        meta_estimator, task_name="child", task_id="id", max_subtasks=1
    )
    outer_root._add_child(outer_child)

    # The root task of the inner estimator is merged with (and effectively replaces)
    # a leaf of the outer estimator because they correspond to the same formal task.
    inner_root = CallbackContext._from_estimator(estimator, task_name="root", task_id=0)
    inner_root._merge_with(outer_child)

    assert inner_root.parent is outer_root
    assert inner_root.task_id == outer_child.task_id
    assert outer_child not in outer_root._children_map.values()
    assert inner_root in outer_root._children_map.values()

    # The name and estimator name of the tasks it was merged with are stored
    assert inner_root.prev_task_name == outer_child.task_name
    assert inner_root.prev_estimator_name == outer_child.estimator_name
