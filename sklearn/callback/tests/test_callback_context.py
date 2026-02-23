# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pytest

from sklearn.callback._callback_context import CallbackContext, get_context_path
from sklearn.callback.tests._utils import (
    MaxIterEstimator,
    MetaEstimator,
    NoCallbackEstimator,
    NoSubtaskEstimator,
    ParentFitEstimator,
    TestingAutoPropagatedCallback,
    TestingCallback,
    ThirdPartyEstimator,
)


def test_propagate_callback_context():
    """Sanity check for the `propagate_callback_context` method."""
    not_propagated_callback = TestingCallback()
    propagated_callback = TestingAutoPropagatedCallback()

    estimator = MaxIterEstimator()
    metaestimator = MetaEstimator(estimator)
    metaestimator.set_callbacks([not_propagated_callback, propagated_callback])

    callback_ctx = CallbackContext._from_estimator(metaestimator, "fit", 0, 0)
    callback_ctx.propagate_callback_context(estimator)

    assert hasattr(estimator, "_parent_callback_ctx")
    assert not_propagated_callback not in estimator._skl_callbacks
    assert propagated_callback in estimator._skl_callbacks


def test_propagate_callback_context_no_callback():
    """Check that no callback is propagated if there's no callback."""
    estimator = MaxIterEstimator()
    metaestimator = MetaEstimator(estimator)

    callback_ctx = CallbackContext._from_estimator(metaestimator, "fit", 0, 0)
    assert len(callback_ctx._callbacks) == 0

    callback_ctx.propagate_callback_context(estimator)

    assert not hasattr(metaestimator, "_skl_callbacks")
    assert not hasattr(estimator, "_skl_callbacks")


def test_auto_propagated_callbacks():
    """Check that it's not possible to set an auto-propagated callback on the
    sub-estimator of a meta-estimator.
    """
    estimator = MaxIterEstimator()
    estimator.set_callbacks(TestingAutoPropagatedCallback())
    meta_estimator = MetaEstimator(estimator=estimator)

    match = (
        r"sub-estimator .*of a meta-estimator .*can't have auto-propagated callbacks"
    )
    with pytest.raises(TypeError, match=match):
        meta_estimator.fit(X=None, y=None)


def _make_task_tree(n_children, n_grandchildren):
    """Helper function to create a tree of tasks with a context for each task."""
    estimator = MaxIterEstimator()
    root = CallbackContext._from_estimator(estimator, "root task", 0, n_children)

    for i in range(n_children):
        child = root.subcontext(
            task_name="child task",
            task_id=i,
            max_subtasks=n_grandchildren,
        )

        for j in range(n_grandchildren):
            grandchild = child.subcontext(
                task_name="grandchild task",
                task_id=j,
                max_subtasks=0,
            )

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

        for grandchild in child._children_map.values():
            assert grandchild.parent is child
            assert len(get_context_path(grandchild)) == 3
            assert len(grandchild._children_map) == 0

    # all _children_map lengths should match the max_subtasks attribute.
    for node in root:
        assert len(node._children_map) == node.max_subtasks

    # 1 root + 1 * 3 children + 1 * 3 * 5 grandchildren
    expected_n_nodes = np.sum(np.cumprod([1, 3, 5]))
    actual_n_nodes = sum(1 for _ in root)
    assert actual_n_nodes == expected_n_nodes

    # None of the nodes should have been merged with another node
    assert all(node.source_estimator_name is None for node in root)
    assert all(node.source_task_name is None for node in root)


def test_add_child():
    """Sanity check for the `_add_child` method."""
    estimator = MaxIterEstimator()
    root = CallbackContext._from_estimator(estimator, "root task", 0, max_subtasks=2)

    first_child = CallbackContext._from_estimator(
        estimator, "child task", task_id=0, max_subtasks=0
    )

    root._add_child(first_child)
    assert root.max_subtasks == 2
    assert len(root._children_map) == 1
    assert first_child.task_id == 0

    # root already has a child with id 0
    second_child = CallbackContext._from_estimator(
        estimator, "child task", task_id=0, max_subtasks=0
    )
    with pytest.raises(
        ValueError, match=r"Callback context .* already has a child with task_id=0"
    ):
        root._add_child(second_child)

    second_child.task_id = 1
    root._add_child(second_child)
    assert len(root._children_map) == 2

    # root can have at most 2 children
    third_child = CallbackContext._from_estimator(
        estimator, "child task", task_id=2, max_subtasks=0
    )

    with pytest.raises(ValueError, match=r"Cannot add child to callback context"):
        root._add_child(third_child)


def test_merge_with():
    """Sanity check for the `_merge_with` method."""
    estimator = MaxIterEstimator()
    meta_estimator = MetaEstimator(estimator)
    outer_root = CallbackContext._from_estimator(
        meta_estimator, "root", 0, max_subtasks=None
    )

    # Add a child task within the same estimator
    outer_child = outer_root.subcontext(task_name="child", task_id=0, max_subtasks=0)

    # The root task of the inner estimator is merged with (and effectively replaces)
    # a leaf of the outer estimator because they correspond to the same formal task.
    inner_root = CallbackContext._from_estimator(estimator, "root", 0, None)
    inner_root._merge_with(outer_child)

    assert inner_root.parent is outer_root
    assert inner_root.task_id == outer_child.task_id
    assert outer_child not in outer_root._children_map.values()
    assert inner_root in outer_root._children_map.values()

    # The name and estimator name of the tasks it was merged with are stored
    assert inner_root.source_task_name == outer_child.task_name
    assert inner_root.source_estimator_name == outer_child.estimator_name


def test_merge_with_error_not_leaf():
    """Check that merging with a non-leaf node raises an error."""
    estimator = MaxIterEstimator()
    inner_root = CallbackContext._from_estimator(estimator, "root", 0, None)

    meta_estimator = MetaEstimator(estimator)
    outer_root = CallbackContext._from_estimator(meta_estimator, "root", 0, None)

    # Add a child task within the same estimator
    outer_root.subcontext(task_name="child", task_id=0, max_subtasks=1)

    with pytest.raises(ValueError, match=r"Cannot merge callback context"):
        inner_root._merge_with(outer_root)


@pytest.mark.parametrize(
    "estimator_class", [MaxIterEstimator, ThirdPartyEstimator, ParentFitEstimator]
)
def test_callback_ctx_removed_after_fit(estimator_class):
    """Check that the _callback_fit_ctx attribute gets removed after fit."""
    estimator = estimator_class().fit()
    assert not hasattr(estimator, "_callback_fit_ctx")


def test_inner_estimator_no_callback_support():
    """Check that meta estimators can have sub estimators without callback support.

    No error is raised when the sub-estimator does not support callbacks. If callbacks
    would be propagated, a warning is raised instead.
    """
    estimator = NoCallbackEstimator()
    meta_estimator = MetaEstimator(estimator)
    meta_estimator.set_callbacks(TestingAutoPropagatedCallback())

    with pytest.warns(
        UserWarning,
        match="The estimator NoCallbackEstimator does not support callbacks.",
    ):
        meta_estimator.fit()


def test_estimator_without_subtask():
    """Check that callback support works for an estimator without subtasks.

    This test is about verifying that an estimator that does not call its callback
    context's `eval_on_fit_task_end` does not cause a problem.
    """
    estimator = NoSubtaskEstimator()
    estimator.set_callbacks([TestingCallback()])
    estimator.fit()


def test_propagate_callback_context_forwards_callbacks_when_cloning():
    """Test propagate_callback_context for callback preservation when cloning."""
    estimator = MaxIterEstimator()
    callback_ctx = CallbackContext._from_estimator(
        estimator, task_name="", task_id=0, max_subtasks=0
    )
    callback = TestingCallback()
    estimator.set_callbacks(callback)
    cloned_estimator = callback_ctx.propagate_callback_context(
        estimator, clone_estimator=True
    )
    assert cloned_estimator is not estimator
    assert cloned_estimator._skl_callbacks == estimator._skl_callbacks == [callback]
