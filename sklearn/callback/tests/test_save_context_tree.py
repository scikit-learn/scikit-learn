# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn.callback import SaveContextTree
from sklearn.callback._callback_context import CallbackContext
from sklearn.callback.tests._utils import MaxIterEstimator, MetaEstimator


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("prefer", ["threads", "processes"])
def test_save_context_tree_single_estimator(n_jobs, prefer):
    """Check that the context tree is correctly saved / reconstructed."""
    n_outer, n_inner, max_iter = 2, 3, 5

    cb = SaveContextTree()
    est = MaxIterEstimator(max_iter=max_iter)
    MetaEstimator(
        est, n_outer=n_outer, n_inner=n_inner, n_jobs=n_jobs, prefer=prefer
    ).set_callbacks(cb).fit()

    assert isinstance(cb.context_tree_, CallbackContext)

    assert cb.context_tree_.task_name == "fit"
    assert cb.context_tree_.max_subtasks == n_outer
    assert len(cb.context_tree_._children_map) == n_outer

    for i, child in enumerate(cb.context_tree_._children_map.values()):
        assert child.task_name == "outer"
        assert child.max_subtasks == n_inner
        assert len(child._children_map) == n_inner

        for j, grandchild in enumerate(child._children_map.values()):
            # Here we're in the subestimator's task tree.
            assert grandchild.task_name == "fit"
            assert grandchild.source_task_name == "inner"
            assert grandchild.task_id == j
            assert grandchild.max_subtasks == max_iter
            assert len(grandchild._children_map) == max_iter

            for k, leaf in enumerate(grandchild._children_map.values()):
                assert leaf.task_name == f"iteration {k}"
                assert leaf.task_id == k
                assert leaf.max_subtasks == 0
                assert len(leaf._children_map) == 0
