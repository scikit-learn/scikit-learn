# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.callback._callback_context import get_context_path
from sklearn.callback._callback_support import get_callback_manager


class SaveContextTree:
    """Callback to save the context tree of the estimator.

    Attributes
    ----------
    context_tree_ : CallbackContext
        The root context of the saved context tree.
    """

    max_propagation_depth = None

    def setup(self, estimator, context):
        self._leaf_contexts = get_callback_manager().list()

    def on_fit_task_begin(self, estimator, context):
        pass

    def on_fit_task_end(self, estimator, context):
        if context.max_subtasks == 0:
            self._leaf_contexts.append(context)

    def teardown(self, estimator, context):
        leaves = list(self._leaf_contexts)
        paths = [get_context_path(c) for c in leaves]

        # Sort paths by task_ids to retrieve sequential tasks ordering.
        paths = sorted(paths, key=lambda p: [ctx.task_id for ctx in p])
        base_path = paths[0]

        for path in paths[1:]:
            for ctx, base_ctx in zip(path, base_path):
                if ctx.task_id != base_ctx.task_id:
                    ctx.parent = base_ctx.parent
                    base_ctx.parent._children_map[ctx.task_id] = ctx
                    break

        self.context_tree_ = base_path[0]
