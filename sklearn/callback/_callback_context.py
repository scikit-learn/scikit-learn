# License: BSD 3 clause
# Authors: the scikit-learn developers

from ._task_tree import TaskNode


class CallbackContext:
    def __init__(self, callbacks):
        self.callbacks = callbacks

    def eval_on_fit_begin(self, *, estimator, max_subtasks=0, data):
        self.task_node = TaskNode(
            estimator_name=estimator.__class__.__name__,
            name="fit",
            max_subtasks=max_subtasks,
        )

        parent_task_node = getattr(estimator, "_parent_task_node", None)
        if parent_task_node is not None:
            self.task_node._merge_with(parent_task_node)

        for callback in self.callbacks:
            # Only call the on_fit_begin method of callbacks that are not
            # propagated from a meta-estimator.
            if not (callback.auto_propagate and parent_task_node is not None):
                callback.on_fit_begin(estimator, data=data)

        return self

    def eval_on_fit_iter_end(self, **kwargs):
        return any(
            callback.on_fit_iter_end(task_node=self.task_node, **kwargs)
            for callback in self.callbacks
        )

    def eval_on_fit_end(self):
        for callback in self.callbacks:
            # Only call the on_fit_end method of callbacks that are not
            # propagated from a meta-estimator.
            if not (callback.auto_propagate and self.task_node.parent is not None):
                callback.on_fit_end(task_node=self.task_node)

    def subcontext(self, task="", max_subtasks=0, idx=0, sub_estimator=None):
        sub_ctx = CallbackContext(callbacks=self.callbacks)

        sub_ctx.task_node = self.task_node._add_child(
            name=task,
            max_subtasks=max_subtasks,
            idx=idx,
        )

        if sub_estimator is not None:
            sub_ctx._propagate_callbacks(sub_estimator=sub_estimator)

        return sub_ctx

    def _propagate_callbacks(self, sub_estimator):
        bad_callbacks = [
            callback.__class__.__name__
            for callback in getattr(sub_estimator, "_skl_callbacks", [])
            if callback.auto_propagate
        ]

        if bad_callbacks:
            raise TypeError(
                f"The sub-estimator ({sub_estimator.__class__.__name__}) of a"
                f" meta-estimator ({self.task_node.estimator_name}) can't have"
                f" auto-propagated callbacks ({bad_callbacks})."
                " Register them directly on the meta-estimator."
            )

        callbacks_to_propagate = [
            callback for callback in self.callbacks if callback.auto_propagate
        ]

        if not callbacks_to_propagate:
            return

        sub_estimator._parent_task_node = self.task_node

        sub_estimator._set_callbacks(
            getattr(sub_estimator, "_skl_callbacks", []) + callbacks_to_propagate
        )
