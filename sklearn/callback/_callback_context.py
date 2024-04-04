# License: BSD 3 clause
# Authors: the scikit-learn developers

from ._task_tree import TaskNode


class CallbackContext:
    def __init__(
        self,
        callbacks,
        estimator_name="",
        task_name="",
        task_id=0,
        max_tasks=1,
        parent_task_node=None,
        parent_estimator_task_node=None,
    ):
        self.callbacks = callbacks
        self.estimator_name = estimator_name

        self.task_node = TaskNode(
            task_name=task_name,
            task_id=task_id,
            max_tasks=max_tasks,
            estimator_name=self.estimator_name,
        )

        if parent_task_node is not None:
            # This task is a subtask of another task of a same estimator
            parent_task_node._add_child(self.task_node)
        elif parent_estimator_task_node is not None:
            # This task is the root task of the estimator which itself corresponds to
            # a leaf task of a meta-estimator. Both tasks actually represent the same
            # task so we merge both task nodes into a single task node, attaching the
            # task tree of the sub-estimator to the task tree of the meta-estimator on
            # the way.
            self.task_node._merge_with(parent_estimator_task_node)

    def subcontext(self, task_name="", task_id=0, max_tasks=1):
        return CallbackContext(
            callbacks=self.callbacks,
            estimator_name=self.estimator_name,
            task_name=task_name,
            task_id=task_id,
            max_tasks=max_tasks,
            parent_task_node=self.task_node,
        )

    def eval_on_fit_begin(self, estimator, *, data):
        for callback in self.callbacks:
            # Only call the on_fit_begin method of callbacks that are not
            # propagated from a meta-estimator.
            if not (callback.auto_propagate and self.task_node.parent is not None):
                callback.on_fit_begin(estimator, data=data)

        return self

    def eval_on_fit_iter_end(self, estimator, **kwargs):
        return any(
            callback.on_fit_iter_end(estimator, self.task_node, **kwargs)
            for callback in self.callbacks
        )

    def eval_on_fit_end(self, estimator):
        for callback in self.callbacks:
            # Only call the on_fit_end method of callbacks that are not
            # propagated from a meta-estimator.
            if not (callback.auto_propagate and self.task_node.parent is not None):
                callback.on_fit_end(estimator, task_node=self.task_node)

    def propagate_callbacks(self, sub_estimator):
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

        sub_estimator._parent_estimator_task_node = self.task_node

        sub_estimator._set_callbacks(
            getattr(sub_estimator, "_skl_callbacks", []) + callbacks_to_propagate
        )

        return self
