# License: BSD 3 clause
# Authors: the scikit-learn developers

from . import AutoPropagatedProtocol
from ._task_tree import TaskNode


class CallbackContext:
    """Task level context for the callbacks.

    This class is responsible for managing the callbacks and task tree of an estimator.
    """

    @classmethod
    def _from_estimator(cls, estimator, *, task_name, task_id, max_tasks=1):
        """Private constructor to create a root context.

        Parameters
        ----------
        estimator : estimator instance
            The estimator this context is responsible for.

        task_name : str
            The name of the task this context is responsible for.

        task_id : int
            The id of the task this context is responsible for.

        max_tasks : int, default=1
            The maximum number of tasks that can be siblings of the task this context is
            responsible for.
        """
        new_ctx = cls.__new__(cls)

        # We don't store the estimator in the context to avoid circular references
        # because the estimator already holds a reference to the context.
        new_ctx._callbacks = getattr(estimator, "_skl_callbacks", [])
        new_ctx._estimator_name = estimator.__class__.__name__

        new_ctx._task_node = TaskNode(
            task_name=task_name,
            task_id=task_id,
            max_tasks=max_tasks,
            estimator_name=new_ctx._estimator_name,
        )

        if hasattr(estimator, "_parent_callback_ctx"):
            # This task is the root task of the estimator which itself corresponds to
            # a leaf task of a meta-estimator. Both tasks actually represent the same
            # task so we merge both task nodes into a single task node, attaching the
            # task tree of the sub-estimator to the task tree of the meta-estimator on
            # the way.
            parent_ctx = estimator._parent_callback_ctx
            new_ctx._task_node._merge_with(parent_ctx._task_node)
            new_ctx._estimator_depth = parent_ctx._estimator_depth + 1
        else:
            new_ctx._estimator_depth = 0

        return new_ctx

    @classmethod
    def _from_parent(cls, parent_context, *, task_name, task_id, max_tasks=1):
        """Private constructor to create a sub-context.

        Parameters
        ----------
        parent_context : CallbackContext instance
            The parent context of the new context.

        task_name : str
            The name of the task this context is responsible for.

        task_id : int
            The id of the task this context is responsible for.

        max_tasks : int, default=1
            The maximum number of tasks that can be siblings of the task this context is
            responsible for.
        """
        new_ctx = cls.__new__(cls)

        new_ctx._callbacks = parent_context._callbacks
        new_ctx._estimator_name = parent_context._estimator_name
        new_ctx._estimator_depth = parent_context._estimator_depth

        new_ctx._task_node = TaskNode(
            task_name=task_name,
            task_id=task_id,
            max_tasks=max_tasks,
            estimator_name=new_ctx._estimator_name,
        )

        # This task is a subtask of another task of a same estimator
        parent_context._task_node._add_child(new_ctx._task_node)

        return new_ctx

    def subcontext(self, task_name="", task_id=0, max_tasks=1):
        """Create a context for a subtask of the current task.

        Parameters
        ----------
        task_name : str, default=""
            The name of the subtask.

        task_id : int, default=0
            An identifier of the subtask. Usually a number between 0 and
            `max_tasks - 1`, but can be any identifier.

        max_tasks : int, default=1
            The maximum number of tasks that can be siblings of the subtask.
        """
        return CallbackContext._from_parent(
            parent_context=self,
            task_name=task_name,
            task_id=task_id,
            max_tasks=max_tasks,
        )

    def eval_on_fit_begin(self, estimator, *, data):
        """Evaluate the on_fit_begin method of the callbacks."""
        for callback in self._callbacks:
            # Only call the on_fit_begin method of callbacks that are not
            # propagated from a meta-estimator.
            if not (
                isinstance(callback, AutoPropagatedProtocol)
                and self._task_node.parent is not None
            ):
                callback.on_fit_begin(estimator, data=data)

        return self

    def eval_on_fit_iter_end(self, estimator, **kwargs):
        """Evaluate the on_fit_iter_end method of the callbacks."""
        return any(
            callback.on_fit_iter_end(estimator, self._task_node, **kwargs)
            for callback in self._callbacks
        )

    def eval_on_fit_end(self, estimator):
        """Evaluate the on_fit_end method of the callbacks."""
        for callback in self._callbacks:
            # Only call the on_fit_end method of callbacks that are not
            # propagated from a meta-estimator.
            if not (
                isinstance(callback, AutoPropagatedProtocol)
                and self._task_node.parent is not None
            ):
                callback.on_fit_end(estimator, task_node=self._task_node)

    def propagate_callbacks(self, sub_estimator):
        """Propagate the callbacks to a sub-estimator."""
        bad_callbacks = [
            callback.__class__.__name__
            for callback in getattr(sub_estimator, "_skl_callbacks", [])
            if isinstance(callback, AutoPropagatedProtocol)
        ]

        if bad_callbacks:
            raise TypeError(
                f"The sub-estimator ({sub_estimator.__class__.__name__}) of a"
                f" meta-estimator ({self._task_node.estimator_name}) can't have"
                f" auto-propagated callbacks ({bad_callbacks})."
                " Register them directly on the meta-estimator."
            )

        callbacks_to_propagate = [
            callback
            for callback in self._callbacks
            if isinstance(callback, AutoPropagatedProtocol)
            and (
                callback.max_estimator_depth is None
                or self._estimator_depth < callback.max_estimator_depth
            )
        ]

        if not callbacks_to_propagate:
            return self

        # We store the parent context in the sub-estimator to be able to merge the
        # task trees of the sub-estimator and the meta-estimator.
        sub_estimator._parent_callback_ctx = self

        sub_estimator.set_callbacks(
            getattr(sub_estimator, "_skl_callbacks", []) + callbacks_to_propagate
        )

        return self
