# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.callback import AutoPropagatedProtocol

# The computations performed by an estimator have an inherent tree structure, with
# each node representing a task. Each loop in the estimator represents a parent task
# node and each iteration of that loop represents a child task node. Usually the root
# task node represents the whole fit task and leaves the innermost loop iterations.

# For instance, a KMeans estimator has two nested loops: the outer loop is controlled
# by `n_init` and the inner loop is controlled by `max_iter`. Its task tree looks like
# this:

# KMeans fit
# ├── init 0
# │   ├── iter 0
# │   ├── iter 1
# │   ├── ...
# │   └── iter n
# ├── init 1
# │   ├── iter 0
# │   ├── ...
# │   └── iter n
# └── init 2
#     ├── iter 0
#     ├── ...
#     └── iter n

# When the estimator is a meta-estimator, a task leaf usually correspond to fitting
# a sub-estimator. Therefore this leaf and the root task of the sub-estimator actually
# represent the same task. In this case the leaf task node of the meta-estimator and
# the root task node of the sub-estimator are merged into a single task node.

# For instance a `Pipeline` would have a task tree that looks like this:

# Pipeline fit
# ├── step 0 | preprocessor fit
# │   └── <insert preprocessor task tree here>
# └── step 1 | estimator fit
#     └── <insert estimator task tree here>


# The task tree is created dynamically, by initializing it with a root task node and
# then adding the child task nodes as the fitting process goes on.
class CallbackContext:
    """Task level context for the callbacks.

    This class is responsible for managing the callbacks and holding the tree structure
    of an estimator's tasks.

    Instances of this class should be created using the `init_callback_context` method
    of the estimator.

    Attributes
    ----------
    children_map : dict
        A mapping from the task_id of a child to the child node
        `{task_id: CallbackContext}`. For a leaf, it's an empty dictionary.

    depth : int
        The depth of this node in the task tree.

    estimator_name : str
        The estimator name of this node.

    max_subtasks : int or None
        The maximum number of subtasks of this node of the task tree. 0 means it's a
        leaf. None means the maximum number of subtasks is not known in advance.

    parent : CallbackContext instance or None
        The parent node of the task tree. None means this is the root.

        Note that it is a dynamic attribute since the root task of an estimator can
        become an intermediate node of a meta-estimator.

    prev_estimator_name : str or None
        The estimator name of the node this node was merged with. None if it was not
        merged with another node.

    prev_task_name : str
        The task name of the node this node was merged with. None if it was not
        merged with another node.

    task_id : int
        The id of this task node.

    task_info : dict
        Dictionary description of the corresponding computation task.

    task_name : str
        The task name of this node.
    """

    @classmethod
    def _from_estimator(cls, estimator, *, task_name, task_id, max_subtasks=None):
        """Private constructor to create a root context.

        Parameters
        ----------
        estimator : estimator instance
            The estimator this context is responsible for.

        task_name : str
            The name of the task this context is responsible for.

        task_id : int
            The id of the task this context is responsible for.

        max_subtasks : int or None, default=None
            The maximum number of subtasks of this node of the task tree. 0 means it's a
            leaf. None means the maximum number of subtasks is not known in advance.
        """
        new_ctx = cls.__new__(cls)

        # We don't store the estimator in the context to avoid circular references
        # because the estimator already holds a reference to the context.
        new_ctx._callbacks = getattr(estimator, "_skl_callbacks", [])
        new_ctx.estimator_name = estimator.__class__.__name__
        new_ctx.task_name = task_name
        new_ctx.task_id = task_id
        new_ctx.parent = None
        new_ctx.children_map = {}
        new_ctx.max_subtasks = max_subtasks
        new_ctx.prev_estimator_name = None
        new_ctx.prev_task_name = None

        if hasattr(estimator, "_parent_callback_ctx"):
            # This task is the root task of the estimator which itself corresponds to
            # a leaf task of a meta-estimator. Both tasks actually represent the same
            # task so we merge both task nodes into a single task node, attaching the
            # task tree of the sub-estimator to the task tree of the meta-estimator on
            # the way.
            parent_ctx = estimator._parent_callback_ctx
            new_ctx._merge_with(parent_ctx)
            new_ctx._estimator_depth = parent_ctx._estimator_depth + 1
        else:
            new_ctx._estimator_depth = 0

        return new_ctx

    @classmethod
    def _from_parent(cls, parent_context, *, task_name, task_id, max_subtasks=None):
        """Private constructor to create a sub-context.

        Parameters
        ----------
        parent_context : CallbackContext instance
            The parent context of the new context.

        task_name : str
            The name of the task this context is responsible for.

        task_id : int
            The id of the task this context is responsible for.

        max_subtasks : int or None, default=None
            The maximum number of tasks that can be children of the task this context is
            responsible for. 0 means it's a leaf. None means the maximum number of
            subtasks is not known in advance.
        """
        new_ctx = cls.__new__(cls)

        new_ctx._callbacks = parent_context._callbacks
        new_ctx.estimator_name = parent_context.estimator_name
        new_ctx._estimator_depth = parent_context._estimator_depth
        new_ctx.task_name = task_name
        new_ctx.task_id = task_id
        new_ctx.parent = None
        new_ctx.children_map = {}
        new_ctx.max_subtasks = max_subtasks
        new_ctx.prev_estimator_name = None
        new_ctx.prev_task_name = None

        # This task is a subtask of another task of a same estimator
        parent_context._add_child(new_ctx)

        return new_ctx

    @property
    def task_info(self):
        """Information about the corresponding computation task. the keys are

        - depth : int
            The depth of the task in the task tree.
        - estimator_name : str
            The name of the estimator.
        - max_subtasks : int
            The maximum number of children tasks for this task.
        - parent : dict or None
            The task_info property of this task's parent, None if this task is the root.
        - prev_estimator_name : str or None
            The name of the meta estimator, in the case that this task was merged to a
            meta estimator's subtask, None otherwise.
        - prev_task_name : str
            The name of the meta estimator's subtask this task was merged to, None
            if it was not merged to a task.
        - task_id : int
            The identifier of the task.
        - task_name : str
            The name of the task.
        """
        return {
            "estimator_name": self.estimator_name,
            "depth": self.depth,
            "task_name": self.task_name,
            "task_id": self.task_id,
            "max_subtasks": self.max_subtasks,
            "prev_estimator_name": self.prev_estimator_name,
            "prev_task_name": self.prev_task_name,
            "parent": None if self.parent is None else self.parent.task_info,
        }

    @property
    def depth(self):
        """The depth of this node in the task tree."""
        return 0 if self.parent is None else self.parent.depth + 1

    def __iter__(self):
        """Pre-order depth-first traversal of the task tree."""
        yield self
        for context in self.children_map.values():
            yield from context

    def _add_child(self, context):
        if context.task_id in self.children_map:
            raise ValueError(
                f"Callback context {self.task_name} of estimator {self.estimator_name} "
                f"already has a child with task_id={context.task_id}."
            )

        if len(self.children_map) == self.max_subtasks:
            raise ValueError(
                f"Cannot add child to callback context {self.task_name} of estimator "
                f"{self.estimator_name} because it already has its maximum "
                f"number of children ({self.max_subtasks})."
            )

        self.children_map[context.task_id] = context
        context.parent = self

    def _merge_with(self, context):
        # Set the parent of the sub-estimator's root context node to the parent
        # of the meta-estimator's leaf context node
        self.parent = context.parent
        self.task_id = context.task_id
        self.max_subtasks = context.max_subtasks
        context.parent.children_map[self.task_id] = self

        # Keep information about the node it was merged with
        self.prev_task_name = context.task_name
        self.prev_estimator_name = context.estimator_name

    def subcontext(self, task_name="", task_id=0, max_subtasks=None):
        """Create a context for a subtask of the current task.

        Parameters
        ----------
        task_name : str, default=""
            The name of the subtask.

        task_id : int, default=0
            An identifier of the subtask. Usually a number between 0 and the maximum
            number of sibling contexts, but can be any identifier.

        max_subtasks : int or None, default=None
            The maximum number of tasks that can be children of the subtask. 0 means
            it's a leaf. None means the maximum number of subtasks is not known in
            advance.
        """
        return CallbackContext._from_parent(
            parent_context=self,
            task_name=task_name,
            task_id=task_id,
            max_subtasks=max_subtasks,
        )

    def eval_on_fit_begin(self, estimator, *, data):
        """Evaluate the _on_fit_begin method of the callbacks.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.

        data : dict
            Dictionary containing the training and validation data. The possible
            keys are "X_train", "y_train", "sample_weight_train", "X_val", "y_val"
            and "sample_weight_val".
        """
        for callback in self._callbacks:
            # Only call the on_fit_begin method of callbacks that are not
            # propagated from a meta-estimator.
            if not (
                isinstance(callback, AutoPropagatedProtocol) and self.parent is not None
            ):
                callback._on_fit_begin(estimator, data=data)

        return self

    def eval_on_fit_iter_end(self, estimator, **kwargs):
        """Evaluate the _on_fit_iter_end method of the callbacks.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.

        **kwargs : dict
            arguments passed to the callback. Possible keys are

            - data: dict
                Dictionary containing the training and validation data. The keys are
                "X_train", "y_train", "sample_weight_train", "X_val", "y_val",
                "sample_weight_val". The values are the corresponding data. If a key is
                missing, the corresponding value is None.

            - stopping_criterion: float
                Usually iterations stop when `stopping_criterion <= tol`.
                This is only provided at the innermost level of iterations.

            - tol: float
                Tolerance for the stopping criterion.
                This is only provided at the innermost level of iterations.

            - from_reconstruction_attributes: estimator instance
                A ready to predict, transform, etc ... estimator as if the fit stopped
                at this node. Usually it's a copy of the caller estimator with the
                necessary attributes set but it can sometimes be an instance of another
                class (e.g. LogisticRegressionCV -> LogisticRegression)

            - fit_state: dict
                Model specific quantities updated during fit. This is not meant to be
                used by generic callbacks but by a callback designed for a specific
                estimator instead.

        Returns
        -------
        stop : bool
            Whether or not to stop the current level of iterations at this task node.
        """
        return any(
            callback._on_fit_iter_end(estimator, self.task_info, **kwargs)
            for callback in self._callbacks
        )

    def eval_on_fit_end(self, estimator):
        """Evaluate the _on_fit_end method of the callbacks.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling this callback hook.
        """
        for callback in self._callbacks:
            # Only call the on_fit_end method of callbacks that are not
            # propagated from a meta-estimator.
            if not (
                isinstance(callback, AutoPropagatedProtocol) and self.parent is not None
            ):
                callback._on_fit_end(estimator, task_info=self.task_info)

    def propagate_callbacks(self, sub_estimator):
        """Propagate the callbacks to a sub-estimator.

        Parameters
        ----------
        sub_estimator : estimator instance
            The estimator to which the callbacks should be propagated.
        """
        bad_callbacks = [
            callback.__class__.__name__
            for callback in getattr(sub_estimator, "_skl_callbacks", [])
            if isinstance(callback, AutoPropagatedProtocol)
        ]

        if bad_callbacks:
            raise TypeError(
                f"The sub-estimator ({sub_estimator.__class__.__name__}) of a"
                f" meta-estimator ({self.estimator_name}) can't have"
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
