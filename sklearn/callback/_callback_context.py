# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import copy
import inspect
import uuid
import warnings
from contextlib import contextmanager

from sklearn.callback._base import AutoPropagatedCallback

# TODO(callbacks): move these explanations into a dedicated user guide.
#
# The computation tasks performed by an estimator during fit have an inherent tree
# structure, where each task can be decomposed into subtasks and so on. The root of the
# tree represents the whole fit task.
#
# Each loop in the estimator represents a parent task and each iteration of that loop
# represents a child task. To allow callbacks to be generic and reusable across
# estimators, the innermost tasks, i.e. the leaves of the task tree, correspond to
# operations on the full input data (or batch for incremental estimators).
#
# For instance, KMeans has two nested loops: the outer loop is controlled by `n_init`
# and the inner loop is controlled by `max_iter`. Its task tree looks like this:
#
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
#
# where each innermost `iter j` task corresponds the computation of the labels and
# centers for the full dataset.
#
# When the estimator is a meta-estimator, a task leaf usually corresponds to fitting
# a sub-estimator. Therefore this leaf and the root task of the sub-estimator actually
# represent the same task. In this case the leaf task of the meta-estimator and the root
# task of the sub-estimator are merged into a single task.
#
# For instance a `Pipeline` would have a task tree that looks like this:
#
# Pipeline fit
# ├── step 0 | preprocessor fit
# │   └── <insert preprocessor task tree here>
# └── step 1 | estimator fit
#     └── <insert estimator task tree here>
#
# Concretely, the tree structure is created dynamically and abstracted in an object
# named `CallbackContext`. There's a context for each task and the context is
# responsible for calling the callback hooks for its task and creating contexts for
# the child tasks.
#
# This `CallbackContext` is the object that has to be used in the implementation of an
# estimator to support callbacks. A context is created at the beginning of fit and
# then sub-contexts are created for each child task.
#
# class MyEstimator(CallbackSupportMixin, BaseEstimator):
#     def __init__(self, max_iter):
#         self.max_iter = max_iter
#
#     @with_callbacks
#     def fit(self, X, y):
#         callback_ctx = self._init_callback_context(max_subtasks=self.max_iter)
#         callback_ctx.call_on_fit_task_begin(X=X, y=y)
#
#         for i in range(self.max_iter):
#             subcontext = callback_ctx.subcontext(task_id=i).call_on_fit_task_begin(
#                X=X, y=y,
#             )
#
#             # Do something
#
#             subcontext.call_on_fit_task_end(X=X, y=y)
#
#         callback_ctx.call_on_fit_task_end(X=X, y=y)
#         return self
#
# It's also an object that is passed to the callback hooks to give them information
# about the task being executed and its position in the task tree.


# List of the parameters expected to be passed to call_on_fit_task_* (IN) and to be in
# the hooks signatures (OUT).
VALID_HOOK_PARAMS_IN = ["X", "y", "metadata", "reconstruction_attributes"]
VALID_HOOK_PARAMS_OUT = ["X", "y", "metadata", "fitted_estimator"]


class CallbackContext:
    """Task level context for the callbacks.

    This class is responsible for managing the callbacks and holding the tree structure
    of an estimator's tasks. Each instance corresponds to a task of the estimator.

    This class should not be instantiated directly, but through the
    `_init_callback_context` method of the estimator to create the root context or using
    the `subcontext` method of this class to create sub-contexts.

    These contexts are passed to the callback hooks to be able to keep track of the
    position of a task in the task tree from within the callbacks.

    Attributes
    ----------
    task_name : str
        The name of the task this context is responsible for.

    task_id : int
        The identifier of the task this context is responsible for. It uniquely
        identifies the task among its siblings.

    max_subtasks : int or None
        The maximum number of children tasks for this task. 0 means it's a leaf.
        None means the maximum number of subtasks is not known in advance.

    sequential_subtasks : bool
        Whether this context's subtasks are sequential. When True, children contexts'
        have consecutive integer task_ids starting from 0.

    estimator_name : str
        The name of the estimator that holds this context.

    parent : CallbackContext or None
        The parent context of this context. None if this context is the root.

    root_uuid : uuid.UUID instance
        The UUID of the root context. All contexts in the same task tree have the same
        root UUID that is used to identify the task tree itself.

    source_estimator_name : str or None
        The name of the estimator that holds the parent task this task was
        merged with. None if it was not merged with another context.

    source_task_name : str or None
        The task name of the parent task this task was merged with. None if it
        was not merged with another context.
    """

    @classmethod
    def _from_estimator(
        cls, estimator, task_name, task_id, max_subtasks, sequential_subtasks
    ):
        """Private constructor to create a root context.

        Parameters
        ----------
        estimator : estimator instance
            The estimator this context is responsible for.

        task_name : str
            The name of the root task.

        task_id : int
            Identifier for the root task.

        max_subtasks : int or None
            The maximum number of subtasks that can be children of the root task. None
            means the maximum number of subtasks is not known in advance. 0 means it's a
            leaf.

        sequential_subtasks : bool
            Whether the root context has sequential subtasks. If True, children contexts
            created via `subcontext` will have automatically assigned consecutive
            integer task_ids starting from 0.
        """
        new_ctx = cls.__new__(cls)

        # We don't store the estimator in the context to avoid circular references
        # because the estimator already holds a reference to the context.
        new_ctx._callbacks = getattr(estimator, "_skl_callbacks", [])
        new_ctx.estimator_name = estimator.__class__.__name__
        new_ctx.task_name = task_name
        new_ctx.task_id = task_id
        new_ctx.max_subtasks = max_subtasks
        new_ctx.sequential_subtasks = sequential_subtasks
        new_ctx.parent = None
        new_ctx.root_uuid = uuid.uuid4()
        new_ctx._children_map = {}
        new_ctx.source_estimator_name = None
        new_ctx.source_task_name = None

        if hasattr(estimator, "_parent_callback_ctx"):
            # This context's task is the root task of the estimator which itself
            # corresponds to a leaf task of a meta-estimator. Both tasks actually
            # represent the same task so we merge both tasks into a single task,
            # attaching the task tree of the sub-estimator to the task tree of
            # the meta-estimator on the way.
            parent_ctx = estimator._parent_callback_ctx
            new_ctx._merge_with(parent_ctx)
            new_ctx._propagation_depth = parent_ctx._propagation_depth + 1
        else:
            new_ctx._propagation_depth = 0

        return new_ctx

    @classmethod
    def _from_parent(
        cls, parent_context, *, task_name, task_id, max_subtasks, sequential_subtasks
    ):
        """Private constructor to create a sub-context.

        Parameters
        ----------
        parent_context : `CallbackContext` instance
            The parent context of the new context.

        task_name : str
            The name of the task this context is responsible for.

        task_id : int
            The identifier of the task this context is responsible for.

        max_subtasks : int or None
            The maximum number of tasks that can be children of the task this context is
            responsible for. 0 means it's a leaf. None means the maximum number of
            subtasks is not known in advance.

        sequential_subtasks : bool
            Whether this context's subtasks are sequential. If True, children contexts
            created via `subcontext` will have automatically assigned consecutive
            integer task_ids starting from 0.
        """
        new_ctx = cls.__new__(cls)

        new_ctx._callbacks = parent_context._callbacks
        new_ctx.estimator_name = parent_context.estimator_name
        new_ctx._propagation_depth = parent_context._propagation_depth
        new_ctx.task_name = task_name
        new_ctx.task_id = task_id
        new_ctx.max_subtasks = max_subtasks
        new_ctx.sequential_subtasks = sequential_subtasks
        new_ctx.root_uuid = parent_context.root_uuid
        new_ctx.parent = None
        new_ctx._children_map = {}
        new_ctx.source_estimator_name = None
        new_ctx.source_task_name = None

        # This task is a subtask of another task of a same estimator
        parent_context._add_child(new_ctx)

        return new_ctx

    def __iter__(self):
        """Pre-order depth-first traversal of the task tree."""
        yield self
        for context in self._children_map.values():
            yield from context

    def _add_child(self, child_context):
        """Add `child_context` as a child of this context."""
        if child_context.task_id in self._children_map:
            raise ValueError(
                f"Callback context {self.task_name} of estimator "
                f"{self.estimator_name} already has a child with "
                f"task_id={child_context.task_id}."
            )

        if (
            self.max_subtasks is not None
            and len(self._children_map) >= self.max_subtasks
        ):
            raise ValueError(
                f"Cannot add child to callback context {self.task_name} of estimator "
                f"{self.estimator_name} because it already has its maximum "
                f"number of children ({self.max_subtasks})."
            )

        self._children_map[child_context.task_id] = child_context
        child_context.parent = self

    def _merge_with(self, other_context):
        """Merge this context with `other_context`.

        This method is called on a sub-estimator's root task to merge it with a
        meta-estimator's leaf task. The sub-estimator's task tree is therefore attached
        to the meta-estimator's task tree. The root node of the sub-estimator's task
        tree and the leaf node of the meta-estimator's task tree are both represented
        by a single node in this combined task tree.
        """
        if other_context.max_subtasks != 0:
            raise ValueError(
                f"Cannot merge callback context (task {self.task_name!r} of estimator "
                f"{self.estimator_name}) with callback context "
                f"(task {other_context.task_name!r} of estimator "
                f"{other_context.estimator_name}) because the latter is not a leaf."
            )

        # Set the parent of the sub-estimator's root context to the parent of the
        # meta-estimator's leaf context
        self.parent = other_context.parent
        self.task_id = other_context.task_id
        self.root_uuid = other_context.root_uuid
        other_context.parent._children_map[self.task_id] = self

        # Keep information about the context it was merged with
        self.source_task_name = other_context.task_name
        self.source_estimator_name = other_context.estimator_name

    def subcontext(
        self, task_name="", task_id=None, max_subtasks=0, sequential_subtasks=True
    ):
        """Create a context for a subtask of the current task.

        Parameters
        ----------
        task_name : str, default=""
            The name of the subtask.

        task_id : int or None, default=None
            An identifier of the subtask. It must be distinct from the task_ids of its
            siblings. If None, task_id is automatically set to the next available
            integer task_id.

        max_subtasks : int or None, default=0
            The maximum number of tasks that can be children of the subtask. 0 means
            it's a leaf. None means the maximum number of subtasks is not known in
            advance.

        sequential_subtasks : bool, default=True
            Whether the new context's subtasks are sequential. If True, children
            contexts of the new context, created via `subcontext`, will have
            automatically assigned consecutive integer task_ids starting from 0.
        """
        if self.sequential_subtasks and task_id is not None:
            raise ValueError(
                f"task_id for {self.estimator_name} {task_name} must be None if "
                f"sequential_subtasks is True for {self.task_name}."
            )
        if not self.sequential_subtasks and task_id is None:
            raise ValueError(
                f"task_id for {self.estimator_name} {task_name} must be provided "
                f"if sequential_subtasks is False for {self.task_name}."
            )
        if task_id is None:
            task_id = len(self._children_map)

        return CallbackContext._from_parent(
            parent_context=self,
            task_name=task_name,
            task_id=task_id,
            max_subtasks=max_subtasks,
            sequential_subtasks=sequential_subtasks,
        )

    def _call_hooks(self, estimator, hook_name, **kwargs):
        """Helper to call the hook of all callbacks with their respective arguments.

        Provide the right arguments to each hook by inspecting their signatures. Any
        value that is a callable is replaced by what it returns to allow lazy loading of
        the arguments.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling the callback hook.

        hook_name : str
            Name of the callback hook to call.

        **kwargs: dict
            Optional keyword arguments passed to the callback context.

        Returns
        -------
        result : bool
            True if any hook call returned True. False otherwise.
        """
        result = False

        # Keep a cache of the evaluated args to evaluate them only once.
        evaluated_args = {}

        for callback in self._callbacks:
            if callback in getattr(self, "_propagated_callbacks", []):
                # Only call the `on_fit_task_end` hook of callbacks that are not
                # propagated. For propagated callbacks, the hook will be called by the
                # sub-estimator's root context (both represent the same task).
                continue

            signature = inspect.signature(getattr(callback, hook_name))
            params_names = {
                p.name
                for p in signature.parameters.values()
                if p.kind == p.KEYWORD_ONLY
            }
            if diff := set(params_names) - set(VALID_HOOK_PARAMS_OUT):
                raise TypeError(
                    f"Hook {hook_name} of the callback {callback.__class__.__name__} "
                    f"has parameters that are not valid: {diff}. The valid parameters "
                    f"are: {VALID_HOOK_PARAMS_OUT}."
                )

            args_to_pass = {}
            for param_name in params_names:
                if param_name not in evaluated_args:
                    # Special case: "reconstruction_attributes" is not directly passed
                    # to the hook. A ready to predict/transform estimator is created
                    # from these attributes and passed to the hook as "fitted_estimator"
                    if param_name == "fitted_estimator":
                        attrs = kwargs.get("reconstruction_attributes", None)
                        attrs = attrs() if callable(attrs) else attrs
                        new_est = (
                            _from_reconstruction_attributes(estimator, attrs)
                            if attrs is not None
                            else None
                        )
                        evaluated_args["fitted_estimator"] = new_est
                    else:
                        val = kwargs.get(param_name, None)
                        val = val() if callable(val) else val
                        evaluated_args[param_name] = val

                args_to_pass[param_name] = evaluated_args[param_name]

            result |= bool(
                getattr(callback, hook_name)(estimator, self, **args_to_pass)
            )

        return result

    def call_on_fit_task_begin(
        self,
        *,
        estimator,
        X=None,
        y=None,
        metadata=None,
        reconstruction_attributes=None,
    ):
        """Call the `on_fit_task_begin` hook of the callbacks.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling the callback hook.

        X : array-like or None, default=None
            The training data of the current task.

        y : array-like or None, default=None
            The training targets of the current task.

        metadata : dict or None, default=None
            A dictionary containing training metadata for the current task.

        reconstruction_attributes : dict or None, default=None
            A dictionary of the sufficient fitted attributes needed to construct a
            `fitted_estimator` from the current state of the estimator, i.e. an
            estimator instance ready to predict, transform, etc ... as if the fit had
            stopped at the beginning of this task. The `fitted_estimator` is the
            object that will be passed to the callbacks, if required.
        """
        self._call_hooks(
            estimator,
            hook_name="on_fit_task_begin",
            X=X,
            y=y,
            metadata=metadata,
            reconstruction_attributes=reconstruction_attributes,
        )
        return self

    def call_on_fit_task_end(
        self,
        *,
        estimator,
        X=None,
        y=None,
        metadata=None,
        reconstruction_attributes=None,
    ):
        """Call the `on_fit_task_end` hook of the callbacks.

        Parameters
        ----------
        estimator : estimator instance
            The estimator calling the callback hook.

        X : array-like or None, default=None
            The training data of the current task.

        y : array-like or None, default=None
            The training targets of the current task.

        metadata : dict or None, default=None
            A dictionary containing training metadata of the current task.

        reconstruction_attributes : dict or None, default=None
            A dictionary of the sufficient fitted attributes needed to construct a
            `fitted_estimator` from the current state of the estimator, i.e. an
            estimator instance ready to predict, transform, etc ... as if the fit had
            stopped at the end of this task. The `fitted_estimator` is the object
            that will be passed to the callbacks, if required.

        Returns
        -------
        stop : bool
            Whether or not to stop the current level of iterations at this end of this
            task.
        """
        return self._call_hooks(
            estimator,
            hook_name="on_fit_task_end",
            X=X,
            y=y,
            metadata=metadata,
            reconstruction_attributes=reconstruction_attributes,
        )

    @contextmanager
    def propagate_callback_context(self, sub_estimator):
        """Propagate the context and callbacks to a sub-estimator.

        Clear the propagated callbacks from the sub-estimator on exit.

        Only auto-propagated callbacks are propagated to the sub-estimator. An error is
        raised if the sub-estimator already holds auto-propagated callbacks.

        The sub-estimator receives this context as an attribute named
        `_parent_callback_ctx` so that the meta-estimator's task tree can be merged with
        the sub-estimator's one.

        Parameters
        ----------
        sub_estimator : estimator instance
            The estimator to propagate the callbacks and context to.
        """
        bad_callbacks = [
            callback.__class__.__name__
            for callback in getattr(sub_estimator, "_skl_callbacks", [])
            if isinstance(callback, AutoPropagatedCallback)
        ]
        if bad_callbacks:
            raise TypeError(
                f"The sub-estimator ({sub_estimator.__class__.__name__}) of a"
                f" meta-estimator ({self.estimator_name}) can't have"
                f" auto-propagated callbacks ({bad_callbacks})."
                " Register them directly on the meta-estimator."
            )

        # We store the parent context in the sub-estimator to be able to merge the task
        # trees of the sub-estimator and the meta-estimator. We want to link the task
        # trees even if there is no callback to propagate, as the sub-estimators might
        # have non auto-propagated callbacks, which would need to have access to the
        # whole tree.
        sub_estimator._parent_callback_ctx = self

        callbacks_to_propagate = [
            callback
            for callback in self._callbacks
            if isinstance(callback, AutoPropagatedCallback)
            and (
                callback.max_propagation_depth is None
                or self._propagation_depth < callback.max_propagation_depth
            )
        ]
        if callbacks_to_propagate and not hasattr(sub_estimator, "set_callbacks"):
            warnings.warn(
                f"The estimator {sub_estimator.__class__.__name__} does not support "
                f"callbacks. The callbacks attached to {self.estimator_name} will not "
                f"be propagated to this estimator."
            )
            callbacks_to_propagate = []

        if callbacks_to_propagate:
            self._propagated_callbacks = callbacks_to_propagate
            curr_callbacks = getattr(sub_estimator, "_skl_callbacks", [])
            sub_estimator.set_callbacks(*(curr_callbacks + callbacks_to_propagate))

        try:
            yield
        finally:
            if callbacks_to_propagate:
                kept_callbacks = [
                    cb
                    for cb in sub_estimator._skl_callbacks
                    if cb not in callbacks_to_propagate
                ]
                sub_estimator.set_callbacks(*kept_callbacks)
            del sub_estimator._parent_callback_ctx


def _from_reconstruction_attributes(estimator, reconstruction_attributes):
    """Return a copy of the estimator as if it was fitted.

    Parameters
    ----------
    estimator : estimator instance
        The estimator from which to make a ready-to-be-evaluated copy.

    reconstruction_attributes : dict
        A dictionary containing the necessary attributes to create a working
        fitted estimator from this instance.

    Returns
    -------
    fitted_estimator : estimator instance
        The fitted copy of this estimator.
    """
    new_estimator = copy.copy(estimator)  # TODO(callbacks) copy / deepcopy / clone ?
    for key, val in reconstruction_attributes.items():
        setattr(new_estimator, key, val)
    return new_estimator


def get_context_path(context):
    """Helper function to get the path from the root context down to a given context.

    Parameters
    ----------
    context : `CallbackContext` instance
        The context to get the path to.

    Returns
    -------
    list of `CallbackContext` instances
        The list of the ancestors (itself included) of the given context. The list is
        ordered from the root context to the given context.
    """
    return (
        [context]
        if context.parent is None
        else get_context_path(context.parent) + [context]
    )
