# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import functools
import inspect
from contextlib import contextmanager

from sklearn.callback._base import AutoPropagatedCallback, FitCallback
from sklearn.callback._callback_context import CallbackContext


def validate_callbacks(callbacks):
    """Validate that callbacks follow the protocol and correctly implement hooks."""
    for callback in callbacks:
        callback_name = callback.__class__.__name__

        # isinstance for a Protocol returns True for classes too (not only instances).
        # Hence the additional check for classes.
        if not isinstance(callback, FitCallback) or isinstance(callback, type):
            raise TypeError(
                f"Callbacks must be instances following the FitCallback protocol. "
                f"Got {callback_name}."
            )

        # isinstance for Protocols does not verify the methods' signatures.
        for hook_name in ("setup", "teardown", "on_fit_task_begin", "on_fit_task_end"):
            hook = getattr(callback, hook_name)
            params = list(inspect.signature(hook).parameters.values())

            positional = [p.name for p in params if p.kind != p.KEYWORD_ONLY]
            expected_positional = ["estimator", "context"]
            if positional != expected_positional:
                raise TypeError(
                    f"Hook {hook_name!r} of the callback {callback_name} must have "
                    f"exactly the positional parameters {expected_positional}, got "
                    f"{positional}."
                )

            if hook_name in ("setup", "teardown"):
                if len(params) > 2:
                    raise TypeError(
                        f"The signature of the {hook_name!r} hook of the callback "
                        f"{callback_name} must be {hook_name}(self, estimator, context)"
                    )
                continue

            kwarg_only = {p.name for p in params if p.kind == p.KEYWORD_ONLY}
            valid_kwargs = {"X", "y", "metadata", "fitted_estimator"}
            if invalid := kwarg_only - valid_kwargs:
                raise TypeError(
                    f"Hook {hook_name!r} of the callback {callback_name} has "
                    f"parameters that are not valid: {invalid}. The valid parameters "
                    f"are: {valid_kwargs}."
                )


class CallbackSupportMixin:
    """Mixin class to add callback support to an estimator.

    .. document the private method
    .. automethod:: _init_callback_context
    """

    def set_callbacks(self, *callbacks):
        """Set callbacks for the estimator.

        Parameters
        ----------
        *callbacks : callback instances
            The callbacks to set.

        Returns
        -------
        self : estimator instance
            The estimator instance itself.
        """
        validate_callbacks(callbacks)
        return self._set_callbacks(callbacks)

    def _set_callbacks(self, callbacks):
        """set_callbacks without validation."""
        if callbacks:
            self._skl_callbacks = list(callbacks)
        else:
            self.__dict__.pop("_skl_callbacks", None)

        return self

    def _init_callback_context(
        self, task_name="fit", task_id=0, max_subtasks=0, sequential_subtasks=True
    ):
        """Initialize the callback context for the estimator and setup its callbacks.

        This method should be called once, at the beginning of the fit method.

        It will only setup callbacks that are not propagated from a meta-estimator.

        Parameters
        ----------
        task_name : str, default="fit"
            The name of the root task.

        task_id : int, default=0
            Identifier for the root task.

        max_subtasks : int or None, default=0
            The maximum number of subtasks that can be children of the root task. None
            means the maximum number of subtasks is not known in advance. 0 means it's a
            leaf.

        sequential_subtasks : bool, default=True
            Whether the root context has sequential subtasks. If True, children contexts
            of the root context, created via `subcontext`, will have automatically
            assigned consecutive integer task_ids starting from 0.

        Returns
        -------
        callback_fit_ctx : CallbackContext
            The root callback context for the estimator.
        """
        self._callback_fit_ctx = CallbackContext._from_estimator(
            estimator=self,
            task_name=task_name,
            task_id=task_id,
            max_subtasks=max_subtasks,
            sequential_subtasks=sequential_subtasks,
        )

        # Setup callbacks. We store callbacks for which setup has started in order to
        # only tear those down after fit.
        self._skl_callbacks_to_teardown = []
        for callback in getattr(self, "_skl_callbacks", []):
            # Only call the setup hook of callbacks that are not propagated from a
            # meta-estimator.
            if not (
                isinstance(callback, AutoPropagatedCallback)
                and hasattr(self, "_parent_callback_ctx")
            ):
                self._skl_callbacks_to_teardown.append(callback)
                callback.setup(estimator=self, context=self._callback_fit_ctx)

        return self._callback_fit_ctx


@contextmanager
def callback_management_context(estimator):
    """Context manager to manage callback lifecycle around estimator fit.

    The context manager is responsible for calling the callbacks `teardown` hook in a
    `try finally` block, which guarantees that callbacks teardown will always be called,
    whether the estimator's fit exits successfully or not.

    Parameters
    ----------
    estimator : estimator instance
        The estimator being fitted.

    Yields
    ------
    None.
    """
    try:
        yield
    finally:
        if hasattr(estimator, "_callback_fit_ctx"):
            teardown_errors = []
            for callback in estimator._skl_callbacks_to_teardown:
                try:
                    callback.teardown(
                        estimator=estimator, context=estimator._callback_fit_ctx
                    )
                except Exception as exc:
                    teardown_errors.append(exc)

            del estimator._skl_callbacks_to_teardown
            del estimator._callback_fit_ctx
            if len(teardown_errors) == 1:
                raise teardown_errors[0]
            if teardown_errors:
                raise ExceptionGroup(
                    "The following callback teardown errors occurred",
                    teardown_errors,
                )


def with_callbacks(method):
    """Decorator to run the method within a callback context manager.

    This decorator is responsible for calling the callbacks `teardown` hooks of
    callbacks in a `try finally` block, which guarantees that callbacks teardown will
    always be called, whether the estimator's method exits successfully or not.

    It will only teardown callbacks that are not propagated from a meta-estimator.

    Parameters
    ----------
    method : method
        The method to decorate.

    Returns
    -------
    decorated_method : method
        The decorated method.
    """

    @functools.wraps(method)
    def callback_managed_method(estimator, *args, **kwargs):
        with callback_management_context(estimator):
            return method(estimator, *args, **kwargs)

    return callback_managed_method
