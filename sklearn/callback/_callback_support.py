# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import functools
from contextlib import contextmanager
from threading import Lock

from joblib.externals.loky.backend import get_context

from sklearn.callback._base import AutoPropagatedCallback, FitCallback
from sklearn.callback._callback_context import CallbackContext


class _CallbackManagerState:
    manager = None
    lock = Lock()


def get_callback_manager():
    """Return the global multiprocessing manager dedicated to callbacks.

    The manager is initialized lazily on first access and reused afterwards.
    """
    if _CallbackManagerState.manager is None:
        with _CallbackManagerState.lock:
            if _CallbackManagerState.manager is None:
                _CallbackManagerState.manager = get_context().Manager()

    return _CallbackManagerState.manager


class CallbackSupportMixin:
    """Mixin class to add callback support to an estimator."""

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
        if not all(
            # isinstance for a Protocol returns True for classes too (not only
            # instances). Hence the additional check for classes.
            isinstance(callback, FitCallback) and not isinstance(callback, type)
            for callback in callbacks
        ):
            raise TypeError(
                "callbacks must be instances following the FitCallback protocol."
            )

        if callbacks:
            self._skl_callbacks = list(callbacks)
        else:
            del self._skl_callbacks

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

        # setup callbacks
        for callback in getattr(self, "_skl_callbacks", []):
            # Only call the setup hook of callbacks that are not propagated from a
            # meta-estimator.
            if not (
                isinstance(callback, AutoPropagatedCallback)
                and hasattr(self, "_parent_callback_ctx")
            ):
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
            for callback in getattr(estimator, "_skl_callbacks", []):
                # Only call the teardown hook of callbacks that are not propagated from
                # a meta-estimator.
                if not (
                    isinstance(callback, AutoPropagatedCallback)
                    and hasattr(estimator, "_parent_callback_ctx")
                ):
                    callback.teardown(
                        estimator=estimator, context=estimator._callback_fit_ctx
                    )

            del estimator._callback_fit_ctx


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
