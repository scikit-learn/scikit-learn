# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import functools
import sys
from contextlib import contextmanager
from threading import Lock

from joblib.externals.loky.backend import get_context

from sklearn.callback._base import AutoPropagatedCallback, FitCallback
from sklearn.callback._callback_context import CallbackContext


class _CallbackManagerState:
    manager = None
    lock = Lock()


def _is_debugger_active():
    """Return True if a debugger (e.g., debugpy or pydevd) is active."""
    return "debugpy" in sys.modules or "pydevd" in sys.modules


def get_callback_manager():
    """Return the global multiprocessing manager dedicated to callbacks.

    The manager is initialized lazily on first access and reused afterwards.
    """
    if _CallbackManagerState.manager is None:
        with _CallbackManagerState.lock:
            if _CallbackManagerState.manager is None:
                if _is_debugger_active():
                    # loky may spawn subprocesses that break debugger sessions
                    # (e.g., debugpy), leading to EOFError. Use standard
                    # multiprocessing "spawn" context when a debugger is active.
                    import multiprocessing

                    ctx = multiprocessing.get_context("spawn")
                else:
                    ctx = get_context()  # default loky context

                _CallbackManagerState.manager = ctx.Manager()

    return _CallbackManagerState.manager


class CallbackSupportMixin:
    """Mixin class to add callback support to an estimator."""

    def set_callbacks(self, callbacks):
        """Set callbacks for the estimator.

        Parameters
        ----------
        callbacks : callback or list of callbacks
            The callbacks to set.

        Returns
        -------
        self : estimator instance
            The estimator instance itself.
        """
        if not isinstance(callbacks, list):
            callbacks = [callbacks]

        if not all(
            isinstance(callback, FitCallback) and not isinstance(callback, type)
            for callback in callbacks
        ):
            raise TypeError(
                "callbacks must be instances following the FitCallback protocol."
            )

        self._skl_callbacks = callbacks
        return self

    def _init_callback_context(self, task_name="fit", task_id=0, max_subtasks=0):
        """Initialize the callback context for the estimator and setup its callbacks."""
        self._callback_fit_ctx = CallbackContext._from_estimator(
            estimator=self,
            task_name=task_name,
            task_id=task_id,
            max_subtasks=max_subtasks,
        )

        for callback in getattr(self, "_skl_callbacks", []):
            if not (
                isinstance(callback, AutoPropagatedCallback)
                and hasattr(self, "_parent_callback_ctx")
            ):
                callback.setup(self._callback_fit_ctx)

        return self._callback_fit_ctx


@contextmanager
def callback_management_context(estimator):
    """Context manager to manage callback lifecycle around estimator fit."""
    try:
        yield
    finally:
        if hasattr(estimator, "_callback_fit_ctx"):
            for callback in getattr(estimator, "_skl_callbacks", []):
                if not (
                    isinstance(callback, AutoPropagatedCallback)
                    and hasattr(estimator, "_parent_callback_ctx")
                ):
                    callback.teardown(estimator._callback_fit_ctx)

            del estimator._callback_fit_ctx
            if hasattr(estimator, "_parent_callback_ctx"):
                del estimator._parent_callback_ctx


def with_callbacks(method):
    """Decorator to run the method within a callback context manager."""

    @functools.wraps(method)
    def callback_managed_method(estimator, *args, **kwargs):
        with callback_management_context(estimator):
            return method(estimator, *args, **kwargs)

    return callback_managed_method
