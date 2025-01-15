"""Customizations of :mod:`joblib` and :mod:`threadpoolctl` tools for scikit-learn
usage.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import functools
import warnings
from functools import update_wrapper

import joblib
from threadpoolctl import ThreadpoolController

from .._config import config_context, get_config

# Global threadpool controller instance that can be used to locally limit the number of
# threads without looping through all shared libraries every time.
# It should not be accessed directly and _get_threadpool_controller should be used
# instead.
_threadpool_controller = None


def _with_config_and_warning_filters(delayed_func, config, warning_filters):
    """Helper function that intends to attach a config to a delayed function."""
    if hasattr(delayed_func, "with_config_and_warning_filters"):
        return delayed_func.with_config_and_warning_filters(config, warning_filters)
    else:
        warnings.warn(
            (
                "`sklearn.utils.parallel.Parallel` needs to be used in "
                "conjunction with `sklearn.utils.parallel.delayed` instead of "
                "`joblib.delayed` to correctly propagate the scikit-learn "
                "configuration to the joblib workers."
            ),
            UserWarning,
        )
        return delayed_func


class Parallel(joblib.Parallel):
    """Tweak of :class:`joblib.Parallel` that propagates the scikit-learn configuration.

    This subclass of :class:`joblib.Parallel` ensures that the active configuration
    (thread-local) of scikit-learn is propagated to the parallel workers for the
    duration of the execution of the parallel tasks.

    The API does not change and you can refer to :class:`joblib.Parallel`
    documentation for more details.

    .. versionadded:: 1.3
    """

    def __call__(self, iterable):
        """Dispatch the tasks and return the results.

        Parameters
        ----------
        iterable : iterable
            Iterable containing tuples of (delayed_function, args, kwargs) that should
            be consumed.

        Returns
        -------
        results : list
            List of results of the tasks.
        """
        # Capture the thread-local scikit-learn configuration at the time
        # Parallel.__call__ is issued since the tasks can be dispatched
        # in a different thread depending on the backend and on the value of
        # pre_dispatch and n_jobs.
        config = get_config()
        warning_filters = warnings.filters
        iterable_with_config_and_warning_filters = (
            (
                _with_config_and_warning_filters(delayed_func, config, warning_filters),
                args,
                kwargs,
            )
            for delayed_func, args, kwargs in iterable
        )
        return super().__call__(iterable_with_config_and_warning_filters)


# remove when https://github.com/joblib/joblib/issues/1071 is fixed
def delayed(function):
    """Decorator used to capture the arguments of a function.

    This alternative to `joblib.delayed` is meant to be used in conjunction
    with `sklearn.utils.parallel.Parallel`. The latter captures the scikit-
    learn configuration by calling `sklearn.get_config()` in the current
    thread, prior to dispatching the first task. The captured configuration is
    then propagated and enabled for the duration of the execution of the
    delayed function in the joblib workers.

    .. versionchanged:: 1.3
       `delayed` was moved from `sklearn.utils.fixes` to `sklearn.utils.parallel`
       in scikit-learn 1.3.

    Parameters
    ----------
    function : callable
        The function to be delayed.

    Returns
    -------
    output: tuple
        Tuple containing the delayed function, the positional arguments, and the
        keyword arguments.
    """

    @functools.wraps(function)
    def delayed_function(*args, **kwargs):
        return _FuncWrapper(function), args, kwargs

    return delayed_function


class _FuncWrapper:
    """Load the global configuration before calling the function."""

    def __init__(self, function):
        self.function = function
        update_wrapper(self, self.function)

    def with_config_and_warning_filters(self, config, warning_filters):
        self.config = config
        self.warning_filters = warning_filters
        return self

    def __call__(self, *args, **kwargs):
        config = getattr(self, "config", {})
        warning_filters = getattr(self, "warning_filters", [])
        if not config or not warning_filters:
            warnings.warn(
                (
                    "`sklearn.utils.parallel.delayed` should be used with"
                    " `sklearn.utils.parallel.Parallel` to make it possible to"
                    " propagate the scikit-learn configuration of the current thread to"
                    " the joblib workers."
                ),
                UserWarning,
            )

        with config_context(**config), warnings.catch_warnings():
            warnings.filters = warning_filters
            return self.function(*args, **kwargs)


def _get_threadpool_controller():
    """Return the global threadpool controller instance."""
    global _threadpool_controller

    if _threadpool_controller is None:
        _threadpool_controller = ThreadpoolController()

    return _threadpool_controller


def _threadpool_controller_decorator(limits=1, user_api="blas"):
    """Decorator to limit the number of threads used at the function level.

    It should be prefered over `threadpoolctl.ThreadpoolController.wrap` because this
    one only loads the shared libraries when the function is called while the latter
    loads them at import time.
    """

    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            controller = _get_threadpool_controller()
            with controller.limit(limits=limits, user_api=user_api):
                return func(*args, **kwargs)

        return wrapper

    return decorator
