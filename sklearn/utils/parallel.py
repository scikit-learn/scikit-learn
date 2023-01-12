"""Module that customize joblib tools for scikit-learn usage."""

import functools
import warnings
from functools import update_wrapper

import joblib

from .._config import config_context, get_config


def _with_config(delayed_func, config):
    """Helper function that intends to attach a config to a delayed function."""
    if hasattr(delayed_func, "with_config"):
        return delayed_func.with_config(config)
    else:
        warnings.warn(
            "You are using `sklearn.utils.fixes.Parallel` that intend to attach a "
            "configuration to a delayed function. However, the function used for "
            "delaying the function does not expose `with_config`. Use "
            "`sklearn.utils.fixes.delayed` to correctly propagate the scikit-learn "
            "configuration to the joblib workers.",
            UserWarning,
        )
        return delayed_func


class Parallel(joblib.Parallel):
    # A tweaked `Parallel` subclass that attaches the configuration of the current
    # thread to each task to be run in parallel.

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
        iterable_with_config = (
            (_with_config(delayed_func, config), args, kwargs)
            for delayed_func, args, kwargs in iterable
        )
        return super().__call__(iterable_with_config)


# remove when https://github.com/joblib/joblib/issues/1071 is fixed
def delayed(function):
    """Decorator used to capture the arguments of a function.

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
    """ "Load the global configuration before calling the function."""

    def __init__(self, function):
        self.function = function
        update_wrapper(self, self.function)

    def with_config(self, config):
        self.config = config
        return self

    def __call__(self, *args, **kwargs):
        config = getattr(self, "config", None)
        if config is None:
            warnings.warn(
                "`sklearn.utils.fixes.delayed` should be used with "
                "`sklearn.utils.fixes.Parallel` to make it possible to propagate "
                "the scikit-learn configuration of the current thread to the "
                "joblib workers.",
                UserWarning,
            )
            config = {}
        with config_context(**config):
            return self.function(*args, **kwargs)
