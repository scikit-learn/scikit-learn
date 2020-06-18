from functools import update_wrapper
from .._config import config_context, get_config
from joblib import delayed as joblib_delayed


def delayed(function, *args, **kwargs):
    """Wrapper around joblib.delayed to pass around the global configuration.
    """
    return joblib_delayed(_FuncWrapper(function), *args, **kwargs)


class _FuncWrapper:
    """"Load the global configuration before calling the function."""
    def __init__(self, function):
        self.function = function
        self.config = get_config()
        update_wrapper(self, self.function)

    def __call__(self, *args, **kwargs):
        with config_context(**self.config):
            return self.function(*args, **kwargs)
