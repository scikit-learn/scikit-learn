from functools import update_wrapper
from .._config import config_context, get_config
from joblib import delayed as joblib_delayed


def delayed(f):
    """Decorator used to capture the arguments of a function while passing
    the global configuration options."""
    return joblib_delayed(_FuncWrapper(f))


class _FuncWrapper:
    """"Load the global configuration before calling the function"""
    def __init__(self, func):
        self.func = func
        self.config = get_config()
        update_wrapper(self, self.func)

    def __call__(self, *args, **kwargs):
        with config_context(**self.config):
            return self.func(*args, **kwargs)
