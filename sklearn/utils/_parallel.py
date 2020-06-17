from functools import wraps
from .._config import config_context, get_config
from joblib import delayed as joblib_delayed


class delayed:
    """Decorator used to capture the arguments of a function while passing
    the global configuration options"""

    def __init__(self, func):
        self.func = func

    def __call__(self, *args, **kwargs):
        config = get_config()

        @wraps(self.func)
        def wrapped(*args, **kwargs):
            with config_context(**config):
                return self.func(*args, **kwargs)

        return joblib_delayed(wrapped)(*args, **kwargs)
