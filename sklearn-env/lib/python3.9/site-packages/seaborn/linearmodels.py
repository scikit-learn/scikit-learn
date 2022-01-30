import warnings
from .regression import *  # noqa

msg = (
    "The `linearmodels` module has been renamed `regression`."
)
warnings.warn(msg)
