import warnings
msg = (
    "The `linearmodels` module has been renamed `regression`."
)
warnings.warn(msg)
from .regression import *
