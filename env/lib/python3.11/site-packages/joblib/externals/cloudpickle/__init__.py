from . import cloudpickle
from .cloudpickle import *  # noqa

__doc__ = cloudpickle.__doc__

__version__ = "3.0.0"

__all__ = [  # noqa
    "__version__",
    "Pickler",
    "CloudPickler",
    "dumps",
    "loads",
    "dump",
    "load",
    "register_pickle_by_value",
    "unregister_pickle_by_value",
]
