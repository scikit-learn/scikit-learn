# tifffile/__init__.py

from .tifffile import *
from .tifffile import __all__, __doc__, __version__, main


def _set_module() -> None:
    """Set __module__ attribute for all public objects."""
    globs = globals()
    for item in __all__:
        obj = globs[item]
        if hasattr(obj, '__module__'):
            obj.__module__ = 'tifffile'
    main.__module__ = 'tifffile'


_set_module()
