from typing import Final

from ..._internal import clone_module

__all__ = clone_module("dask.array", globals())

# These imports may overwrite names from the import * above.
from . import _aliases
from ._aliases import *  # type: ignore[assignment] # noqa: F403
from ._info import __array_namespace_info__  # noqa: F401

__array_api_version__: Final = "2024.12"
del Final

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')
__import__(__package__ + '.fft')

__all__ = sorted(
    set(__all__)
    | set(_aliases.__all__)
    | {"__array_api_version__", "__array_namespace_info__", "linalg", "fft"}
)

def __dir__() -> list[str]:
    return __all__
