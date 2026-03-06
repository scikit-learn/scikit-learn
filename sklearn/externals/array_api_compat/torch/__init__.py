from typing import Final

from .._internal import clone_module

__all__ = clone_module("torch", globals())

# These imports may overwrite names from the import * above.
from . import _aliases
from ._aliases import * # noqa: F403
from ._info import __array_namespace_info__  # noqa: F401

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')
__import__(__package__ + '.fft')

__array_api_version__: Final = '2024.12'

__all__ = sorted(
    set(__all__)
    | set(_aliases.__all__)
    | {"__array_api_version__", "__array_namespace_info__", "linalg", "fft"}
)

def __dir__() -> list[str]:
    return __all__
