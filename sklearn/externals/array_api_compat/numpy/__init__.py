# ruff: noqa: PLC0414
from typing import Final

from .._internal import clone_module

# This needs to be loaded explicitly before cloning
import numpy.typing  # noqa: F401

__all__ = clone_module("numpy", globals())

# These imports may overwrite names from the import * above.
from . import _aliases
from ._aliases import *  # type: ignore[assignment,no-redef] # noqa: F403
from ._info import __array_namespace_info__  # noqa: F401

# Don't know why, but we have to do an absolute import to import linalg. If we
# instead do
#
# from . import linalg
#
# It doesn't overwrite np.linalg from above. The import is generated
# dynamically so that the library can be vendored.
__import__(__package__ + ".linalg")

__import__(__package__ + ".fft")

from .linalg import matrix_transpose, vecdot  # type: ignore[no-redef]  # noqa: F401

__array_api_version__: Final = "2024.12"

__all__ = sorted(
    set(__all__) 
    | set(_aliases.__all__) 
    | {"__array_api_version__", "__array_namespace_info__", "linalg", "fft"}
)

def __dir__() -> list[str]:
    return __all__
