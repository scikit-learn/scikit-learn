from typing import Final
from cupy import * # noqa: F403

# from cupy import * doesn't overwrite these builtin names
from cupy import abs, max, min, round # noqa: F401

# These imports may overwrite names from the import * above.
from ._aliases import * # noqa: F403
from ._info import __array_namespace_info__  # noqa: F401

# See the comment in the numpy __init__.py
__import__(__package__ + '.linalg')
__import__(__package__ + '.fft')

__array_api_version__: Final = '2024.12'

__all__ = sorted(
    {name for name in globals() if not name.startswith("__")}
    - {"Final", "_aliases", "_info", "_typing"}
    | {"__array_api_version__", "__array_namespace_info__", "linalg", "fft"}
)

def __dir__() -> list[str]:
    return __all__
