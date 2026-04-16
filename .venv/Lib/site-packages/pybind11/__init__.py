from __future__ import annotations

import sys

if sys.version_info < (3, 8):  # noqa: UP036
    msg = "pybind11 does not support Python < 3.8. v2.13 was the last release supporting Python 3.7."
    raise ImportError(msg)


from ._version import __version__, version_info
from .commands import get_cmake_dir, get_include, get_pkgconfig_dir

__all__ = (
    "version_info",
    "__version__",
    "get_include",
    "get_cmake_dir",
    "get_pkgconfig_dir",
)
