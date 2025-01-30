from __future__ import annotations

from typing import Any

from polars._utils.polars_version import get_polars_version

try:
    from polars.polars import __build__
except ImportError:
    __build__ = {}

__build__["version"] = get_polars_version() or "<missing>"


def build_info() -> dict[str, Any]:
    """
    Return detailed Polars build information.

    The dictionary with build information contains the following keys:

    - `"compiler"`
    - `"time"`
    - `"dependencies"`
    - `"features"`
    - `"host"`
    - `"target"`
    - `"git"`
    - `"version"`

    If Polars was compiled without the `build_info` feature flag, only the `"version"`
    key is included.
    """
    return __build__
