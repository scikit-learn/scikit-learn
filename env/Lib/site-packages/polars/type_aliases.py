"""
Deprecated module - do not use.

Used to contain private type aliases. These are now in the `polars._typing` module.
"""

from typing import Any

import polars._typing as plt
from polars._utils.deprecation import issue_deprecation_warning


def __getattr__(name: str) -> Any:
    if name in dir(plt):
        issue_deprecation_warning(
            "The `polars.type_aliases` module is deprecated."
            " The type aliases have moved to the `polars._typing` module to explicitly mark them as private."
            " Please define your own type aliases, or temporarily import from the `polars._typing` module."
            " A public `polars.typing` module will be added in the future.",
            version="1.0.0",
        )
        return getattr(plt, name)

    msg = f"module {__name__!r} has no attribute {name!r}"
    raise AttributeError(msg)
