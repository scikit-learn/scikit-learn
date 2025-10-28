from polars._dependencies import _HYPOTHESIS_AVAILABLE

if not _HYPOTHESIS_AVAILABLE:
    msg = (
        "polars.testing.parametric requires the 'hypothesis' module\n"
        "Please install it using the command: pip install hypothesis"
    )
    raise ModuleNotFoundError(msg)

from polars.testing.parametric.profiles import load_profile, set_profile
from polars.testing.parametric.strategies import (
    column,
    columns,
    create_list_strategy,
    dataframes,
    dtypes,
    lists,
    series,
)

__all__ = [
    # strategies
    "dataframes",
    "series",
    "column",
    "columns",
    "dtypes",
    "lists",
    "create_list_strategy",
    # profiles
    "load_profile",
    "set_profile",
]
