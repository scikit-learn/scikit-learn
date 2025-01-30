from __future__ import annotations

from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from polars import DataFrame


def _execute_from_rust(
    function: Any, with_columns: list[str] | None, *args: Any
) -> DataFrame:
    """
    Deserialize and execute the given function for the projected columns.

    Called from polars-lazy. Polars-lazy provides the bytes of the pickled function and
    the projected columns.

    Parameters
    ----------
    function
        function object
    with_columns
        Columns that are projected
    *args
        Additional function arguments.
    """
    return function(with_columns, *args)
