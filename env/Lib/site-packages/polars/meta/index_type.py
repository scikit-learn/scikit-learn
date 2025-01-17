from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars.polars as plr

if TYPE_CHECKING:
    from polars.datatypes import DataType


def get_index_type() -> DataType:
    """
    Return the data type used for Polars indexing.

    Returns
    -------
    DataType
        :class:`UInt32` in regular Polars, :class:`UInt64` in bigidx Polars.

    Examples
    --------
    >>> pl.get_index_type()
    UInt32
    """
    return plr.get_index_type()
