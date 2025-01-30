from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING

import polars._reexport as pl
from polars.dependencies import pyarrow as pa

if TYPE_CHECKING:
    from polars import DataFrame, LazyFrame


def _scan_pyarrow_dataset(
    ds: pa.dataset.Dataset,
    *,
    allow_pyarrow_filter: bool = True,
    batch_size: int | None = None,
) -> LazyFrame:
    """
    Pickle the partially applied function `_scan_pyarrow_dataset_impl`.

    The bytes are then sent to the polars logical plan. It can be deserialized once
    executed and ran.

    Parameters
    ----------
    ds
        pyarrow dataset
    allow_pyarrow_filter
        Allow predicates to be pushed down to pyarrow. This can lead to different
        results if comparisons are done with null values as pyarrow handles this
        different than polars does.
    batch_size
        The maximum row count for scanned pyarrow record batches.
    """
    func = partial(_scan_pyarrow_dataset_impl, ds, batch_size=batch_size)
    return pl.LazyFrame._scan_python_function(
        ds.schema, func, pyarrow=allow_pyarrow_filter
    )


def _scan_pyarrow_dataset_impl(
    ds: pa.dataset.Dataset,
    with_columns: list[str] | None,
    predicate: str | None,
    n_rows: int | None,
    batch_size: int | None,
) -> DataFrame:
    """
    Take the projected columns and materialize an arrow table.

    Parameters
    ----------
    ds
        pyarrow dataset
    with_columns
        Columns that are projected
    predicate
        pyarrow expression that can be evaluated with eval
    n_rows:
        Materialize only n rows from the arrow dataset
    batch_size
        The maximum row count for scanned pyarrow record batches.

    Returns
    -------
    DataFrame
    """
    from polars import from_arrow

    _filter = None

    if predicate:
        from polars._utils.convert import (
            to_py_date,
            to_py_datetime,
            to_py_time,
            to_py_timedelta,
        )
        from polars.datatypes import Date, Datetime, Duration

        _filter = eval(
            predicate,
            {
                "pa": pa,
                "Date": Date,
                "Datetime": Datetime,
                "Duration": Duration,
                "to_py_date": to_py_date,
                "to_py_datetime": to_py_datetime,
                "to_py_time": to_py_time,
                "to_py_timedelta": to_py_timedelta,
            },
        )

    common_params = {"columns": with_columns, "filter": _filter}
    if batch_size is not None:
        common_params["batch_size"] = batch_size

    if n_rows:
        return from_arrow(ds.head(n_rows, **common_params))  # type: ignore[return-value]

    return from_arrow(ds.to_table(**common_params))  # type: ignore[return-value]
