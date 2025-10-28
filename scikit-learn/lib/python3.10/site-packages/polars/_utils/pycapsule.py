from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING, Any

from polars._utils.construction.dataframe import dataframe_to_pydf
from polars._utils.wrap import wrap_df, wrap_s

with contextlib.suppress(ImportError):
    from polars._plr import PySeries

if TYPE_CHECKING:
    from polars import DataFrame
    from polars._typing import SchemaDefinition, SchemaDict


def is_pycapsule(obj: Any) -> bool:
    """Check if object supports the PyCapsule interface."""
    return hasattr(obj, "__arrow_c_stream__") or hasattr(obj, "__arrow_c_array__")


def pycapsule_to_frame(
    obj: Any,
    *,
    schema: SchemaDefinition | None = None,
    schema_overrides: SchemaDict | None = None,
    rechunk: bool = False,
) -> DataFrame:
    """Convert PyCapsule object to DataFrame."""
    if hasattr(obj, "__arrow_c_array__"):
        # This uses the fact that PySeries.from_arrow_c_array will create a
        # struct-typed Series. Then we unpack that to a DataFrame.
        tmp_col_name = ""
        s = wrap_s(PySeries.from_arrow_c_array(obj))
        df = s.to_frame(tmp_col_name).unnest(tmp_col_name)

    elif hasattr(obj, "__arrow_c_stream__"):
        # This uses the fact that PySeries.from_arrow_c_stream will create a
        # struct-typed Series. Then we unpack that to a DataFrame.
        tmp_col_name = ""
        s = wrap_s(PySeries.from_arrow_c_stream(obj))
        df = s.to_frame(tmp_col_name).unnest(tmp_col_name)
    else:
        msg = f"object does not support PyCapsule interface; found {obj!r} "
        raise TypeError(msg)

    if rechunk:
        df = df.rechunk()
    if schema or schema_overrides:
        df = wrap_df(
            dataframe_to_pydf(df, schema=schema, schema_overrides=schema_overrides)
        )
    return df
