from typing import Any

from polars import DataFrame
from polars._dependencies import numpy as np
from polars._typing import IndexOrder
from polars.datatypes import Array, List


def frame_to_numpy(
    df: DataFrame,
    *,
    writable: bool,
    target: str,
    order: IndexOrder = "fortran",
) -> np.ndarray[Any, Any]:
    """Convert a DataFrame to a NumPy array for use with Jax or PyTorch."""
    for nm, tp in df.schema.items():
        if tp == List:
            msg = f"cannot convert List column {nm!r} to {target} (use Array dtype instead)"
            raise TypeError(msg) from None

    if df.width == 1 and df.schema.dtypes()[0] == Array:
        arr = df[df.columns[0]].to_numpy(writable=writable)
    else:
        arr = df.to_numpy(writable=writable, order=order)

    if arr.dtype == object:
        msg = f"cannot convert DataFrame to {target} (mixed type columns result in `object` dtype)\n{df.schema!r}"
        raise TypeError(msg)
    return arr
