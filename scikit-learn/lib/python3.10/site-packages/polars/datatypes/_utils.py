"""Utility functions for handling and processing of datatypes."""

from polars._typing import PolarsDataType
from polars.datatypes.classes import Array, List, Struct


def dtype_to_init_repr(dtype: PolarsDataType, prefix: str = "pl.") -> str:
    """Convert a Polars dtype to a prefixed string representation."""
    if isinstance(dtype, List):
        init_repr = _dtype_to_init_repr_list(dtype, prefix)
    elif isinstance(dtype, Array):
        init_repr = _dtype_to_init_repr_array(dtype, prefix)
    elif isinstance(dtype, Struct):
        init_repr = _dtype_to_init_repr_struct(dtype, prefix)
    else:
        init_repr = f"{prefix}{dtype!r}"
    return init_repr


def _dtype_to_init_repr_list(dtype: List, prefix: str) -> str:
    class_name = dtype.__class__.__name__
    if dtype.inner is not None:
        inner_repr = dtype_to_init_repr(dtype.inner, prefix)
    else:
        inner_repr = ""
    init_repr = f"{prefix}{class_name}({inner_repr})"
    return init_repr


def _dtype_to_init_repr_array(dtype: Array, prefix: str) -> str:
    class_name = dtype.__class__.__name__
    if dtype.inner is not None:
        inner_repr = dtype_to_init_repr(dtype.inner, prefix)
    else:
        inner_repr = ""
    init_repr = f"{prefix}{class_name}({inner_repr}, shape={dtype.shape})"
    return init_repr


def _dtype_to_init_repr_struct(dtype: Struct, prefix: str) -> str:
    class_name = dtype.__class__.__name__
    inner_list = [
        f"{field_name!r}: {dtype_to_init_repr(inner_dtype, prefix)}"
        for field_name, inner_dtype in dict(dtype).items()
    ]
    inner_repr = "{" + ", ".join(inner_list) + "}"
    init_repr = f"{prefix}{class_name}({inner_repr})"
    return init_repr
