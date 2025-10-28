from __future__ import annotations

from collections.abc import Sequence
from functools import lru_cache
from typing import TYPE_CHECKING, Any, Callable, get_type_hints

from polars._dependencies import _check_for_pydantic, pydantic

if TYPE_CHECKING:
    import pandas as pd

PANDAS_SIMPLE_NUMPY_DTYPES = {
    "int64",
    "int32",
    "int16",
    "int8",
    "uint64",
    "uint32",
    "uint16",
    "uint8",
    "float64",
    "float32",
    "datetime64[ms]",
    "datetime64[us]",
    "datetime64[ns]",
    "timedelta64[ms]",
    "timedelta64[us]",
    "timedelta64[ns]",
    "bool",
}


def _get_annotations(obj: type) -> dict[str, Any]:
    return getattr(obj, "__annotations__", {})


def try_get_type_hints(obj: type) -> dict[str, Any]:
    try:
        # often the same as obj.__annotations__, but handles forward references
        # encoded as string literals, adds Optional[t] if a default value equal
        # to None is set and recursively replaces 'Annotated[T, ...]' with 'T'.
        return get_type_hints(obj)
    except TypeError:
        # fallback on edge-cases (eg: InitVar inference on python 3.10).
        return _get_annotations(obj)


@lru_cache(64)
def is_namedtuple(cls: Any, *, annotated: bool = False) -> bool:
    """Check if given class derives from NamedTuple."""
    if all(hasattr(cls, attr) for attr in ("_fields", "_field_defaults", "_replace")):
        if not isinstance(cls._fields, property):
            if not annotated or len(cls.__annotations__) == len(cls._fields):
                return all(isinstance(fld, str) for fld in cls._fields)
    return False


def is_pydantic_model(value: Any) -> bool:
    """Check if value derives from a pydantic.BaseModel."""
    return _check_for_pydantic(value) and isinstance(value, pydantic.BaseModel)


def is_sqlalchemy_row(value: Any) -> bool:
    """Check if value is an instance of a SQLAlchemy sequence or mapping object."""
    return getattr(value, "__module__", "").startswith("sqlalchemy.") and isinstance(
        value, Sequence
    )


def get_first_non_none(values: Sequence[Any | None]) -> Any:
    """
    Return the first value from a sequence that isn't None.

    If sequence doesn't contain non-None values, return None.
    """
    if values is not None:
        return next((v for v in values if v is not None), None)


def nt_unpack(obj: Any) -> Any:
    """Recursively unpack a nested NamedTuple."""
    if isinstance(obj, dict):
        return {key: nt_unpack(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [nt_unpack(value) for value in obj]
    elif is_namedtuple(obj.__class__):
        return {key: nt_unpack(value) for key, value in obj._asdict().items()}
    elif isinstance(obj, tuple):
        return tuple(nt_unpack(value) for value in obj)
    else:
        return obj


def contains_nested(value: Any, is_nested: Callable[[Any], bool]) -> bool:
    """Determine if value contains (or is) nested structured data."""
    if is_nested(value):
        return True
    elif isinstance(value, dict):
        return any(contains_nested(v, is_nested) for v in value.values())
    elif isinstance(value, (list, tuple)):
        return any(contains_nested(v, is_nested) for v in value)
    return False


def is_simple_numpy_backed_pandas_series(
    series: pd.Series[Any] | pd.Index[Any] | pd.DatetimeIndex,
) -> bool:
    if len(series.shape) > 1:
        # Pandas Series is actually a Pandas DataFrame when the original DataFrame
        # contains duplicated columns and a duplicated column is requested with df["a"].
        msg = f"duplicate column names found: {series.columns.tolist()!s}"  # type: ignore[union-attr]
        raise ValueError(msg)
    return (str(series.dtype) in PANDAS_SIMPLE_NUMPY_DTYPES) or (
        series.dtype == "object"
        and not series.hasnans
        and not series.empty
        and isinstance(next(iter(series)), str)
    )
