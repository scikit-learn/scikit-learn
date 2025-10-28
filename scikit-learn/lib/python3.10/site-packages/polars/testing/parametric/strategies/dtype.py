from __future__ import annotations

from typing import TYPE_CHECKING

import hypothesis.strategies as st
from hypothesis.errors import InvalidArgument

from polars.datatypes import (
    Array,
    Binary,
    Boolean,
    Categorical,
    DataType,
    Date,
    Datetime,
    Decimal,
    Duration,
    Enum,
    Field,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    List,
    Null,
    String,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
)

if TYPE_CHECKING:
    from collections.abc import Collection, Sequence

    from hypothesis.strategies import DrawFn, SearchStrategy

    from polars._typing import PolarsDataType, TimeUnit
    from polars.datatypes import DataTypeClass


# Supported data type classes which do not take any arguments
_SIMPLE_DTYPES: list[DataTypeClass] = [
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    Float64,
    Float32,
    Boolean,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    String,
    Binary,
    Date,
    Time,
    Null,
    # TODO: Enable Object types by default when various issues are solved.
    # Object,
]
# Supported data type classes with arguments
_COMPLEX_DTYPES: list[DataTypeClass] = [
    Datetime,
    Duration,
    Categorical,
    Decimal,
    Enum,
]
# Supported data type classes that contain other data types
_NESTED_DTYPES: list[DataTypeClass] = [
    # TODO: Enable nested types by default when various issues are solved.
    # List,
    # Array,
    Struct,
]
# Supported data type classes that do not contain other data types
_FLAT_DTYPES = _SIMPLE_DTYPES + _COMPLEX_DTYPES

_DEFAULT_ARRAY_WIDTH_LIMIT = 3
_DEFAULT_STRUCT_FIELDS_LIMIT = 3
_DEFAULT_ENUM_CATEGORIES_LIMIT = 3


def dtypes(
    *,
    allowed_dtypes: Collection[PolarsDataType] | None = None,
    excluded_dtypes: Sequence[PolarsDataType] | None = None,
    allow_time_zones: bool = True,
    nesting_level: int = 3,
) -> SearchStrategy[DataType]:
    """
    Create a strategy for generating Polars :class:`DataType` objects.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    Parameters
    ----------
    allowed_dtypes
        Data types the strategy will pick from. If set to `None` (default),
        all supported data types are included.
    excluded_dtypes
        Data types the strategy will *not* pick from. This takes priority over
        data types specified in `allowed_dtypes`.
    allow_time_zones
        Allow generating `Datetime` data types with a time zone.
    nesting_level
        The complexity of nested data types. If set to 0, nested data types are
        disabled.
    """
    flat_dtypes, nested_dtypes, excluded_dtypes = _parse_dtype_restrictions(
        allowed_dtypes, excluded_dtypes
    )

    if nesting_level > 0 and nested_dtypes:
        if not flat_dtypes:
            return _nested_dtypes(
                inner=st.just(Null()),
                allowed_dtypes=nested_dtypes,
                excluded_dtypes=excluded_dtypes,
                allow_time_zones=allow_time_zones,
            )
        return st.recursive(
            base=_flat_dtypes(
                allowed_dtypes=flat_dtypes,
                excluded_dtypes=excluded_dtypes,
                allow_time_zones=allow_time_zones,
            ),
            extend=lambda s: _nested_dtypes(
                s,
                allowed_dtypes=nested_dtypes,
                excluded_dtypes=excluded_dtypes,
                allow_time_zones=allow_time_zones,
            ),
            max_leaves=nesting_level,
        )
    else:
        return _flat_dtypes(
            allowed_dtypes=flat_dtypes,
            excluded_dtypes=excluded_dtypes,
            allow_time_zones=allow_time_zones,
        )


def _parse_dtype_restrictions(
    allowed_dtypes: Collection[PolarsDataType] | None = None,
    excluded_dtypes: Sequence[PolarsDataType] | None = None,
) -> tuple[list[PolarsDataType], list[PolarsDataType], list[DataType]]:
    """
    Parse data type restrictions.

    Splits allowed data types into flat and nested data types.
    Filters the allowed data types by excluded data type classes.
    Excluded instantiated data types are returned to be filtered later.
    """
    # Split excluded dtypes into instances and classes
    excluded_dtypes_instance = []
    excluded_dtypes_class = []
    if excluded_dtypes:
        for dt in excluded_dtypes:
            if isinstance(dt, DataType):
                excluded_dtypes_instance.append(dt)
            else:
                excluded_dtypes_class.append(dt)

    # Split allowed dtypes into flat and nested, excluding certain dtype classes
    allowed_dtypes_flat: list[PolarsDataType]
    allowed_dtypes_nested: list[PolarsDataType]
    if allowed_dtypes is None:
        allowed_dtypes_flat = [
            dt for dt in _FLAT_DTYPES if dt not in excluded_dtypes_class
        ]
        allowed_dtypes_nested = [
            dt for dt in _NESTED_DTYPES if dt not in excluded_dtypes_class
        ]
    else:
        allowed_dtypes_flat = []
        allowed_dtypes_nested = []
        for dt in allowed_dtypes:
            if dt in excluded_dtypes_class:
                continue
            elif dt.is_nested():
                allowed_dtypes_nested.append(dt)
            else:
                allowed_dtypes_flat.append(dt)

    return allowed_dtypes_flat, allowed_dtypes_nested, excluded_dtypes_instance


@st.composite
def _flat_dtypes(
    draw: DrawFn,
    allowed_dtypes: Sequence[PolarsDataType] | None = None,
    excluded_dtypes: Sequence[PolarsDataType] | None = None,
    *,
    allow_time_zones: bool = True,
) -> DataType:
    """Create a strategy for generating non-nested Polars :class:`DataType` objects."""
    if allowed_dtypes is None:
        allowed_dtypes = _FLAT_DTYPES
    if excluded_dtypes is None:
        excluded_dtypes = []

    dtype = draw(st.sampled_from(allowed_dtypes))
    return draw(
        _instantiate_flat_dtype(dtype, allow_time_zones=allow_time_zones).filter(
            lambda x: x not in excluded_dtypes
        )
    )


@st.composite
def _instantiate_flat_dtype(
    draw: DrawFn, dtype: PolarsDataType, *, allow_time_zones: bool = True
) -> DataType:
    """Take a flat data type and instantiate it."""
    if isinstance(dtype, DataType):
        return dtype
    elif dtype in _SIMPLE_DTYPES:
        return dtype()
    elif dtype == Datetime:
        time_unit = draw(_time_units())
        time_zone = draw(st.none() | _time_zones()) if allow_time_zones else None
        return Datetime(time_unit, time_zone)
    elif dtype == Duration:
        time_unit = draw(_time_units())
        return Duration(time_unit)
    elif dtype == Categorical:
        return Categorical()
    elif dtype == Enum:
        n_categories = draw(
            st.integers(min_value=1, max_value=_DEFAULT_ENUM_CATEGORIES_LIMIT)
        )
        categories = [f"c{i}" for i in range(n_categories)]
        return Enum(categories)
    elif dtype == Decimal:
        precision = draw(st.integers(min_value=1, max_value=38) | st.none())
        scale = draw(st.integers(min_value=0, max_value=precision or 38))
        return Decimal(precision, scale)
    else:
        msg = f"unsupported data type: {dtype}"
        raise InvalidArgument(msg)


@st.composite
def _nested_dtypes(
    draw: DrawFn,
    inner: SearchStrategy[DataType],
    allowed_dtypes: Sequence[PolarsDataType] | None = None,
    excluded_dtypes: Sequence[PolarsDataType] | None = None,
    *,
    allow_time_zones: bool = True,
) -> DataType:
    """Create a strategy for generating nested Polars :class:`DataType` objects."""
    if allowed_dtypes is None:
        allowed_dtypes = _NESTED_DTYPES
    if excluded_dtypes is None:
        excluded_dtypes = []

    dtype = draw(st.sampled_from(allowed_dtypes))
    return draw(
        _instantiate_nested_dtype(
            dtype, inner, allow_time_zones=allow_time_zones
        ).filter(lambda x: x not in excluded_dtypes)
    )


@st.composite
def _instantiate_nested_dtype(
    draw: DrawFn,
    dtype: PolarsDataType,
    inner: SearchStrategy[DataType],
    *,
    allow_time_zones: bool = True,
) -> DataType:
    """Take a nested data type and instantiate it."""

    def instantiate_inner(inner_dtype: PolarsDataType | None) -> DataType:
        if inner_dtype is None:
            return draw(inner)
        elif inner_dtype.is_nested():
            return draw(
                _instantiate_nested_dtype(
                    inner_dtype, inner, allow_time_zones=allow_time_zones
                )
            )
        else:
            return draw(
                _instantiate_flat_dtype(inner_dtype, allow_time_zones=allow_time_zones)
            )

    if dtype == List:
        inner_dtype = instantiate_inner(getattr(dtype, "inner", None))
        return List(inner_dtype)
    elif dtype == Array:
        inner_dtype = instantiate_inner(getattr(dtype, "inner", None))
        size = getattr(
            dtype,
            "size",
            draw(st.integers(min_value=1, max_value=_DEFAULT_ARRAY_WIDTH_LIMIT)),
        )
        return Array(inner_dtype, size)
    elif dtype == Struct:
        if isinstance(dtype, Struct):
            fields = [Field(f.name, instantiate_inner(f.dtype)) for f in dtype.fields]
        else:
            n_fields = draw(
                st.integers(min_value=1, max_value=_DEFAULT_STRUCT_FIELDS_LIMIT)
            )
            fields = [Field(f"f{i}", draw(inner)) for i in range(n_fields)]
        return Struct(fields)
    else:
        msg = f"unsupported data type: {dtype}"
        raise InvalidArgument(msg)


def _time_units() -> SearchStrategy[TimeUnit]:
    """Create a strategy for generating valid units of time."""
    return st.sampled_from(["us", "ns", "ms"])


def _time_zones() -> SearchStrategy[str]:
    """Create a strategy for generating valid time zones."""
    # Not available when building docs, so just import here.
    from polars._plr import _known_timezones

    chrono_known_tz = set(_known_timezones())
    return st.timezone_keys(allow_prefix=False).filter(
        lambda tz: tz not in {"Factory", "localtime"} and tz in chrono_known_tz
    )


@st.composite
def _instantiate_dtype(
    draw: DrawFn,
    dtype: PolarsDataType,
    *,
    allowed_dtypes: Collection[PolarsDataType] | None = None,
    excluded_dtypes: Sequence[PolarsDataType] | None = None,
    nesting_level: int = 3,
    allow_time_zones: bool = True,
) -> DataType:
    """Take a data type and instantiate it."""
    if not dtype.is_nested():
        if isinstance(dtype, DataType):
            return dtype

        if allowed_dtypes is None:
            allowed_dtypes = [dtype]
        else:
            same_dtypes = [dt for dt in allowed_dtypes if dt == dtype]
            allowed_dtypes = same_dtypes if same_dtypes else [dtype]

        return draw(
            _flat_dtypes(
                allowed_dtypes=allowed_dtypes,
                excluded_dtypes=excluded_dtypes,
                allow_time_zones=allow_time_zones,
            )
        )

    def draw_inner(dtype: PolarsDataType | None) -> DataType:
        if dtype is None:
            return draw(
                dtypes(
                    allowed_dtypes=allowed_dtypes,
                    excluded_dtypes=excluded_dtypes,
                    nesting_level=nesting_level - 1,
                    allow_time_zones=allow_time_zones,
                )
            )
        else:
            return draw(
                _instantiate_dtype(
                    dtype,
                    allowed_dtypes=allowed_dtypes,
                    excluded_dtypes=excluded_dtypes,
                    nesting_level=nesting_level - 1,
                    allow_time_zones=allow_time_zones,
                )
            )

    if dtype == List:
        inner = draw_inner(getattr(dtype, "inner", None))
        return List(inner)
    elif dtype == Array:
        inner = draw_inner(getattr(dtype, "inner", None))
        size = getattr(
            dtype,
            "size",
            draw(st.integers(min_value=1, max_value=_DEFAULT_ARRAY_WIDTH_LIMIT)),
        )
        return Array(inner, size)
    elif dtype == Struct:
        if isinstance(dtype, Struct):
            fields = [
                Field(
                    name=f.name,
                    dtype=draw(
                        _instantiate_dtype(
                            f.dtype,
                            allowed_dtypes=allowed_dtypes,
                            excluded_dtypes=excluded_dtypes,
                            nesting_level=nesting_level - 1,
                            allow_time_zones=allow_time_zones,
                        )
                    ),
                )
                for f in dtype.fields
            ]
        else:
            n_fields = draw(
                st.integers(min_value=1, max_value=_DEFAULT_STRUCT_FIELDS_LIMIT)
            )
            inner_strategy = dtypes(
                allowed_dtypes=allowed_dtypes,
                excluded_dtypes=excluded_dtypes,
                nesting_level=nesting_level - 1,
                allow_time_zones=allow_time_zones,
            )
            fields = [Field(f"f{i}", draw(inner_strategy)) for i in range(n_fields)]
        return Struct(fields)
    else:
        msg = f"unsupported data type: {dtype}"
        raise InvalidArgument(msg)
