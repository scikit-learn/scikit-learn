from __future__ import annotations

from typing import TYPE_CHECKING, Any

from polars.datatypes.classes import (
    Array,
    DataType,
    DataTypeClass,
    Date,
    Datetime,
    Decimal,
    Duration,
    Float32,
    Float64,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    List,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
)

if TYPE_CHECKING:
    import sys
    from collections.abc import Iterable

    from polars._typing import (
        PolarsDataType,
        PolarsIntegerType,
        PolarsTemporalType,
    )

    if sys.version_info >= (3, 11):
        from typing import Self
    else:
        from typing_extensions import Self


class DataTypeGroup(frozenset):  # type: ignore[type-arg]
    """Group of data types."""

    _match_base_type: bool

    def __new__(
        cls, items: Iterable[DataType | DataTypeClass], *, match_base_type: bool = True
    ) -> Self:
        """
        Construct a DataTypeGroup.

        Parameters
        ----------
        items :
            iterable of data types
        match_base_type:
            match the base type
        """
        for it in items:
            if not isinstance(it, (DataType, DataTypeClass)):
                msg = f"DataTypeGroup items must be dtypes; found {type(it).__name__!r}"
                raise TypeError(msg)
        dtype_group = super().__new__(cls, items)
        dtype_group._match_base_type = match_base_type
        return dtype_group

    def __contains__(self, item: Any) -> bool:
        if self._match_base_type and isinstance(item, (DataType, DataTypeClass)):
            item = item.base_type()
        return super().__contains__(item)


SIGNED_INTEGER_DTYPES: frozenset[PolarsIntegerType] = DataTypeGroup(
    [
        Int8,
        Int16,
        Int32,
        Int64,
        Int128,
    ]
)
UNSIGNED_INTEGER_DTYPES: frozenset[PolarsIntegerType] = DataTypeGroup(
    [
        UInt8,
        UInt16,
        UInt32,
        UInt64,
    ]
)
INTEGER_DTYPES: frozenset[PolarsIntegerType] = (
    SIGNED_INTEGER_DTYPES | UNSIGNED_INTEGER_DTYPES
)
FLOAT_DTYPES: frozenset[PolarsDataType] = DataTypeGroup([Float32, Float64])
NUMERIC_DTYPES: frozenset[PolarsDataType] = DataTypeGroup(
    FLOAT_DTYPES | INTEGER_DTYPES | frozenset([Decimal])
)

DATETIME_DTYPES: frozenset[PolarsDataType] = DataTypeGroup(
    [
        Datetime,
        Datetime("ms"),
        Datetime("us"),
        Datetime("ns"),
        Datetime("ms", "*"),
        Datetime("us", "*"),
        Datetime("ns", "*"),
    ]
)
DURATION_DTYPES: frozenset[PolarsDataType] = DataTypeGroup(
    [
        Duration,
        Duration("ms"),
        Duration("us"),
        Duration("ns"),
    ]
)
TEMPORAL_DTYPES: frozenset[PolarsTemporalType] = DataTypeGroup(
    frozenset([Date, Time]) | DATETIME_DTYPES | DURATION_DTYPES
)

NESTED_DTYPES: frozenset[PolarsDataType] = DataTypeGroup([List, Struct, Array])
