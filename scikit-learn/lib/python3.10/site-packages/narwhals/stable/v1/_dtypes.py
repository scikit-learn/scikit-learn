from __future__ import annotations

from typing import TYPE_CHECKING

from narwhals._utils import inherit_doc
from narwhals.dtypes import (
    Array,
    Binary,
    Boolean,
    Categorical,
    Date,
    Datetime as NwDatetime,
    Decimal,
    DType,
    Duration as NwDuration,
    Enum as NwEnum,
    Field,
    Float32,
    Float64,
    FloatType,
    Int8,
    Int16,
    Int32,
    Int64,
    Int128,
    IntegerType,
    List,
    NestedType,
    NumericType,
    Object,
    SignedIntegerType,
    String,
    Struct,
    Time,
    UInt8,
    UInt16,
    UInt32,
    UInt64,
    UInt128,
    Unknown,
    UnsignedIntegerType,
)

if TYPE_CHECKING:
    from datetime import timezone

    from narwhals.typing import TimeUnit


class Datetime(NwDatetime):
    __slots__ = NwDatetime.__slots__

    @inherit_doc(NwDatetime)
    def __init__(
        self, time_unit: TimeUnit = "us", time_zone: str | timezone | None = None
    ) -> None:
        super().__init__(time_unit, time_zone)

    def __hash__(self) -> int:
        return hash(self.__class__)


class Duration(NwDuration):
    __slots__ = NwDuration.__slots__

    @inherit_doc(NwDuration)
    def __init__(self, time_unit: TimeUnit = "us") -> None:
        super().__init__(time_unit)

    def __hash__(self) -> int:
        return hash(self.__class__)


class Enum(NwEnum):
    """A fixed categorical encoding of a unique set of strings.

    Polars has an Enum data type, while pandas and PyArrow do not.

    Examples:
       >>> import polars as pl
       >>> import narwhals.stable.v1 as nw
       >>> data = ["beluga", "narwhal", "orca"]
       >>> s_native = pl.Series(data, dtype=pl.Enum(data))
       >>> nw.from_native(s_native, series_only=True).dtype
       Enum
    """

    __slots__ = ()

    def __init__(self) -> None:
        super(NwEnum, self).__init__()

    def __eq__(self, other: DType | type[DType]) -> bool:  # type: ignore[override]
        if type(other) is type:
            return other in {type(self), NwEnum}
        return isinstance(other, type(self))

    def __hash__(self) -> int:  # pragma: no cover
        return super(NwEnum, self).__hash__()

    def __repr__(self) -> str:  # pragma: no cover
        return super(NwEnum, self).__repr__()


__all__ = [
    "Array",
    "Binary",
    "Boolean",
    "Categorical",
    "DType",
    "Date",
    "Datetime",
    "Decimal",
    "Duration",
    "Enum",
    "Field",
    "Float32",
    "Float64",
    "FloatType",
    "Int8",
    "Int16",
    "Int32",
    "Int64",
    "Int128",
    "IntegerType",
    "List",
    "NestedType",
    "NumericType",
    "Object",
    "SignedIntegerType",
    "String",
    "Struct",
    "Time",
    "UInt8",
    "UInt16",
    "UInt32",
    "UInt64",
    "UInt128",
    "Unknown",
    "UnsignedIntegerType",
]
