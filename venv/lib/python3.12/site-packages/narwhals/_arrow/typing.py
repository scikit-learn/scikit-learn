from __future__ import annotations  # pragma: no cover

from typing import (
    TYPE_CHECKING,  # pragma: no cover
    Any,  # pragma: no cover
    TypeVar,  # pragma: no cover
)

if TYPE_CHECKING:
    import sys
    from typing import Generic, Literal

    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

    import pyarrow as pa
    from pyarrow.__lib_pxi.table import (
        AggregateOptions,  # noqa: F401
        Aggregation,  # noqa: F401
    )
    from pyarrow._stubs_typing import (  # pyright: ignore[reportMissingModuleSource]  # pyright: ignore[reportMissingModuleSource]  # pyright: ignore[reportMissingModuleSource]
        Indices,  # noqa: F401
        Mask,  # noqa: F401
        Order,  # noqa: F401
    )

    from narwhals._arrow.expr import ArrowExpr
    from narwhals._arrow.series import ArrowSeries

    IntoArrowExpr: TypeAlias = "ArrowExpr | ArrowSeries"
    TieBreaker: TypeAlias = Literal["min", "max", "first", "dense"]
    NullPlacement: TypeAlias = Literal["at_start", "at_end"]
    NativeIntervalUnit: TypeAlias = Literal[
        "year",
        "quarter",
        "month",
        "week",
        "day",
        "hour",
        "minute",
        "second",
        "millisecond",
        "microsecond",
        "nanosecond",
    ]

    ChunkedArrayAny: TypeAlias = pa.ChunkedArray[Any]
    ArrayAny: TypeAlias = pa.Array[Any]
    ArrayOrChunkedArray: TypeAlias = "ArrayAny | ChunkedArrayAny"
    ScalarAny: TypeAlias = pa.Scalar[Any]
    ArrayOrScalar: TypeAlias = "ArrayOrChunkedArray | ScalarAny"
    ArrayOrScalarT1 = TypeVar("ArrayOrScalarT1", ArrayAny, ChunkedArrayAny, ScalarAny)
    ArrayOrScalarT2 = TypeVar("ArrayOrScalarT2", ArrayAny, ChunkedArrayAny, ScalarAny)
    _AsPyType = TypeVar("_AsPyType")

    class _BasicDataType(pa.DataType, Generic[_AsPyType]): ...


Incomplete: TypeAlias = Any  # pragma: no cover
"""
Marker for working code that fails on the stubs.

Common issues:
- Annotated for `Array`, but not `ChunkedArray`
- Relies on typing information that the stubs don't provide statically
- Missing attributes
- Incorrect return types
- Inconsistent use of generic/concrete types
- `_clone_signature` used on signatures that are not identical
"""
