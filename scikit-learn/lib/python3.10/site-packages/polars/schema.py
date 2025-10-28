from __future__ import annotations

import contextlib
import sys
from collections import OrderedDict
from collections.abc import Mapping
from typing import TYPE_CHECKING, Literal, Union, overload

from polars._typing import PythonDataType
from polars._utils.unstable import unstable
from polars.datatypes import DataType, DataTypeClass, is_polars_dtype
from polars.datatypes._parse import parse_into_dtype
from polars.exceptions import DuplicateError
from polars.interchange.protocol import CompatLevel

with contextlib.suppress(ImportError):  # Module not available when building docs
    from polars._plr import (
        init_polars_schema_from_arrow_c_schema,
        polars_schema_field_from_arrow_c_schema,
        polars_schema_to_pycapsule,
    )

if TYPE_CHECKING:
    from collections.abc import Iterable

    from polars import DataFrame, LazyFrame
    from polars._typing import ArrowSchemaExportable

    if sys.version_info >= (3, 10):
        from typing import TypeAlias
    else:
        from typing_extensions import TypeAlias

if sys.version_info >= (3, 10):

    def _required_init_args(tp: DataTypeClass) -> bool:
        # note: this check is ~20% faster than the check for a
        # custom "__init__", below, but is not available on py39
        return bool(tp.__annotations__)
else:

    def _required_init_args(tp: DataTypeClass) -> bool:
        # indicates override of the default __init__
        # (eg: this type requires specific args)
        return "__init__" in tp.__dict__


BaseSchema = OrderedDict[str, DataType]
SchemaInitDataType: TypeAlias = Union[DataType, DataTypeClass, PythonDataType]

__all__ = ["Schema"]


def _check_dtype(tp: DataType | DataTypeClass) -> DataType:
    if not isinstance(tp, DataType):
        # note: if nested/decimal, or has signature params, this implies required args
        if tp.is_nested() or tp.is_decimal() or _required_init_args(tp):
            msg = f"dtypes must be fully-specified, got: {tp!r}"
            raise TypeError(msg)
        tp = tp()
    return tp  # type: ignore[return-value]


class Schema(BaseSchema):
    """
    Ordered mapping of column names to their data type.

    Parameters
    ----------
    schema
        The schema definition given by column names and their associated
        Polars data type. Accepts a mapping, or an iterable of tuples, or any
        object implementing the  `__arrow_c_schema__` PyCapsule interface
        (e.g. pyarrow schemas).

    Examples
    --------
    Define a schema by passing instantiated data types.

    >>> schema = pl.Schema(
    ...     {
    ...         "foo": pl.String(),
    ...         "bar": pl.Duration("us"),
    ...         "baz": pl.Array(pl.Int8, 4),
    ...     }
    ... )
    >>> schema
    Schema({'foo': String, 'bar': Duration(time_unit='us'), 'baz': Array(Int8, shape=(4,))})

    Access the data type associated with a specific column name.

    >>> schema["baz"]
    Array(Int8, shape=(4,))

    Access various schema properties using the `names`, `dtypes`, and `len` methods.

    >>> schema.names()
    ['foo', 'bar', 'baz']
    >>> schema.dtypes()
    [String, Duration(time_unit='us'), Array(Int8, shape=(4,))]
    >>> schema.len()
    3

    Import a pyarrow schema.

    >>> import pyarrow as pa
    >>> pl.Schema(pa.schema([pa.field("x", pa.int32())]))
    Schema({'x': Int32})

    Export a schema to pyarrow.

    >>> pa.schema(pl.Schema({"x": pl.Int32}))
    x: int32
    """  # noqa: W505

    def __init__(
        self,
        schema: (
            Mapping[str, SchemaInitDataType]
            | Iterable[tuple[str, SchemaInitDataType] | ArrowSchemaExportable]
            | ArrowSchemaExportable
            | None
        ) = None,
        *,
        check_dtypes: bool = True,
    ) -> None:
        if hasattr(schema, "__arrow_c_schema__") and not isinstance(schema, Schema):
            init_polars_schema_from_arrow_c_schema(self, schema)
            return

        input = schema.items() if isinstance(schema, Mapping) else (schema or ())
        for v in input:
            name, tp = (
                polars_schema_field_from_arrow_c_schema(v)
                if hasattr(v, "__arrow_c_schema__") and not isinstance(v, DataType)
                else v
            )

            if name in self:
                msg = f"iterable passed to pl.Schema contained duplicate name '{name}'"
                raise DuplicateError(msg)

            if not check_dtypes:
                super().__setitem__(name, tp)  # type: ignore[assignment]
            elif is_polars_dtype(tp):
                super().__setitem__(name, _check_dtype(tp))
            else:
                self[name] = tp

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Mapping):
            return False
        if len(self) != len(other):
            return False
        for (nm1, tp1), (nm2, tp2) in zip(self.items(), other.items()):
            if nm1 != nm2 or not tp1.is_(tp2):
                return False
        return True

    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def __setitem__(
        self, name: str, dtype: DataType | DataTypeClass | PythonDataType
    ) -> None:
        dtype = _check_dtype(parse_into_dtype(dtype))
        super().__setitem__(name, dtype)

    @unstable()
    def __arrow_c_schema__(self) -> object:
        """
        Export a Schema via the Arrow PyCapsule Interface.

        https://arrow.apache.org/docs/dev/format/CDataInterface/PyCapsuleInterface.html
        """
        return polars_schema_to_pycapsule(self, CompatLevel.newest()._version)

    def names(self) -> list[str]:
        """
        Get the column names of the schema.

        Examples
        --------
        >>> s = pl.Schema({"x": pl.Float64(), "y": pl.Datetime(time_zone="UTC")})
        >>> s.names()
        ['x', 'y']
        """
        return list(self.keys())

    def dtypes(self) -> list[DataType]:
        """
        Get the data types of the schema.

        Examples
        --------
        >>> s = pl.Schema({"x": pl.UInt8(), "y": pl.List(pl.UInt8)})
        >>> s.dtypes()
        [UInt8, List(UInt8)]
        """
        return list(self.values())

    @overload
    def to_frame(self, *, eager: Literal[False]) -> LazyFrame: ...

    @overload
    def to_frame(self, *, eager: Literal[True] = ...) -> DataFrame: ...

    def to_frame(self, *, eager: bool = True) -> DataFrame | LazyFrame:
        """
        Create an empty DataFrame (or LazyFrame) from this Schema.

        Parameters
        ----------
        eager
            If True, create a DataFrame; otherwise, create a LazyFrame.

        Examples
        --------
        >>> s = pl.Schema({"x": pl.Int32(), "y": pl.String()})
        >>> s.to_frame()
        shape: (0, 2)
        ┌─────┬─────┐
        │ x   ┆ y   │
        │ --- ┆ --- │
        │ i32 ┆ str │
        ╞═════╪═════╡
        └─────┴─────┘
        >>> s.to_frame(eager=False)  # doctest: +IGNORE_RESULT
        <LazyFrame at 0x11BC0AD80>
        """
        from polars import DataFrame, LazyFrame

        return DataFrame(schema=self) if eager else LazyFrame(schema=self)

    def len(self) -> int:
        """
        Get the number of schema entries.

        Examples
        --------
        >>> s = pl.Schema({"x": pl.Int32(), "y": pl.List(pl.String)})
        >>> s.len()
        2
        >>> len(s)
        2
        """
        return len(self)

    def to_python(self) -> dict[str, type]:
        """
        Return a dictionary of column names and Python types.

        Examples
        --------
        >>> s = pl.Schema(
        ...     {
        ...         "x": pl.Int8(),
        ...         "y": pl.String(),
        ...         "z": pl.Duration("us"),
        ...     }
        ... )
        >>> s.to_python()
        {'x': <class 'int'>, 'y':  <class 'str'>, 'z': <class 'datetime.timedelta'>}
        """
        return {name: tp.to_python() for name, tp in self.items()}
