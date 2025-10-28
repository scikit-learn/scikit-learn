from __future__ import annotations

import abc
import ast
import contextlib
from _ast import GtE, Lt, LtE
from ast import (
    Attribute,
    BinOp,
    BitAnd,
    BitOr,
    Call,
    Compare,
    Constant,
    Eq,
    Gt,
    Invert,
    List,
    Name,
    UnaryOp,
)
from dataclasses import dataclass
from functools import cache, singledispatch
from typing import TYPE_CHECKING, Any, Callable

import polars._reexport as pl
from polars._utils.convert import to_py_date, to_py_datetime
from polars._utils.logging import eprint
from polars._utils.wrap import wrap_s
from polars.exceptions import ComputeError

if TYPE_CHECKING:
    from collections.abc import Sequence
    from datetime import date, datetime

    import pyiceberg
    import pyiceberg.schema
    from pyiceberg.manifest import DataFile
    from pyiceberg.table import Table
    from pyiceberg.types import IcebergType

    from polars import DataFrame, Series
else:
    from polars._dependencies import pyiceberg

_temporal_conversions: dict[str, Callable[..., datetime | date]] = {
    "to_py_date": to_py_date,
    "to_py_datetime": to_py_datetime,
}

ICEBERG_TIME_TO_NS: int = 1000


def _scan_pyarrow_dataset_impl(
    tbl: Table,
    with_columns: list[str] | None = None,
    iceberg_table_filter: Any | None = None,
    n_rows: int | None = None,
    snapshot_id: int | None = None,
    **kwargs: Any,
) -> DataFrame | Series:
    """
    Take the projected columns and materialize an arrow table.

    Parameters
    ----------
    tbl
        pyarrow dataset
    with_columns
        Columns that are projected
    iceberg_table_filter
        PyIceberg filter expression
    n_rows:
        Materialize only n rows from the arrow dataset.
    snapshot_id:
        The snapshot ID to scan from.
    batch_size
        The maximum row count for scanned pyarrow record batches.
    kwargs:
        For backward compatibility

    Returns
    -------
    DataFrame
    """
    from polars import from_arrow

    scan = tbl.scan(limit=n_rows, snapshot_id=snapshot_id)

    if with_columns is not None:
        scan = scan.select(*with_columns)

    if iceberg_table_filter is not None:
        scan = scan.filter(iceberg_table_filter)

    return from_arrow(scan.to_arrow())


def try_convert_pyarrow_predicate(pyarrow_predicate: str) -> Any | None:
    with contextlib.suppress(Exception):
        expr_ast = _to_ast(pyarrow_predicate)
        return _convert_predicate(expr_ast)

    return None


def _to_ast(expr: str) -> ast.expr:
    """
    Converts a Python string to an AST.

    This will take the Python Arrow expression (as a string), and it will
    be converted into a Python AST that can be traversed to convert it to a PyIceberg
    expression.

    The reason to convert it to an AST is because the PyArrow expression
    itself doesn't have any methods/properties to traverse the expression.
    We need this to convert it into a PyIceberg expression.

    Parameters
    ----------
    expr
        The string expression

    Returns
    -------
    The AST representing the Arrow expression
    """
    return ast.parse(expr, mode="eval").body


@singledispatch
def _convert_predicate(a: Any) -> Any:
    """Walks the AST to convert the PyArrow expression to a PyIceberg expression."""
    msg = f"Unexpected symbol: {a}"
    raise ValueError(msg)


@_convert_predicate.register(Constant)
def _(a: Constant) -> Any:
    return a.value


@_convert_predicate.register(Name)
def _(a: Name) -> Any:
    return a.id


@_convert_predicate.register(UnaryOp)
def _(a: UnaryOp) -> Any:
    if isinstance(a.op, Invert):
        return pyiceberg.expressions.Not(_convert_predicate(a.operand))
    else:
        msg = f"Unexpected UnaryOp: {a}"
        raise TypeError(msg)


@_convert_predicate.register(Call)
def _(a: Call) -> Any:
    args = [_convert_predicate(arg) for arg in a.args]
    f = _convert_predicate(a.func)
    if f == "field":
        return args
    elif f == "scalar":
        return args[0]
    elif f in _temporal_conversions:
        # convert from polars-native i64 to ISO8601 string
        return _temporal_conversions[f](*args).isoformat()
    else:
        ref = _convert_predicate(a.func.value)[0]  # type: ignore[attr-defined]
        if f == "isin":
            return pyiceberg.expressions.In(ref, args[0])
        elif f == "is_null":
            return pyiceberg.expressions.IsNull(ref)
        elif f == "is_nan":
            return pyiceberg.expressions.IsNaN(ref)

    msg = f"Unknown call: {f!r}"
    raise ValueError(msg)


@_convert_predicate.register(Attribute)
def _(a: Attribute) -> Any:
    return a.attr


@_convert_predicate.register(BinOp)
def _(a: BinOp) -> Any:
    lhs = _convert_predicate(a.left)
    rhs = _convert_predicate(a.right)

    op = a.op
    if isinstance(op, BitAnd):
        return pyiceberg.expressions.And(lhs, rhs)
    if isinstance(op, BitOr):
        return pyiceberg.expressions.Or(lhs, rhs)
    else:
        msg = f"Unknown: {lhs} {op} {rhs}"
        raise TypeError(msg)


@_convert_predicate.register(Compare)
def _(a: Compare) -> Any:
    op = a.ops[0]
    lhs = _convert_predicate(a.left)[0]
    rhs = _convert_predicate(a.comparators[0])

    if isinstance(op, Gt):
        return pyiceberg.expressions.GreaterThan(lhs, rhs)
    if isinstance(op, GtE):
        return pyiceberg.expressions.GreaterThanOrEqual(lhs, rhs)
    if isinstance(op, Eq):
        return pyiceberg.expressions.EqualTo(lhs, rhs)
    if isinstance(op, Lt):
        return pyiceberg.expressions.LessThan(lhs, rhs)
    if isinstance(op, LtE):
        return pyiceberg.expressions.LessThanOrEqual(lhs, rhs)
    else:
        msg = f"Unknown comparison: {op}"
        raise TypeError(msg)


@_convert_predicate.register(List)
def _(a: List) -> Any:
    return [_convert_predicate(e) for e in a.elts]


class IdentityTransformedPartitionValuesBuilder:
    def __init__(
        self,
        table: Table,
        projected_schema: pyiceberg.schema.Schema,
    ) -> None:
        import pyiceberg.schema
        from pyiceberg.io.pyarrow import schema_to_pyarrow
        from pyiceberg.transforms import IdentityTransform
        from pyiceberg.types import (
            DoubleType,
            FloatType,
            IntegerType,
            LongType,
        )

        projected_ids: set[int] = projected_schema.field_ids

        # {source_field_id: [values] | error_message}
        self.partition_values: dict[int, list[Any] | str] = {}
        # Logical types will have length-2 list [<constructor type>, <cast type>].
        # E.g. for Datetime it will be [Int64, Datetime]
        self.partition_values_dtypes: dict[int, pl.DataType] = {}

        # {spec_id: [partition_value_index, source_field_id]}
        self.partition_spec_id_to_identity_transforms: dict[
            int, list[tuple[int, int]]
        ] = {}

        partition_specs = table.specs()

        for spec_id, spec in partition_specs.items():
            out = []

            for field_index, field in enumerate(spec.fields):
                if field.source_id in projected_ids and isinstance(
                    field.transform, IdentityTransform
                ):
                    out.append((field_index, field.source_id))
                    self.partition_values[field.source_id] = []

            self.partition_spec_id_to_identity_transforms[spec_id] = out

        for field_id in self.partition_values:
            projected_field = projected_schema.find_field(field_id)
            projected_type = projected_field.field_type

            _, output_dtype = pl.Schema(
                schema_to_pyarrow(pyiceberg.schema.Schema(projected_field))
            ).popitem()

            self.partition_values_dtypes[field_id] = output_dtype

            if not projected_type.is_primitive or output_dtype.is_nested():
                self.partition_values[field_id] = (
                    f"non-primitive type: {projected_type = } {output_dtype = }"
                )

            for schema in table.schemas().values():
                try:
                    type_this_schema = schema.find_field(field_id).field_type
                except ValueError:
                    continue

                if not (
                    projected_type == type_this_schema
                    or (
                        isinstance(projected_type, LongType)
                        and isinstance(type_this_schema, IntegerType)
                    )
                    or (
                        isinstance(projected_type, (DoubleType, FloatType))
                        and isinstance(type_this_schema, (DoubleType, FloatType))
                    )
                ):
                    self.partition_values[field_id] = (
                        f"unsupported type change: from: {type_this_schema}, "
                        f"to: {projected_type}"
                    )

    def push_partition_values(
        self,
        *,
        current_index: int,
        partition_spec_id: int,
        partition_values: pyiceberg.typedef.Record,
    ) -> None:
        try:
            identity_transforms = self.partition_spec_id_to_identity_transforms[
                partition_spec_id
            ]
        except KeyError:
            self.partition_values = {
                k: f"partition spec ID not found: {partition_spec_id}"
                for k in self.partition_values
            }
            return

        for i, source_field_id in identity_transforms:
            partition_value = partition_values[i]

            if isinstance(values := self.partition_values[source_field_id], list):
                # extend() - there can be gaps from partitions being
                # added/removed/re-added
                values.extend(None for _ in range(current_index - len(values)))
                values.append(partition_value)

    def finish(self) -> dict[int, pl.Series | str]:
        from polars.datatypes import Date, Datetime, Duration, Int32, Int64, Time

        out: dict[int, pl.Series | str] = {}

        for field_id, v in self.partition_values.items():
            if isinstance(v, str):
                out[field_id] = v
            else:
                try:
                    output_dtype = self.partition_values_dtypes[field_id]

                    constructor_dtype = (
                        Int64
                        if isinstance(output_dtype, (Datetime, Duration, Time))
                        else Int32
                        if isinstance(output_dtype, Date)
                        else output_dtype
                    )

                    s = pl.Series(v, dtype=constructor_dtype)

                    assert not s.dtype.is_nested()

                    if isinstance(output_dtype, Time):
                        # Physical from PyIceberg is in microseconds, physical
                        # used by polars is in nanoseconds.
                        s = s * ICEBERG_TIME_TO_NS

                    s = s.cast(output_dtype)

                    out[field_id] = s

                except Exception as e:
                    out[field_id] = f"failed to load partition values: {e}"

        return out


class IcebergStatisticsLoader:
    def __init__(
        self,
        table: Table,
        projected_filter_schema: pyiceberg.schema.Schema,
    ) -> None:
        import pyiceberg.schema
        from pyiceberg.io.pyarrow import schema_to_pyarrow

        import polars as pl
        import polars._utils.logging

        verbose = polars._utils.logging.verbose()

        self.file_column_statistics: dict[int, IcebergColumnStatisticsLoader] = {}
        self.load_as_empty_statistics: list[str] = []
        self.file_lengths: list[int] = []
        self.projected_filter_schema = projected_filter_schema

        for field in projected_filter_schema.fields:
            field_all_types = set()

            for schema in table.schemas().values():
                with contextlib.suppress(ValueError):
                    field_all_types.add(schema.find_field(field.field_id).field_type)

            _, field_polars_dtype = pl.Schema(
                schema_to_pyarrow(pyiceberg.schema.Schema(field))
            ).popitem()

            load_from_bytes_impl = LoadFromBytesImpl.init_for_field_type(
                field.field_type,
                field_all_types,
                field_polars_dtype,
            )

            if verbose:
                _load_from_bytes_impl = (
                    type(load_from_bytes_impl).__name__
                    if load_from_bytes_impl is not None
                    else "None"
                )

                eprint(
                    "IcebergStatisticsLoader: "
                    f"{field.name = }, "
                    f"{field.field_id = }, "
                    f"{field.field_type = }, "
                    f"{field_all_types = }, "
                    f"{field_polars_dtype = }, "
                    f"{_load_from_bytes_impl = }"
                )

            self.file_column_statistics[field.field_id] = IcebergColumnStatisticsLoader(
                field_id=field.field_id,
                column_name=field.name,
                column_dtype=field_polars_dtype,
                load_from_bytes_impl=load_from_bytes_impl,
                min_values=[],
                max_values=[],
                null_count=[],
            )

    def push_file_statistics(self, file: DataFile) -> None:
        self.file_lengths.append(file.record_count)

        for stats in self.file_column_statistics.values():
            stats.push_file_statistics(file)

    def finish(
        self,
        expected_height: int,
        identity_transformed_values: dict[int, pl.Series | str],
    ) -> pl.DataFrame:
        import polars as pl

        out: list[pl.DataFrame] = [
            pl.Series("len", self.file_lengths, dtype=pl.UInt32).to_frame()
        ]

        for field_id, stat_builder in self.file_column_statistics.items():
            if (p := identity_transformed_values.get(field_id)) is not None:
                if isinstance(p, str):
                    msg = f"statistics load failure for filter column: {p}"
                    raise ComputeError(msg)

            column_stats_df = stat_builder.finish(expected_height, p)
            out.append(column_stats_df)

        return pl.concat(out, how="horizontal")


@dataclass
class IcebergColumnStatisticsLoader:
    column_name: str
    column_dtype: pl.DataType
    field_id: int
    load_from_bytes_impl: LoadFromBytesImpl | None
    null_count: list[int | None]
    min_values: list[bytes | None]
    max_values: list[bytes | None]

    def push_file_statistics(self, file: DataFile) -> None:
        self.null_count.append(file.null_value_counts.get(self.field_id))

        if self.load_from_bytes_impl is not None:
            self.min_values.append(file.lower_bounds.get(self.field_id))
            self.max_values.append(file.upper_bounds.get(self.field_id))

    def finish(
        self,
        expected_height: int,
        identity_transformed_values: pl.Series | None,
    ) -> pl.DataFrame:
        import polars as pl

        c = self.column_name
        assert len(self.null_count) == expected_height

        out = pl.Series(f"{c}_nc", self.null_count, dtype=pl.UInt32).to_frame()

        if self.load_from_bytes_impl is None:
            s = (
                identity_transformed_values
                if identity_transformed_values is not None
                else pl.repeat(None, expected_height, dtype=self.column_dtype)
            )

            return out.with_columns(s.alias(f"{c}_min"), s.alias(f"{c}_max"))

        assert len(self.min_values) == expected_height
        assert len(self.max_values) == expected_height

        if self.column_dtype.is_nested():
            raise NotImplementedError

        min_values = self.load_from_bytes_impl.load_from_bytes(self.min_values)
        max_values = self.load_from_bytes_impl.load_from_bytes(self.max_values)

        if identity_transformed_values is not None:
            assert identity_transformed_values.dtype == self.column_dtype

            identity_transformed_values = identity_transformed_values.extend_constant(
                None, expected_height - identity_transformed_values.len()
            )

            min_values = identity_transformed_values.fill_null(min_values)
            max_values = identity_transformed_values.fill_null(max_values)

        return out.with_columns(
            min_values.alias(f"{c}_min"), max_values.alias(f"{c}_max")
        )


# Lazy init instead of global const as PyIceberg is an optional dependency
@cache
def _bytes_loader_lookup() -> dict[
    type[IcebergType],
    tuple[type[LoadFromBytesImpl], type[IcebergType] | Sequence[type[IcebergType]]],
]:
    from pyiceberg.types import (
        BinaryType,
        BooleanType,
        DateType,
        DecimalType,
        FixedType,
        IntegerType,
        LongType,
        StringType,
        TimestampType,
        TimestamptzType,
        TimeType,
    )

    # TODO: Float statistics
    return {
        BooleanType: (LoadBooleanFromBytes, BooleanType),
        DateType: (LoadDateFromBytes, DateType),
        TimeType: (LoadTimeFromBytes, TimeType),
        TimestampType: (LoadTimestampFromBytes, TimestampType),
        TimestamptzType: (LoadTimestamptzFromBytes, TimestamptzType),
        IntegerType: (LoadInt32FromBytes, IntegerType),
        LongType: (LoadInt64FromBytes, (LongType, IntegerType)),
        StringType: (LoadStringFromBytes, StringType),
        BinaryType: (LoadBinaryFromBytes, BinaryType),
        DecimalType: (LoadDecimalFromBytes, DecimalType),
        FixedType: (LoadFixedFromBytes, FixedType),
    }


class LoadFromBytesImpl(abc.ABC):
    def __init__(self, polars_dtype: pl.DataType) -> None:
        self.polars_dtype = polars_dtype

    @staticmethod
    def init_for_field_type(
        current_field_type: IcebergType,
        # All types that this field ID has been set to across schema changes.
        all_field_types: set[IcebergType],
        field_polars_dtype: pl.DataType,
    ) -> LoadFromBytesImpl | None:
        if (v := _bytes_loader_lookup().get(type(current_field_type))) is None:
            return None

        loader_impl, allowed_field_types = v

        return (
            loader_impl(field_polars_dtype)
            if all(isinstance(x, allowed_field_types) for x in all_field_types)  # type: ignore[arg-type]
            else None
        )

    @abc.abstractmethod
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        """`bytes_values` should be of binary type."""


class LoadBinaryFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return pl.Series(byte_values, dtype=pl.Binary)


class LoadDateFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return (
            pl.Series(byte_values, dtype=pl.Binary)
            .bin.reinterpret(dtype=pl.Int32, endianness="little")
            .cast(pl.Date)
        )


class LoadTimeFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return (
            pl.Series(byte_values, dtype=pl.Binary).bin.reinterpret(
                dtype=pl.Int64, endianness="little"
            )
            * ICEBERG_TIME_TO_NS
        ).cast(pl.Time)


class LoadTimestampFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return (
            pl.Series(byte_values, dtype=pl.Binary)
            .bin.reinterpret(dtype=pl.Int64, endianness="little")
            .cast(pl.Datetime("us"))
        )


class LoadTimestamptzFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return (
            pl.Series(byte_values, dtype=pl.Binary)
            .bin.reinterpret(dtype=pl.Int64, endianness="little")
            .cast(pl.Datetime("us", time_zone="UTC"))
        )


class LoadBooleanFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return (
            pl.Series(byte_values, dtype=pl.Binary)
            .bin.reinterpret(dtype=pl.UInt8, endianness="little")
            .cast(pl.Boolean)
        )


class LoadDecimalFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl
        from polars._plr import PySeries

        dtype = self.polars_dtype
        assert isinstance(dtype, pl.Decimal)
        assert dtype.precision is not None

        return wrap_s(
            PySeries._import_decimal_from_iceberg_binary_repr(
                bytes_list=byte_values,
                precision=dtype.precision,
                scale=dtype.scale,
            )
        )


class LoadFixedFromBytes(LoadBinaryFromBytes): ...


class LoadInt32FromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return pl.Series(byte_values, dtype=pl.Binary).bin.reinterpret(
            dtype=pl.Int32, endianness="little"
        )


class LoadInt64FromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        s = pl.Series(byte_values, dtype=pl.Binary)

        return s.bin.reinterpret(dtype=pl.Int64, endianness="little").fill_null(
            s.bin.reinterpret(dtype=pl.Int32, endianness="little").cast(pl.Int64)
        )


class LoadStringFromBytes(LoadFromBytesImpl):
    def load_from_bytes(self, byte_values: list[bytes | None]) -> pl.Series:
        import polars as pl

        return pl.Series(byte_values, dtype=pl.Binary).cast(pl.String)
