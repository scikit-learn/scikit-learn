from __future__ import annotations

from collections.abc import Collection, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Literal, cast, overload

import pyarrow as pa
import pyarrow.compute as pc

from narwhals._arrow.series import ArrowSeries
from narwhals._arrow.utils import concat_tables, native_to_narwhals_dtype, repeat
from narwhals._compliant import EagerDataFrame
from narwhals._expression_parsing import ExprKind
from narwhals._utils import (
    Implementation,
    Version,
    check_column_names_are_unique,
    convert_str_slice_to_int_slice,
    generate_temporary_column_name,
    not_implemented,
    parse_columns_to_drop,
    scale_bytes,
    supports_arrow_c_stream,
    zip_strict,
)
from narwhals.dependencies import is_numpy_array_1d
from narwhals.exceptions import ShapeError

if TYPE_CHECKING:
    from collections.abc import Iterable
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import pandas as pd
    import polars as pl
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._arrow.expr import ArrowExpr
    from narwhals._arrow.group_by import ArrowGroupBy
    from narwhals._arrow.namespace import ArrowNamespace
    from narwhals._arrow.typing import (  # type: ignore[attr-defined]
        ChunkedArrayAny,
        Order,
    )
    from narwhals._compliant.typing import CompliantDataFrameAny, CompliantLazyFrameAny
    from narwhals._spark_like.utils import SparkSession
    from narwhals._translate import IntoArrowTable
    from narwhals._typing import _EagerAllowedImpl, _LazyAllowedImpl
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dtypes import DType
    from narwhals.typing import (
        IntoSchema,
        JoinStrategy,
        SizedMultiIndexSelector,
        SizedMultiNameSelector,
        SizeUnit,
        UniqueKeepStrategy,
        _1DArray,
        _2DArray,
        _SliceIndex,
        _SliceName,
    )

    JoinType: TypeAlias = Literal[
        "left semi",
        "right semi",
        "left anti",
        "right anti",
        "inner",
        "left outer",
        "right outer",
        "full outer",
    ]


class ArrowDataFrame(
    EagerDataFrame["ArrowSeries", "ArrowExpr", "pa.Table", "ChunkedArrayAny"]
):
    _implementation = Implementation.PYARROW

    def __init__(
        self,
        native_dataframe: pa.Table,
        *,
        version: Version,
        validate_column_names: bool,
        validate_backend_version: bool = False,
    ) -> None:
        if validate_column_names:
            check_column_names_are_unique(native_dataframe.column_names)
        if validate_backend_version:
            self._validate_backend_version()
        self._native_frame = native_dataframe
        self._version = version

    @classmethod
    def from_arrow(cls, data: IntoArrowTable, /, *, context: _LimitedContext) -> Self:
        backend_version = context._implementation._backend_version()
        if cls._is_native(data):
            native = data
        elif backend_version >= (14,) or isinstance(data, Collection):
            native = pa.table(data)
        elif supports_arrow_c_stream(data):  # pragma: no cover
            msg = f"'pyarrow>=14.0.0' is required for `from_arrow` for object of type {type(data).__name__!r}."
            raise ModuleNotFoundError(msg)
        else:  # pragma: no cover
            msg = f"`from_arrow` is not supported for object of type {type(data).__name__!r}."
            raise TypeError(msg)
        return cls.from_native(native, context=context)

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | None,
    ) -> Self:
        from narwhals.schema import Schema

        pa_schema = Schema(schema).to_arrow() if schema is not None else schema
        if pa_schema and not data:
            native = pa_schema.empty_table()
        else:
            native = pa.Table.from_pydict(data, schema=pa_schema)
        return cls.from_native(native, context=context)

    @classmethod
    def from_dicts(
        cls,
        data: Sequence[Mapping[str, Any]],
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | None,
    ) -> Self:
        from narwhals.schema import Schema

        pa_schema = Schema(schema).to_arrow() if schema is not None else schema
        if pa_schema and not data:
            native = pa_schema.empty_table()
        else:
            native = pa.Table.from_pylist(data, schema=pa_schema)
        return cls.from_native(native, context=context)

    @staticmethod
    def _is_native(obj: pa.Table | Any) -> TypeIs[pa.Table]:
        return isinstance(obj, pa.Table)

    @classmethod
    def from_native(cls, data: pa.Table, /, *, context: _LimitedContext) -> Self:
        return cls(data, version=context._version, validate_column_names=True)

    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        /,
        *,
        context: _LimitedContext,
        schema: IntoSchema | Sequence[str] | None,
    ) -> Self:
        from narwhals.schema import Schema

        arrays = [pa.array(val) for val in data.T]
        if isinstance(schema, (Mapping, Schema)):
            native = pa.Table.from_arrays(arrays, schema=Schema(schema).to_arrow())
        else:
            native = pa.Table.from_arrays(arrays, cls._numpy_column_names(data, schema))
        return cls.from_native(native, context=context)

    def __narwhals_namespace__(self) -> ArrowNamespace:
        from narwhals._arrow.namespace import ArrowNamespace

        return ArrowNamespace(version=self._version)

    def __native_namespace__(self) -> ModuleType:
        if self._implementation is Implementation.PYARROW:
            return self._implementation.to_native_namespace()

        msg = f"Expected pyarrow, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    def __narwhals_dataframe__(self) -> Self:
        return self

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, version=version, validate_column_names=False)

    def _with_native(self, df: pa.Table, *, validate_column_names: bool = True) -> Self:
        return self.__class__(
            df, version=self._version, validate_column_names=validate_column_names
        )

    @property
    def shape(self) -> tuple[int, int]:
        return self.native.shape

    def __len__(self) -> int:
        return len(self.native)

    def row(self, index: int) -> tuple[Any, ...]:
        return tuple(col[index] for col in self.native.itercolumns())

    @overload
    def rows(self, *, named: Literal[True]) -> list[dict[str, Any]]: ...

    @overload
    def rows(self, *, named: Literal[False]) -> list[tuple[Any, ...]]: ...

    @overload
    def rows(self, *, named: bool) -> list[tuple[Any, ...]] | list[dict[str, Any]]: ...

    def rows(self, *, named: bool) -> list[tuple[Any, ...]] | list[dict[str, Any]]:
        if not named:
            return list(self.iter_rows(named=False, buffer_size=512))  # type: ignore[return-value]
        return self.native.to_pylist()

    def iter_columns(self) -> Iterator[ArrowSeries]:
        for name, series in zip_strict(self.columns, self.native.itercolumns()):
            yield ArrowSeries.from_native(series, context=self, name=name)

    _iter_columns = iter_columns

    def iter_rows(
        self, *, named: bool, buffer_size: int
    ) -> Iterator[tuple[Any, ...]] | Iterator[dict[str, Any]]:
        df = self.native
        num_rows = df.num_rows

        if not named:
            for i in range(0, num_rows, buffer_size):
                rows = df[i : i + buffer_size].to_pydict().values()
                yield from zip_strict(*rows)
        else:
            for i in range(0, num_rows, buffer_size):
                yield from df[i : i + buffer_size].to_pylist()

    def get_column(self, name: str) -> ArrowSeries:
        if not isinstance(name, str):
            msg = f"Expected str, got: {type(name)}"
            raise TypeError(msg)
        return ArrowSeries.from_native(self.native[name], context=self, name=name)

    def __array__(self, dtype: Any, *, copy: bool | None) -> _2DArray:
        return self.native.__array__(dtype, copy=copy)

    def _gather(self, rows: SizedMultiIndexSelector[ChunkedArrayAny]) -> Self:
        if len(rows) == 0:
            return self._with_native(self.native.slice(0, 0))
        if self._backend_version < (18,) and isinstance(rows, tuple):
            rows = list(rows)
        return self._with_native(self.native.take(rows))

    def _gather_slice(self, rows: _SliceIndex | range) -> Self:
        start = rows.start or 0
        stop = rows.stop if rows.stop is not None else len(self.native)
        if start < 0:
            start = len(self.native) + start
        if stop < 0:
            stop = len(self.native) + stop
        if rows.step is not None and rows.step != 1:
            msg = "Slicing with step is not supported on PyArrow tables"
            raise NotImplementedError(msg)
        return self._with_native(self.native.slice(start, stop - start))

    def _select_slice_name(self, columns: _SliceName) -> Self:
        start, stop, step = convert_str_slice_to_int_slice(columns, self.columns)
        return self._with_native(self.native.select(self.columns[start:stop:step]))

    def _select_slice_index(self, columns: _SliceIndex | range) -> Self:
        return self._with_native(
            self.native.select(self.columns[columns.start : columns.stop : columns.step])
        )

    def _select_multi_index(
        self, columns: SizedMultiIndexSelector[ChunkedArrayAny]
    ) -> Self:
        selector: Sequence[int]
        if isinstance(columns, pa.ChunkedArray):
            # TODO @dangotbanned: Fix upstream with `pa.ChunkedArray.to_pylist(self) -> list[Any]:`
            selector = cast("Sequence[int]", columns.to_pylist())
        # TODO @dangotbanned: Fix upstream, it is actually much narrower
        # **Doesn't accept `ndarray`**
        elif is_numpy_array_1d(columns):
            selector = columns.tolist()
        else:
            selector = columns
        return self._with_native(self.native.select(selector))

    def _select_multi_name(
        self, columns: SizedMultiNameSelector[ChunkedArrayAny]
    ) -> Self:
        selector: Sequence[str] | _1DArray
        if isinstance(columns, pa.ChunkedArray):
            # TODO @dangotbanned: Fix upstream with `pa.ChunkedArray.to_pylist(self) -> list[Any]:`
            selector = cast("Sequence[str]", columns.to_pylist())
        else:
            selector = columns
        # NOTE: Fixed in https://github.com/zen-xu/pyarrow-stubs/pull/221
        return self._with_native(self.native.select(selector))  # pyright: ignore[reportArgumentType]

    @property
    def schema(self) -> dict[str, DType]:
        return {
            field.name: native_to_narwhals_dtype(field.type, self._version)
            for field in self.native.schema
        }

    def collect_schema(self) -> dict[str, DType]:
        return self.schema

    def estimated_size(self, unit: SizeUnit) -> int | float:
        sz = self.native.nbytes
        return scale_bytes(sz, unit)

    explode = not_implemented()

    @property
    def columns(self) -> list[str]:
        return self.native.column_names

    def simple_select(self, *column_names: str) -> Self:
        return self._with_native(
            self.native.select(list(column_names)), validate_column_names=False
        )

    def select(self, *exprs: ArrowExpr) -> Self:
        new_series = self._evaluate_into_exprs(*exprs)
        if not new_series:
            # return empty dataframe, like Polars does
            return self._with_native(
                self.native.__class__.from_arrays([]), validate_column_names=False
            )
        names = [s.name for s in new_series]
        align = new_series[0]._align_full_broadcast
        reshaped = align(*new_series)
        df = pa.Table.from_arrays([s.native for s in reshaped], names=names)
        return self._with_native(df, validate_column_names=True)

    def _extract_comparand(self, other: ArrowSeries) -> ChunkedArrayAny:
        length = len(self)
        if not other._broadcast:
            if (len_other := len(other)) != length:
                msg = f"Expected object of length {length}, got: {len_other}."
                raise ShapeError(msg)
            return other.native

        value = other.native[0]
        return pa.chunked_array([pa.repeat(value, length)])

    def with_columns(self, *exprs: ArrowExpr) -> Self:
        # NOTE: We use a faux-mutable variable and repeatedly "overwrite" (native_frame)
        # All `pyarrow` data is immutable, so this is fine
        native_frame = self.native
        new_columns = self._evaluate_into_exprs(*exprs)
        columns = self.columns

        for col_value in new_columns:
            col_name = col_value.name
            column = self._extract_comparand(col_value)
            native_frame = (
                native_frame.set_column(columns.index(col_name), col_name, column=column)
                if col_name in columns
                else native_frame.append_column(col_name, column=column)
            )

        return self._with_native(native_frame, validate_column_names=False)

    def group_by(
        self, keys: Sequence[str] | Sequence[ArrowExpr], *, drop_null_keys: bool
    ) -> ArrowGroupBy:
        from narwhals._arrow.group_by import ArrowGroupBy

        return ArrowGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def join(
        self,
        other: Self,
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        how_to_join_map: dict[str, JoinType] = {
            "anti": "left anti",
            "semi": "left semi",
            "inner": "inner",
            "left": "left outer",
            "full": "full outer",
        }

        if how == "cross":
            plx = self.__narwhals_namespace__()
            key_token = generate_temporary_column_name(
                n_bytes=8, columns=[*self.columns, *other.columns]
            )

            return self._with_native(
                self.with_columns(
                    plx.lit(0, None).alias(key_token).broadcast(ExprKind.LITERAL)
                )
                .native.join(
                    other.with_columns(
                        plx.lit(0, None).alias(key_token).broadcast(ExprKind.LITERAL)
                    ).native,
                    keys=key_token,
                    right_keys=key_token,
                    join_type="inner",
                    right_suffix=suffix,
                )
                .drop([key_token])
            )

        coalesce_keys = how != "full"  # polars full join does not coalesce keys
        return self._with_native(
            self.native.join(
                other.native,
                keys=left_on or [],  # type: ignore[arg-type]
                right_keys=right_on,  # type: ignore[arg-type]
                join_type=how_to_join_map[how],
                right_suffix=suffix,
                coalesce_keys=coalesce_keys,
            )
        )

    join_asof = not_implemented()

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        to_drop = parse_columns_to_drop(self, columns, strict=strict)
        return self._with_native(self.native.drop(to_drop), validate_column_names=False)

    def drop_nulls(self, subset: Sequence[str] | None) -> Self:
        if subset is None:
            return self._with_native(self.native.drop_null(), validate_column_names=False)
        plx = self.__narwhals_namespace__()
        mask = ~plx.any_horizontal(plx.col(*subset).is_null(), ignore_nulls=True)
        return self.filter(mask)

    def sort(self, *by: str, descending: bool | Sequence[bool], nulls_last: bool) -> Self:
        if isinstance(descending, bool):
            order: Order = "descending" if descending else "ascending"
            sorting: list[tuple[str, Order]] = [(key, order) for key in by]
        else:
            sorting = [
                (key, "descending" if is_descending else "ascending")
                for key, is_descending in zip_strict(by, descending)
            ]

        null_placement = "at_end" if nulls_last else "at_start"

        return self._with_native(
            self.native.sort_by(sorting, null_placement=null_placement),
            validate_column_names=False,
        )

    def top_k(self, k: int, *, by: Iterable[str], reverse: bool | Sequence[bool]) -> Self:
        if isinstance(reverse, bool):
            order: Order = "ascending" if reverse else "descending"
            sorting: list[tuple[str, Order]] = [(key, order) for key in by]
        else:
            sorting = [
                (key, "ascending" if is_ascending else "descending")
                for key, is_ascending in zip_strict(by, reverse)
            ]
        return self._with_native(
            self.native.take(pc.select_k_unstable(self.native, k, sorting)),  # type: ignore[call-overload]
            validate_column_names=False,
        )

    def to_pandas(self) -> pd.DataFrame:
        return self.native.to_pandas()

    def to_polars(self) -> pl.DataFrame:
        import polars as pl  # ignore-banned-import

        return pl.from_arrow(self.native)  # type: ignore[return-value]

    def to_numpy(self, dtype: Any = None, *, copy: bool | None = None) -> _2DArray:
        import numpy as np  # ignore-banned-import

        arr: Any = np.column_stack([col.to_numpy() for col in self.native.columns])
        return arr

    @overload
    def to_dict(self, *, as_series: Literal[True]) -> dict[str, ArrowSeries]: ...

    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...

    def to_dict(
        self, *, as_series: bool
    ) -> dict[str, ArrowSeries] | dict[str, list[Any]]:
        it = self.iter_columns()
        if as_series:
            return {ser.name: ser for ser in it}
        return {ser.name: ser.to_list() for ser in it}

    def with_row_index(self, name: str, order_by: Sequence[str] | None) -> Self:
        plx = self.__narwhals_namespace__()
        if order_by is None:
            import numpy as np  # ignore-banned-import

            data = pa.array(np.arange(len(self), dtype=np.int64))
            row_index = plx._expr._from_series(
                plx._series.from_iterable(data, context=self, name=name)
            )
        else:
            rank = plx.col(order_by[0]).rank("ordinal", descending=False)
            row_index = (rank.over(partition_by=[], order_by=order_by) - 1).alias(name)
        return self.select(row_index, plx.all())

    def filter(self, predicate: ArrowExpr) -> Self:
        # `[0]` is safe as the predicate's expression only returns a single column
        mask_native = self._evaluate_into_exprs(predicate)[0].native
        return self._with_native(
            self.native.filter(mask_native), validate_column_names=False
        )

    def head(self, n: int) -> Self:
        df = self.native
        if n >= 0:
            return self._with_native(df.slice(0, n), validate_column_names=False)
        num_rows = df.num_rows
        return self._with_native(
            df.slice(0, max(0, num_rows + n)), validate_column_names=False
        )

    def tail(self, n: int) -> Self:
        df = self.native
        if n >= 0:
            num_rows = df.num_rows
            return self._with_native(
                df.slice(max(0, num_rows - n)), validate_column_names=False
            )
        return self._with_native(df.slice(abs(n)), validate_column_names=False)

    def lazy(
        self,
        backend: _LazyAllowedImpl | None = None,
        *,
        session: SparkSession | None = None,
    ) -> CompliantLazyFrameAny:
        if backend is None:
            return self
        if backend is Implementation.DUCKDB:
            import duckdb  # ignore-banned-import

            from narwhals._duckdb.dataframe import DuckDBLazyFrame

            _df = self.native
            return DuckDBLazyFrame(
                duckdb.table("_df"), validate_backend_version=True, version=self._version
            )
        if backend is Implementation.POLARS:
            import polars as pl  # ignore-banned-import

            from narwhals._polars.dataframe import PolarsLazyFrame

            return PolarsLazyFrame(
                cast("pl.DataFrame", pl.from_arrow(self.native)).lazy(),
                validate_backend_version=True,
                version=self._version,
            )
        if backend is Implementation.DASK:
            import dask.dataframe as dd  # ignore-banned-import

            from narwhals._dask.dataframe import DaskLazyFrame

            return DaskLazyFrame(
                dd.from_pandas(self.native.to_pandas()),
                validate_backend_version=True,
                version=self._version,
            )
        if backend is Implementation.IBIS:
            import ibis  # ignore-banned-import

            from narwhals._ibis.dataframe import IbisLazyFrame

            return IbisLazyFrame(
                ibis.memtable(self.native, columns=self.columns),
                validate_backend_version=True,
                version=self._version,
            )

        if backend.is_spark_like():
            from narwhals._spark_like.dataframe import SparkLikeLazyFrame

            if session is None:
                msg = "Spark like backends require `session` to be not None."
                raise ValueError(msg)

            return SparkLikeLazyFrame._from_compliant_dataframe(
                self, session=session, implementation=backend, version=self._version
            )

        raise AssertionError  # pragma: no cover

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        if backend is Implementation.PYARROW or backend is None:
            from narwhals._arrow.dataframe import ArrowDataFrame

            return ArrowDataFrame(
                self.native, version=self._version, validate_column_names=False
            )

        if backend is Implementation.PANDAS:
            from narwhals._pandas_like.dataframe import PandasLikeDataFrame

            return PandasLikeDataFrame(
                self.native.to_pandas(),
                implementation=Implementation.PANDAS,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=False,
            )

        if backend is Implementation.POLARS:
            import polars as pl  # ignore-banned-import

            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame(
                cast("pl.DataFrame", pl.from_arrow(self.native)),
                validate_backend_version=True,
                version=self._version,
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise AssertionError(msg)  # pragma: no cover

    def clone(self) -> Self:
        return self._with_native(self.native, validate_column_names=False)

    def item(self, row: int | None, column: int | str | None) -> Any:
        from narwhals._arrow.series import maybe_extract_py_scalar

        if row is None and column is None:
            if self.shape != (1, 1):
                msg = (
                    "can only call `.item()` if the dataframe is of shape (1, 1),"
                    " or if explicit row/col values are provided;"
                    f" frame has shape {self.shape!r}"
                )
                raise ValueError(msg)
            return maybe_extract_py_scalar(self.native[0][0], return_py_scalar=True)

        if row is None or column is None:
            msg = "cannot call `.item()` with only one of `row` or `column`"
            raise ValueError(msg)

        _col = self.columns.index(column) if isinstance(column, str) else column
        return maybe_extract_py_scalar(self.native[_col][row], return_py_scalar=True)

    def rename(self, mapping: Mapping[str, str]) -> Self:
        names: dict[str, str] | list[str]
        if self._backend_version >= (17,):
            names = cast("dict[str, str]", mapping)
        else:  # pragma: no cover
            names = [mapping.get(c, c) for c in self.columns]
        return self._with_native(self.native.rename_columns(names))

    def write_parquet(self, file: str | Path | BytesIO) -> None:
        import pyarrow.parquet as pp

        pp.write_table(self.native, file)

    @overload
    def write_csv(self, file: None) -> str: ...

    @overload
    def write_csv(self, file: str | Path | BytesIO) -> None: ...

    def write_csv(self, file: str | Path | BytesIO | None) -> str | None:
        import pyarrow.csv as pa_csv

        if file is None:
            csv_buffer = pa.BufferOutputStream()
            pa_csv.write_csv(self.native, csv_buffer)
            return csv_buffer.getvalue().to_pybytes().decode()
        pa_csv.write_csv(self.native, file)
        return None

    def is_unique(self) -> ArrowSeries:
        import numpy as np  # ignore-banned-import

        col_token = generate_temporary_column_name(n_bytes=8, columns=self.columns)
        row_index = pa.array(np.arange(len(self)))
        keep_idx = (
            self.native.append_column(col_token, row_index)
            .group_by(self.columns)
            .aggregate([(col_token, "min"), (col_token, "max")])
        )
        native = pa.chunked_array(
            pc.and_(
                pc.is_in(row_index, keep_idx[f"{col_token}_min"]),
                pc.is_in(row_index, keep_idx[f"{col_token}_max"]),
            )
        )
        return ArrowSeries.from_native(native, context=self)

    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        maintain_order: bool | None = None,
        order_by: Sequence[str] | None,
    ) -> Self:
        # The param `maintain_order` is only here for compatibility with the Polars API
        # and has no effect on the output.
        import numpy as np  # ignore-banned-import

        if subset and (error := self._check_columns_exist(subset)):
            raise error
        subset = list(subset or self.columns)

        if keep in {"any", "first", "last"}:
            from narwhals._arrow.group_by import ArrowGroupBy

            agg_func = ArrowGroupBy._REMAP_UNIQUE[keep]
            col_token = generate_temporary_column_name(n_bytes=8, columns=self.columns)
            if order_by and maintain_order:
                idx_token = generate_temporary_column_name(
                    n_bytes=8, columns=[*self.columns, col_token]
                )
                df = (
                    self.with_row_index(idx_token, order_by=None)
                    .sort(*order_by, nulls_last=False, descending=False)
                    .unique(subset=subset, keep=keep, maintain_order=False, order_by=None)
                )
                return df.sort(idx_token, descending=False, nulls_last=False).drop(
                    [idx_token], strict=False
                )
            if order_by:
                native = self.sort(*order_by, nulls_last=False, descending=False).native
            else:
                native = self.native
            keep_idx_native = (
                native.append_column(col_token, pa.array(np.arange(len(self))))
                .group_by(subset)
                .aggregate([(col_token, agg_func)])
                .column(f"{col_token}_{agg_func}")
            )
            return self._with_native(
                native.take(keep_idx_native), validate_column_names=False
            )

        keep_idx = self.simple_select(*subset).is_unique()
        plx = self.__narwhals_namespace__()
        return self.filter(plx._expr._from_series(keep_idx))

    def gather_every(self, n: int, offset: int) -> Self:
        return self._with_native(self.native[offset::n], validate_column_names=False)

    def to_arrow(self) -> pa.Table:
        return self.native

    def sample(
        self,
        n: int | None,
        *,
        fraction: float | None,
        with_replacement: bool,
        seed: int | None,
    ) -> Self:
        import numpy as np  # ignore-banned-import

        num_rows = len(self)
        if n is None and fraction is not None:
            n = int(num_rows * fraction)
        rng = np.random.default_rng(seed=seed)
        idx = np.arange(num_rows)
        mask = rng.choice(idx, size=n, replace=with_replacement)
        return self._with_native(self.native.take(mask), validate_column_names=False)

    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        # TODO(Unassigned): Even with promote_options="permissive", pyarrow does not
        # upcast numeric to non-numeric (e.g. string) datatypes
        n = len(self)
        index = [] if index is None else list(index)
        on_ = (c for c in self.columns if c not in index) if on is None else iter(on)
        index_cols = self.native.select(index)
        column = self.native.column
        tables = (
            index_cols.append_column(variable_name, repeat(name, n)).append_column(
                value_name, column(name)
            )
            for name in on_
        )
        return self._with_native(concat_tables(tables, "permissive"))

    pivot = not_implemented()
