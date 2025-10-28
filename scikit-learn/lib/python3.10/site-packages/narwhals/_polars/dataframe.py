from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence, Sized
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar, cast, overload

import polars as pl

from narwhals._polars.namespace import PolarsNamespace
from narwhals._polars.series import PolarsSeries
from narwhals._polars.utils import (
    FROM_DICTS_ACCEPTS_MAPPINGS,
    catch_polars_exception,
    extract_args_kwargs,
    native_to_narwhals_dtype,
)
from narwhals._utils import (
    Implementation,
    _into_arrow_table,
    convert_str_slice_to_int_slice,
    generate_temporary_column_name,
    is_compliant_series,
    is_index_selector,
    is_range,
    is_sequence_like,
    is_slice_index,
    is_slice_none,
    parse_columns_to_drop,
    requires,
)
from narwhals.dependencies import is_numpy_array_1d
from narwhals.exceptions import ColumnNotFoundError

if TYPE_CHECKING:
    from collections.abc import Iterable
    from types import ModuleType
    from typing import Callable

    import pandas as pd
    import pyarrow as pa
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._compliant.typing import CompliantDataFrameAny, CompliantLazyFrameAny
    from narwhals._polars.expr import PolarsExpr
    from narwhals._polars.group_by import PolarsGroupBy, PolarsLazyGroupBy
    from narwhals._spark_like.utils import SparkSession
    from narwhals._translate import IntoArrowTable
    from narwhals._typing import _EagerAllowedImpl, _LazyAllowedImpl
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dataframe import DataFrame, LazyFrame
    from narwhals.dtypes import DType
    from narwhals.typing import (
        IntoSchema,
        JoinStrategy,
        MultiColSelector,
        MultiIndexSelector,
        PivotAgg,
        SingleIndexSelector,
        UniqueKeepStrategy,
        _2DArray,
    )

    T = TypeVar("T")
    R = TypeVar("R")

Method: TypeAlias = "Callable[..., R]"
"""Generic alias representing all methods implemented via `__getattr__`.

Where `R` is the return type.
"""

# DataFrame methods where PolarsDataFrame just defers to Polars.DataFrame directly.
INHERITED_METHODS = frozenset(
    [
        "clone",
        "drop_nulls",
        "estimated_size",
        "explode",
        "filter",
        "gather_every",
        "head",
        "is_unique",
        "item",
        "iter_rows",
        "join_asof",
        "rename",
        "row",
        "rows",
        "sample",
        "select",
        "sink_parquet",
        "sort",
        "tail",
        "to_arrow",
        "to_pandas",
        "with_columns",
        "write_csv",
        "write_parquet",
    ]
)

NativePolarsFrame = TypeVar("NativePolarsFrame", pl.DataFrame, pl.LazyFrame)


class PolarsBaseFrame(Generic[NativePolarsFrame]):
    drop_nulls: Method[Self]
    explode: Method[Self]
    filter: Method[Self]
    gather_every: Method[Self]
    head: Method[Self]
    join_asof: Method[Self]
    rename: Method[Self]
    select: Method[Self]
    sort: Method[Self]
    tail: Method[Self]
    with_columns: Method[Self]

    _native_frame: NativePolarsFrame
    _implementation = Implementation.POLARS
    _version: Version

    def __init__(
        self,
        df: NativePolarsFrame,
        *,
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        self._native_frame = df
        self._version = version
        if validate_backend_version:
            self._validate_backend_version()

    def _validate_backend_version(self) -> None:
        """Raise if installed version below `nw._utils.MIN_VERSIONS`.

        **Only use this when moving between backends.**
        Otherwise, the validation will have taken place already.
        """
        _ = self._implementation._backend_version()

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    @property
    def native(self) -> NativePolarsFrame:
        return self._native_frame

    @property
    def columns(self) -> list[str]:
        return self.native.columns

    def __narwhals_namespace__(self) -> PolarsNamespace:
        return PolarsNamespace(version=self._version)

    def __native_namespace__(self) -> ModuleType:
        if self._implementation is Implementation.POLARS:
            return self._implementation.to_native_namespace()

        msg = f"Expected polars, got: {type(self._implementation)}"  # pragma: no cover
        raise AssertionError(msg)

    def _with_native(self, df: NativePolarsFrame) -> Self:
        return self.__class__(df, version=self._version)

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, version=version)

    @classmethod
    def from_native(cls, data: NativePolarsFrame, /, *, context: _LimitedContext) -> Self:
        return cls(data, version=context._version)

    def simple_select(self, *column_names: str) -> Self:
        return self._with_native(self.native.select(*column_names))

    def aggregate(self, *exprs: Any) -> Self:
        return self.select(*exprs)

    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        maintain_order: bool | None = None,
        order_by: Sequence[str] | None = None,
    ) -> Self:
        if order_by and maintain_order:
            token = generate_temporary_column_name(8, self.columns, prefix="row_index_")
            res = (
                self.native.with_row_index(token)
                .sort(order_by, nulls_last=False)
                .unique(subset or self.columns, keep=keep)
                .sort(token)
                .drop(token)
            )
        elif order_by:
            res = self.native.sort(order_by).unique(subset, keep=keep)
        else:
            res = self.native.unique(
                subset, keep=keep, maintain_order=maintain_order or False
            )
        return self._with_native(res)

    @property
    def schema(self) -> dict[str, DType]:
        return self.collect_schema()

    def join(
        self,
        other: PolarsBaseFrame[NativePolarsFrame],
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        how_native = (
            "outer" if (self._backend_version < (0, 20, 29) and how == "full") else how
        )
        return self._with_native(
            self.native.join(
                other=other.native,
                how=how_native,  # type: ignore[arg-type]
                left_on=left_on,
                right_on=right_on,
                suffix=suffix,
            )
        )

    def top_k(
        self, k: int, *, by: str | Iterable[str], reverse: bool | Sequence[bool]
    ) -> Self:
        if self._backend_version < (1, 0, 0):
            return self._with_native(
                self.native.top_k(
                    k=k,
                    by=by,
                    descending=reverse,  # type: ignore[call-arg]
                )
            )
        return self._with_native(self.native.top_k(k=k, by=by, reverse=reverse))

    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        if self._backend_version < (1, 0, 0):
            return self._with_native(
                self.native.melt(
                    id_vars=index,
                    value_vars=on,
                    variable_name=variable_name,
                    value_name=value_name,
                )
            )
        return self._with_native(
            self.native.unpivot(
                on=on, index=index, variable_name=variable_name, value_name=value_name
            )
        )

    def collect_schema(self) -> dict[str, DType]:
        df = self.native
        schema = df.schema if self._backend_version < (1,) else df.collect_schema()
        return {
            name: native_to_narwhals_dtype(dtype, self._version)
            for name, dtype in schema.items()
        }

    def with_row_index(self, name: str, order_by: Sequence[str] | None) -> Self:
        frame = self.native
        if order_by is None:
            result = frame.with_row_index(name)
        else:
            end = pl.count() if self._backend_version < (0, 20, 5) else pl.len()
            result = frame.select(
                pl.int_range(start=0, end=end).sort_by(order_by).alias(name), pl.all()
            )

        return self._with_native(result)


class PolarsDataFrame(PolarsBaseFrame[pl.DataFrame]):
    clone: Method[Self]
    collect: Method[CompliantDataFrameAny]
    estimated_size: Method[int | float]
    gather_every: Method[Self]
    item: Method[Any]
    iter_rows: Method[Iterator[tuple[Any, ...]] | Iterator[Mapping[str, Any]]]
    is_unique: Method[PolarsSeries]
    row: Method[tuple[Any, ...]]
    rows: Method[Sequence[tuple[Any, ...]] | Sequence[Mapping[str, Any]]]
    sample: Method[Self]
    to_arrow: Method[pa.Table]
    to_pandas: Method[pd.DataFrame]
    # NOTE: `write_csv` requires an `@overload` for `str | None`
    # Can't do that here 😟
    write_csv: Method[Any]
    write_parquet: Method[None]

    @classmethod
    def from_arrow(cls, data: IntoArrowTable, /, *, context: _LimitedContext) -> Self:
        if context._implementation._backend_version() >= (1, 3):
            native = pl.DataFrame(data)
        else:  # pragma: no cover
            native = cast("pl.DataFrame", pl.from_arrow(_into_arrow_table(data, context)))
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

        pl_schema = Schema(schema).to_polars() if schema is not None else schema
        return cls.from_native(pl.from_dict(data, pl_schema), context=context)

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

        pl_schema = Schema(schema).to_polars() if schema is not None else schema
        if not data:
            native = pl.DataFrame(schema=pl_schema)
        elif FROM_DICTS_ACCEPTS_MAPPINGS or isinstance(data[0], dict):
            native = pl.from_dicts(data, pl_schema)
        else:  # pragma: no cover
            columns = pl_schema or tuple(data[0])
            native = pl.DataFrame(
                (tuple(row.values()) for row in data), schema=columns, orient="row"
            )

        return cls.from_native(native, context=context)

    @staticmethod
    def _is_native(obj: pl.DataFrame | Any) -> TypeIs[pl.DataFrame]:
        return isinstance(obj, pl.DataFrame)

    @classmethod
    def from_numpy(
        cls,
        data: _2DArray,
        /,
        *,
        context: _LimitedContext,  # NOTE: Maybe only `Implementation`?
        schema: IntoSchema | Sequence[str] | None,
    ) -> Self:
        from narwhals.schema import Schema

        pl_schema = (
            Schema(schema).to_polars()
            if isinstance(schema, (Mapping, Schema))
            else schema
        )
        return cls.from_native(pl.from_numpy(data, pl_schema), context=context)

    def to_narwhals(self) -> DataFrame[pl.DataFrame]:
        return self._version.dataframe(self, level="full")

    def __repr__(self) -> str:  # pragma: no cover
        return "PolarsDataFrame"

    def __narwhals_dataframe__(self) -> Self:
        return self

    @overload
    def _from_native_object(self, obj: pl.Series) -> PolarsSeries: ...

    @overload
    def _from_native_object(self, obj: pl.DataFrame) -> Self: ...

    @overload
    def _from_native_object(self, obj: T) -> T: ...

    def _from_native_object(
        self, obj: pl.Series | pl.DataFrame | T
    ) -> Self | PolarsSeries | T:
        if isinstance(obj, pl.Series):
            return PolarsSeries.from_native(obj, context=self)
        if self._is_native(obj):
            return self._with_native(obj)
        # scalar
        return obj

    def __len__(self) -> int:
        return len(self.native)

    def __getattr__(self, attr: str) -> Any:
        if attr not in INHERITED_METHODS:  # pragma: no cover
            msg = f"{self.__class__.__name__} has not attribute '{attr}'."
            raise AttributeError(msg)

        def func(*args: Any, **kwargs: Any) -> Any:
            pos, kwds = extract_args_kwargs(args, kwargs)
            try:
                return self._from_native_object(getattr(self.native, attr)(*pos, **kwds))
            except pl.exceptions.ColumnNotFoundError as e:  # pragma: no cover
                msg = f"{e!s}\n\nHint: Did you mean one of these columns: {self.columns}?"
                raise ColumnNotFoundError(msg) from e
            except Exception as e:  # noqa: BLE001
                raise catch_polars_exception(e) from None

        return func

    def __array__(
        self, dtype: Any | None = None, *, copy: bool | None = None
    ) -> _2DArray:
        if self._backend_version < (0, 20, 28) and copy is not None:
            msg = "`copy` in `__array__` is only supported for 'polars>=0.20.28'"
            raise NotImplementedError(msg)
        if self._backend_version < (0, 20, 28):
            return self.native.__array__(dtype)
        return self.native.__array__(dtype)

    def to_numpy(self, dtype: Any = None, *, copy: bool | None = None) -> _2DArray:
        return self.native.to_numpy()

    @property
    def shape(self) -> tuple[int, int]:
        return self.native.shape

    def __getitem__(  # noqa: C901, PLR0912
        self,
        item: tuple[
            SingleIndexSelector | MultiIndexSelector[PolarsSeries],
            MultiColSelector[PolarsSeries],
        ],
    ) -> Any:
        rows, columns = item
        if self._backend_version > (0, 20, 30):
            rows_native = rows.native if is_compliant_series(rows) else rows
            columns_native = columns.native if is_compliant_series(columns) else columns
            selector = rows_native, columns_native
            selected = self.native.__getitem__(selector)  # type: ignore[index]
            return self._from_native_object(selected)
        else:  # pragma: no cover # noqa: RET505
            # TODO(marco): we can delete this branch after Polars==0.20.30 becomes the minimum
            # Polars version we support
            # This mostly mirrors the logic in `EagerDataFrame.__getitem__`.
            rows = list(rows) if isinstance(rows, tuple) else rows
            columns = list(columns) if isinstance(columns, tuple) else columns
            if is_numpy_array_1d(columns):
                columns = columns.tolist()

            native = self.native
            if not is_slice_none(columns):
                if isinstance(columns, Sized) and len(columns) == 0:
                    return self.select()
                if is_index_selector(columns):
                    if is_slice_index(columns) or is_range(columns):
                        native = native.select(
                            self.columns[slice(columns.start, columns.stop, columns.step)]
                        )
                    # NOTE: `mypy` loses track of `PolarsSeries` when `is_compliant_series` is used here
                    # `pyright` is fine
                    elif isinstance(columns, PolarsSeries):
                        native = native[:, columns.native.to_list()]
                    else:
                        native = native[:, columns]
                elif isinstance(columns, slice):
                    native = native.select(
                        self.columns[
                            slice(*convert_str_slice_to_int_slice(columns, self.columns))
                        ]
                    )
                elif is_compliant_series(columns):
                    native = native.select(columns.native.to_list())
                elif is_sequence_like(columns):
                    native = native.select(columns)
                else:
                    msg = f"Unreachable code, got unexpected type: {type(columns)}"
                    raise AssertionError(msg)

            if not is_slice_none(rows):
                if isinstance(rows, int):
                    native = native[[rows], :]
                elif isinstance(rows, (slice, range)):
                    native = native[rows, :]
                elif is_compliant_series(rows):
                    native = native[rows.native, :]
                elif is_sequence_like(rows):
                    native = native[rows, :]
                else:
                    msg = f"Unreachable code, got unexpected type: {type(rows)}"
                    raise AssertionError(msg)

            return self._with_native(native)

    def get_column(self, name: str) -> PolarsSeries:
        return PolarsSeries.from_native(self.native.get_column(name), context=self)

    def iter_columns(self) -> Iterator[PolarsSeries]:
        for series in self.native.iter_columns():
            yield PolarsSeries.from_native(series, context=self)

    def lazy(
        self,
        backend: _LazyAllowedImpl | None = None,
        *,
        session: SparkSession | None = None,
    ) -> CompliantLazyFrameAny:
        if backend is None or backend is Implementation.POLARS:
            return PolarsLazyFrame.from_native(self.native.lazy(), context=self)
        if backend is Implementation.DUCKDB:
            import duckdb  # ignore-banned-import

            from narwhals._duckdb.dataframe import DuckDBLazyFrame

            _df = self.native
            return DuckDBLazyFrame(
                duckdb.table("_df"), validate_backend_version=True, version=self._version
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
                self,  # pyright: ignore[reportArgumentType]
                session=session,
                implementation=backend,
                version=self._version,
            )

        raise AssertionError  # pragma: no cover

    @overload
    def to_dict(self, *, as_series: Literal[True]) -> dict[str, PolarsSeries]: ...

    @overload
    def to_dict(self, *, as_series: Literal[False]) -> dict[str, list[Any]]: ...

    def to_dict(
        self, *, as_series: bool
    ) -> dict[str, PolarsSeries] | dict[str, list[Any]]:
        if as_series:
            return {
                name: PolarsSeries.from_native(col, context=self)
                for name, col in self.native.to_dict().items()
            }
        return self.native.to_dict(as_series=False)

    def group_by(
        self, keys: Sequence[str] | Sequence[PolarsExpr], *, drop_null_keys: bool
    ) -> PolarsGroupBy:
        from narwhals._polars.group_by import PolarsGroupBy

        return PolarsGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        to_drop = parse_columns_to_drop(self, columns, strict=strict)
        return self._with_native(self.native.drop(to_drop))

    @requires.backend_version((1,))
    def pivot(
        self,
        on: Sequence[str],
        *,
        index: Sequence[str] | None,
        values: Sequence[str] | None,
        aggregate_function: PivotAgg | None,
        sort_columns: bool,
        separator: str,
    ) -> Self:
        try:
            result = self.native.pivot(
                on,
                index=index,
                values=values,
                aggregate_function=aggregate_function,
                sort_columns=sort_columns,
                separator=separator,
            )
        except Exception as e:  # noqa: BLE001
            raise catch_polars_exception(e) from None
        return self._from_native_object(result)

    def to_polars(self) -> pl.DataFrame:
        return self.native

    def join(
        self,
        other: PolarsBaseFrame[pl.DataFrame],
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        try:
            return super().join(
                other=other, how=how, left_on=left_on, right_on=right_on, suffix=suffix
            )
        except Exception as e:  # noqa: BLE001
            raise catch_polars_exception(e) from None

    def top_k(
        self, k: int, *, by: str | Iterable[str], reverse: bool | Sequence[bool]
    ) -> Self:
        try:
            return super().top_k(k=k, by=by, reverse=reverse)
        except Exception as e:  # noqa: BLE001  # pragma: no cover
            raise catch_polars_exception(e) from None


class PolarsLazyFrame(PolarsBaseFrame[pl.LazyFrame]):
    sink_parquet: Method[None]

    @staticmethod
    def _is_native(obj: pl.LazyFrame | Any) -> TypeIs[pl.LazyFrame]:
        return isinstance(obj, pl.LazyFrame)

    def to_narwhals(self) -> LazyFrame[pl.LazyFrame]:
        return self._version.lazyframe(self, level="lazy")

    def __repr__(self) -> str:  # pragma: no cover
        return "PolarsLazyFrame"

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def __getattr__(self, attr: str) -> Any:
        if attr not in INHERITED_METHODS:  # pragma: no cover
            msg = f"{self.__class__.__name__} has not attribute '{attr}'."
            raise AttributeError(msg)

        def func(*args: Any, **kwargs: Any) -> Any:
            pos, kwds = extract_args_kwargs(args, kwargs)
            try:
                return self._with_native(getattr(self.native, attr)(*pos, **kwds))
            except pl.exceptions.ColumnNotFoundError as e:  # pragma: no cover
                raise ColumnNotFoundError(str(e)) from e

        return func

    def _iter_columns(self) -> Iterator[PolarsSeries]:  # pragma: no cover
        yield from self.collect(Implementation.POLARS).iter_columns()

    def collect_schema(self) -> dict[str, DType]:
        try:
            return super().collect_schema()
        except Exception as e:  # noqa: BLE001
            raise catch_polars_exception(e) from None

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        try:
            result = self.native.collect(**kwargs)
        except Exception as e:  # noqa: BLE001
            raise catch_polars_exception(e) from None

        if backend is None or backend is Implementation.POLARS:
            return PolarsDataFrame.from_native(result, context=self)

        if backend is Implementation.PANDAS:
            from narwhals._pandas_like.dataframe import PandasLikeDataFrame

            return PandasLikeDataFrame(
                result.to_pandas(),
                implementation=Implementation.PANDAS,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=False,
            )

        if backend is Implementation.PYARROW:
            from narwhals._arrow.dataframe import ArrowDataFrame

            return ArrowDataFrame(
                result.to_arrow(),
                validate_backend_version=True,
                version=self._version,
                validate_column_names=False,
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise ValueError(msg)  # pragma: no cover

    def group_by(
        self, keys: Sequence[str] | Sequence[PolarsExpr], *, drop_null_keys: bool
    ) -> PolarsLazyGroupBy:
        from narwhals._polars.group_by import PolarsLazyGroupBy

        return PolarsLazyGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        if self._backend_version < (1, 0, 0):
            return self._with_native(self.native.drop(columns))
        return self._with_native(self.native.drop(columns, strict=strict))
