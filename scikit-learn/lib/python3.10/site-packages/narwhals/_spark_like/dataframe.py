from __future__ import annotations

from functools import reduce
from operator import and_
from typing import TYPE_CHECKING, Any

from narwhals._exceptions import issue_warning
from narwhals._native import is_native_spark_like
from narwhals._spark_like.utils import (
    catch_pyspark_connect_exception,
    catch_pyspark_sql_exception,
    evaluate_exprs,
    import_functions,
    import_native_dtypes,
    import_window,
    native_to_narwhals_dtype,
)
from narwhals._sql.dataframe import SQLLazyFrame
from narwhals._utils import (
    Implementation,
    ValidateBackendVersion,
    extend_bool,
    generate_temporary_column_name,
    not_implemented,
    parse_columns_to_drop,
    to_pyarrow_table,
    zip_strict,
)
from narwhals.exceptions import InvalidOperationError

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import pyarrow as pa
    from sqlframe.base.column import Column
    from sqlframe.base.dataframe import BaseDataFrame
    from sqlframe.base.window import Window
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._compliant.typing import CompliantDataFrameAny
    from narwhals._spark_like.expr import SparkLikeExpr
    from narwhals._spark_like.group_by import SparkLikeLazyGroupBy
    from narwhals._spark_like.namespace import SparkLikeNamespace
    from narwhals._spark_like.utils import SparkSession
    from narwhals._typing import _EagerAllowedImpl
    from narwhals._utils import Version, _LimitedContext
    from narwhals.dataframe import LazyFrame
    from narwhals.dtypes import DType
    from narwhals.typing import JoinStrategy, UniqueKeepStrategy

    SQLFrameDataFrame = BaseDataFrame[Any, Any, Any, Any, Any]

Incomplete: TypeAlias = Any  # pragma: no cover
"""Marker for working code that fails type checking."""


class SparkLikeLazyFrame(
    SQLLazyFrame["SparkLikeExpr", "SQLFrameDataFrame", "LazyFrame[SQLFrameDataFrame]"],
    ValidateBackendVersion,
):
    def __init__(
        self,
        native_dataframe: SQLFrameDataFrame,
        *,
        version: Version,
        implementation: Implementation,
        validate_backend_version: bool = False,
    ) -> None:
        self._native_frame: SQLFrameDataFrame = native_dataframe
        self._implementation = implementation
        self._version = version
        self._cached_schema: dict[str, DType] | None = None
        self._cached_columns: list[str] | None = None
        if validate_backend_version:  # pragma: no cover
            self._validate_backend_version()

    @property
    def _backend_version(self) -> tuple[int, ...]:  # pragma: no cover
        return self._implementation._backend_version()

    @property
    def _F(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            from sqlframe.base import functions

            return functions
        return import_functions(self._implementation)

    @property
    def _native_dtypes(self):  # type: ignore[no-untyped-def] # noqa: ANN202
        if TYPE_CHECKING:
            from sqlframe.base import types

            return types
        return import_native_dtypes(self._implementation)

    @property
    def _Window(self) -> type[Window]:
        if TYPE_CHECKING:
            from sqlframe.base.window import Window

            return Window
        return import_window(self._implementation)

    @staticmethod
    def _is_native(obj: SQLFrameDataFrame | Any) -> TypeIs[SQLFrameDataFrame]:
        return is_native_spark_like(obj)

    @classmethod
    def from_native(cls, data: SQLFrameDataFrame, /, *, context: _LimitedContext) -> Self:
        return cls(data, version=context._version, implementation=context._implementation)

    def to_narwhals(self) -> LazyFrame[SQLFrameDataFrame]:
        return self._version.lazyframe(self, level="lazy")

    def __native_namespace__(self) -> ModuleType:  # pragma: no cover
        return self._implementation.to_native_namespace()

    def __narwhals_namespace__(self) -> SparkLikeNamespace:
        from narwhals._spark_like.namespace import SparkLikeNamespace

        return SparkLikeNamespace(
            version=self._version, implementation=self._implementation
        )

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def _with_version(self, version: Version) -> Self:
        return self.__class__(
            self.native, version=version, implementation=self._implementation
        )

    def _with_native(self, df: SQLFrameDataFrame) -> Self:
        return self.__class__(
            df, version=self._version, implementation=self._implementation
        )

    def _to_arrow_schema(self) -> pa.Schema:  # pragma: no cover
        import pyarrow as pa  # ignore-banned-import

        from narwhals._arrow.utils import narwhals_to_native_dtype

        schema: list[tuple[str, pa.DataType]] = []
        nw_schema = self.collect_schema()
        native_schema = self.native.schema
        for key, value in nw_schema.items():
            try:
                native_dtype = narwhals_to_native_dtype(value, self._version)
            except Exception as exc:  # noqa: BLE001,PERF203
                native_spark_dtype = native_schema[key].dataType  # type: ignore[index]
                # If we can't convert the type, just set it to `pa.null`, and warn.
                # Avoid the warning if we're starting from PySpark's void type.
                # We can avoid the check when we introduce `nw.Null` dtype.
                null_type = self._native_dtypes.NullType  # pyright: ignore[reportAttributeAccessIssue]
                if not isinstance(native_spark_dtype, null_type):
                    issue_warning(
                        f"Could not convert dtype {native_spark_dtype} to PyArrow dtype, {exc!r}",
                        UserWarning,
                    )
                schema.append((key, pa.null()))
            else:
                schema.append((key, native_dtype))
        return pa.schema(schema)

    def _collect_to_arrow(self) -> pa.Table:
        if self._implementation.is_pyspark() and self._backend_version < (4,):
            import pyarrow as pa  # ignore-banned-import

            try:
                return pa.Table.from_batches(self.native._collect_as_arrow())
            except ValueError as exc:
                if "at least one RecordBatch" in str(exc):
                    # Empty dataframe

                    data: dict[str, list[Any]] = {k: [] for k in self.columns}
                    pa_schema = self._to_arrow_schema()
                    return pa.Table.from_pydict(data, schema=pa_schema)
                raise  # pragma: no cover
        elif self._implementation.is_pyspark_connect() and self._backend_version < (4,):
            import pyarrow as pa  # ignore-banned-import

            pa_schema = self._to_arrow_schema()
            return pa.Table.from_pandas(self.native.toPandas(), schema=pa_schema)
        else:
            return to_pyarrow_table(self.native.toArrow())

    def _iter_columns(self) -> Iterator[Column]:
        for col in self.columns:
            yield self._F.col(col)

    @property
    def columns(self) -> list[str]:
        if self._cached_columns is None:
            self._cached_columns = (
                list(self.schema)
                if self._cached_schema is not None
                else self.native.columns
            )
        return self._cached_columns

    def _collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        if backend is Implementation.PANDAS:
            from narwhals._pandas_like.dataframe import PandasLikeDataFrame

            return PandasLikeDataFrame(
                self.native.toPandas(),
                implementation=Implementation.PANDAS,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is None or backend is Implementation.PYARROW:
            from narwhals._arrow.dataframe import ArrowDataFrame

            return ArrowDataFrame(
                self._collect_to_arrow(),
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is Implementation.POLARS:
            import polars as pl  # ignore-banned-import

            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame(
                pl.from_arrow(self._collect_to_arrow()),  # type: ignore[arg-type]
                validate_backend_version=True,
                version=self._version,
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise ValueError(msg)  # pragma: no cover

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        if self._implementation.is_pyspark_connect():
            try:
                return self._collect(backend, **kwargs)
            except Exception as e:  # noqa: BLE001
                raise catch_pyspark_connect_exception(e) from None
        return self._collect(backend, **kwargs)

    def simple_select(self, *column_names: str) -> Self:
        return self._with_native(self.native.select(*column_names))

    def aggregate(self, *exprs: SparkLikeExpr) -> Self:
        new_columns = evaluate_exprs(self, *exprs)

        new_columns_list = [col.alias(col_name) for col_name, col in new_columns]
        if self._implementation.is_pyspark():
            try:
                return self._with_native(self.native.agg(*new_columns_list))
            except Exception as e:  # noqa: BLE001
                raise catch_pyspark_sql_exception(e, self) from None
        return self._with_native(self.native.agg(*new_columns_list))

    def select(self, *exprs: SparkLikeExpr) -> Self:
        new_columns = evaluate_exprs(self, *exprs)
        new_columns_list = [col.alias(col_name) for (col_name, col) in new_columns]
        if self._implementation.is_pyspark():  # pragma: no cover
            try:
                return self._with_native(self.native.select(*new_columns_list))
            except Exception as e:  # noqa: BLE001
                raise catch_pyspark_sql_exception(e, self) from None
        return self._with_native(self.native.select(*new_columns_list))

    def with_columns(self, *exprs: SparkLikeExpr) -> Self:
        new_columns = evaluate_exprs(self, *exprs)
        if self._implementation.is_pyspark():  # pragma: no cover
            try:
                return self._with_native(self.native.withColumns(dict(new_columns)))
            except Exception as e:  # noqa: BLE001
                raise catch_pyspark_sql_exception(e, self) from None

        return self._with_native(self.native.withColumns(dict(new_columns)))

    def filter(self, predicate: SparkLikeExpr) -> Self:
        # `[0]` is safe as the predicate's expression only returns a single column
        condition = predicate._call(self)[0]
        if self._implementation.is_pyspark():
            try:
                return self._with_native(self.native.where(condition))
            except Exception as e:  # noqa: BLE001
                raise catch_pyspark_sql_exception(e, self) from None
        return self._with_native(self.native.where(condition))

    @property
    def schema(self) -> dict[str, DType]:
        if self._cached_schema is None:
            self._cached_schema = {
                field.name: native_to_narwhals_dtype(
                    field.dataType,
                    self._version,
                    self._native_dtypes,
                    self.native.sparkSession,
                )
                for field in self.native.schema
            }
        return self._cached_schema

    def collect_schema(self) -> dict[str, DType]:
        return self.schema

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        columns_to_drop = parse_columns_to_drop(self, columns, strict=strict)
        return self._with_native(self.native.drop(*columns_to_drop))

    def head(self, n: int) -> Self:
        return self._with_native(self.native.limit(n))

    def group_by(
        self, keys: Sequence[str] | Sequence[SparkLikeExpr], *, drop_null_keys: bool
    ) -> SparkLikeLazyGroupBy:
        from narwhals._spark_like.group_by import SparkLikeLazyGroupBy

        return SparkLikeLazyGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def sort(self, *by: str, descending: bool | Sequence[bool], nulls_last: bool) -> Self:
        descending = extend_bool(descending, len(by))
        if nulls_last:
            sort_funcs = (
                self._F.desc_nulls_last if d else self._F.asc_nulls_last
                for d in descending
            )
        else:
            sort_funcs = (
                self._F.desc_nulls_first if d else self._F.asc_nulls_first
                for d in descending
            )

        sort_cols = [sort_f(col) for col, sort_f in zip_strict(by, sort_funcs)]
        return self._with_native(self.native.sort(*sort_cols))

    def top_k(self, k: int, *, by: Iterable[str], reverse: bool | Sequence[bool]) -> Self:
        by = tuple(by)
        reverse = extend_bool(reverse, len(by))
        sort_funcs = (
            self._F.desc_nulls_last if not d else self._F.asc_nulls_last for d in reverse
        )
        sort_cols = [sort_f(col) for col, sort_f in zip_strict(by, sort_funcs)]
        return self._with_native(self.native.sort(*sort_cols).limit(k))

    def drop_nulls(self, subset: Sequence[str] | None) -> Self:
        subset = list(subset) if subset else None
        return self._with_native(self.native.dropna(subset=subset))

    def rename(self, mapping: Mapping[str, str]) -> Self:
        rename_mapping = {
            colname: mapping.get(colname, colname) for colname in self.columns
        }
        return self._with_native(
            self.native.select(
                [self._F.col(old).alias(new) for old, new in rename_mapping.items()]
            )
        )

    def unique(
        self,
        subset: Sequence[str] | None,
        *,
        keep: UniqueKeepStrategy,
        order_by: Sequence[str] | None,
    ) -> Self:
        subset_ = subset or self.columns
        if error := self._check_columns_exist(subset_):
            raise error
        tmp_name = generate_temporary_column_name(8, self.columns, prefix="row_index_")
        window = self._Window.partitionBy(subset_)
        if order_by and keep == "last":
            window = window.orderBy(*[self._F.desc_nulls_last(x) for x in order_by])
        elif order_by:
            window = window.orderBy(*[self._F.asc_nulls_first(x) for x in order_by])
        else:
            window = window.orderBy(self._F.lit(1))
        if keep == "none":
            expr = self._F.count("*").over(window)
        else:
            expr = self._F.row_number().over(window)
        df = (
            self.native.withColumn(tmp_name, expr)
            .filter(self._F.col(tmp_name) == self._F.lit(1))
            .drop(tmp_name)
        )
        return self._with_native(df)

    def join(
        self,
        other: Self,
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        left_columns = self.columns
        right_columns = other.columns

        right_on_: list[str] = list(right_on) if right_on is not None else []
        left_on_: list[str] = list(left_on) if left_on is not None else []

        # create a mapping for columns on other
        # `right_on` columns will be renamed as `left_on`
        # the remaining columns will be either added the suffix or left unchanged.
        right_cols_to_rename = (
            [c for c in right_columns if c not in right_on_]
            if how != "full"
            else right_columns
        )

        rename_mapping = {
            **dict(zip(right_on_, left_on_)),
            **{
                colname: f"{colname}{suffix}" if colname in left_columns else colname
                for colname in right_cols_to_rename
            },
        }
        other_native = other.native.select(
            [self._F.col(old).alias(new) for old, new in rename_mapping.items()]
        )

        # If how in {"semi", "anti"}, then resulting columns are same as left columns
        # Otherwise, we add the right columns with the new mapping, while keeping the
        # original order of right_columns.
        col_order = left_columns.copy()

        if how in {"inner", "left", "cross"}:
            col_order.extend(
                rename_mapping[colname]
                for colname in right_columns
                if colname not in right_on_
            )
        elif how == "full":
            col_order.extend(rename_mapping.values())

        right_on_remapped = [rename_mapping[c] for c in right_on_]
        on_ = (
            reduce(
                and_,
                (
                    getattr(self.native, left_key) == getattr(other_native, right_key)
                    for left_key, right_key in zip_strict(left_on_, right_on_remapped)
                ),
            )
            if how == "full"
            else None
            if how == "cross"
            else left_on_
        )
        how_native = "full_outer" if how == "full" else how
        return self._with_native(
            self.native.join(other_native, on=on_, how=how_native).select(col_order)
        )

    def explode(self, columns: Sequence[str]) -> Self:
        dtypes = self._version.dtypes

        schema = self.collect_schema()
        for col_to_explode in columns:
            dtype = schema[col_to_explode]

            if dtype != dtypes.List:
                msg = (
                    f"`explode` operation not supported for dtype `{dtype}`, "
                    "expected List type"
                )
                raise InvalidOperationError(msg)

        column_names = self.columns

        if len(columns) != 1:
            msg = (
                "Exploding on multiple columns is not supported with SparkLike backend since "
                "we cannot guarantee that the exploded columns have matching element counts."
            )
            raise NotImplementedError(msg)

        if self._implementation.is_pyspark() or self._implementation.is_pyspark_connect():
            return self._with_native(
                self.native.select(
                    *[
                        self._F.col(col_name).alias(col_name)
                        if col_name != columns[0]
                        else self._F.explode_outer(col_name).alias(col_name)
                        for col_name in column_names
                    ]
                )
            )
        if self._implementation.is_sqlframe():
            # Not every sqlframe dialect supports `explode_outer` function
            # (see https://github.com/eakmanrq/sqlframe/blob/3cb899c515b101ff4c197d84b34fae490d0ed257/sqlframe/base/functions.py#L2288-L2289)
            # therefore we simply explode the array column which will ignore nulls and
            # zero sized arrays, and append these specific condition with nulls (to
            # match polars behavior).

            def null_condition(col_name: str) -> Column:
                return self._F.isnull(col_name) | (self._F.array_size(col_name) == 0)

            return self._with_native(
                self.native.select(
                    *[
                        self._F.col(col_name).alias(col_name)
                        if col_name != columns[0]
                        else self._F.explode(col_name).alias(col_name)
                        for col_name in column_names
                    ]
                ).union(
                    self.native.filter(null_condition(columns[0])).select(
                        *[
                            self._F.col(col_name).alias(col_name)
                            if col_name != columns[0]
                            else self._F.lit(None).alias(col_name)
                            for col_name in column_names
                        ]
                    )
                )
            )
        msg = "Unreachable code, please report an issue at https://github.com/narwhals-dev/narwhals/issues"  # pragma: no cover
        raise AssertionError(msg)

    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        if self._implementation.is_sqlframe():
            if variable_name == "":
                msg = "`variable_name` cannot be empty string for sqlframe backend."
                raise NotImplementedError(msg)

            if value_name == "":
                msg = "`value_name` cannot be empty string for sqlframe backend."
                raise NotImplementedError(msg)
        else:  # pragma: no cover
            pass

        ids = tuple(index) if index else ()
        values = (
            tuple(set(self.columns).difference(set(ids))) if on is None else tuple(on)
        )
        unpivoted_native_frame = self.native.unpivot(
            ids=ids,
            values=values,
            variableColumnName=variable_name,
            valueColumnName=value_name,
        )
        if index is None:
            unpivoted_native_frame = unpivoted_native_frame.drop(*ids)
        return self._with_native(unpivoted_native_frame)

    def with_row_index(self, name: str, order_by: Sequence[str]) -> Self:
        if order_by is None:
            msg = "Cannot pass `order_by` to `with_row_index` for PySpark-like"
            raise TypeError(msg)
        row_index_expr = (
            self._F.row_number().over(
                self._Window.partitionBy(self._F.lit(1)).orderBy(*order_by)
            )
            - 1
        ).alias(name)
        return self._with_native(self.native.select(row_index_expr, *self.columns))

    def sink_parquet(self, file: str | Path | BytesIO) -> None:
        self.native.write.parquet(file)

    @classmethod
    def _from_compliant_dataframe(
        cls,
        frame: CompliantDataFrameAny,
        /,
        *,
        session: SparkSession,
        implementation: Implementation,
        version: Version,
    ) -> SparkLikeLazyFrame:
        from importlib.util import find_spec

        impl = implementation
        is_spark_v4 = (not impl.is_sqlframe()) and impl._backend_version() >= (4, 0, 0)
        if is_spark_v4:  # pragma: no cover
            # pyspark.sql requires pyarrow to be installed from v4.0.0
            # and since v4.0.0 the input to `createDataFrame` can be a PyArrow Table.
            data: Any = frame.to_arrow()
        elif find_spec("pandas"):
            data = frame.to_pandas()
        else:  # pragma: no cover
            data = tuple(frame.iter_rows(named=True, buffer_size=512))

        return cls(
            session.createDataFrame(data),
            version=version,
            implementation=implementation,
            validate_backend_version=True,
        )

    gather_every = not_implemented.deprecated(
        "`LazyFrame.gather_every` is deprecated and will be removed in a future version."
    )
    join_asof = not_implemented()
    tail = not_implemented.deprecated(
        "`LazyFrame.tail` is deprecated and will be removed in a future version."
    )
