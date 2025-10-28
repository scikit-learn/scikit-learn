from __future__ import annotations

from functools import reduce
from operator import and_
from typing import TYPE_CHECKING, Any

import duckdb
from duckdb import StarExpression

from narwhals._duckdb.utils import (
    DeferredTimeZone,
    F,
    catch_duckdb_exception,
    col,
    evaluate_exprs,
    join_column_names,
    lit,
    native_to_narwhals_dtype,
    window_expression,
)
from narwhals._sql.dataframe import SQLLazyFrame
from narwhals._utils import (
    Implementation,
    ValidateBackendVersion,
    Version,
    extend_bool,
    generate_temporary_column_name,
    not_implemented,
    parse_columns_to_drop,
    requires,
    zip_strict,
)
from narwhals.dependencies import get_duckdb
from narwhals.exceptions import InvalidOperationError

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from io import BytesIO
    from pathlib import Path
    from types import ModuleType

    import pandas as pd
    import pyarrow as pa
    from duckdb import Expression
    from typing_extensions import Self, TypeIs

    from narwhals._compliant.typing import CompliantDataFrameAny
    from narwhals._duckdb.expr import DuckDBExpr
    from narwhals._duckdb.group_by import DuckDBGroupBy
    from narwhals._duckdb.namespace import DuckDBNamespace
    from narwhals._duckdb.series import DuckDBInterchangeSeries
    from narwhals._duckdb.utils import duckdb_dtypes
    from narwhals._typing import _EagerAllowedImpl
    from narwhals._utils import _LimitedContext
    from narwhals.dataframe import LazyFrame
    from narwhals.dtypes import DType
    from narwhals.stable.v1 import DataFrame as DataFrameV1
    from narwhals.typing import AsofJoinStrategy, JoinStrategy, UniqueKeepStrategy


class DuckDBLazyFrame(
    SQLLazyFrame[
        "DuckDBExpr",
        "duckdb.DuckDBPyRelation",
        "LazyFrame[duckdb.DuckDBPyRelation] | DataFrameV1[duckdb.DuckDBPyRelation]",
    ],
    ValidateBackendVersion,
):
    _implementation = Implementation.DUCKDB

    def __init__(
        self,
        df: duckdb.DuckDBPyRelation,
        *,
        version: Version,
        validate_backend_version: bool = False,
    ) -> None:
        self._native_frame: duckdb.DuckDBPyRelation = df
        self._version = version
        self._cached_native_schema: dict[str, duckdb_dtypes.DuckDBPyType] | None = None
        self._cached_columns: list[str] | None = None
        if validate_backend_version:
            self._validate_backend_version()

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    @staticmethod
    def _is_native(obj: duckdb.DuckDBPyRelation | Any) -> TypeIs[duckdb.DuckDBPyRelation]:
        return isinstance(obj, duckdb.DuckDBPyRelation)

    @classmethod
    def from_native(
        cls, data: duckdb.DuckDBPyRelation, /, *, context: _LimitedContext
    ) -> Self:
        return cls(data, version=context._version)

    def to_narwhals(
        self, *args: Any, **kwds: Any
    ) -> LazyFrame[duckdb.DuckDBPyRelation] | DataFrameV1[duckdb.DuckDBPyRelation]:
        if self._version is Version.V1:
            from narwhals.stable.v1 import DataFrame as DataFrameV1

            return DataFrameV1(self, level="interchange")  # type: ignore[no-any-return]
        return self._version.lazyframe(self, level="lazy")

    def __narwhals_dataframe__(self) -> Self:  # pragma: no cover
        # Keep around for backcompat.
        if self._version is not Version.V1:
            msg = "__narwhals_dataframe__ is not implemented for DuckDBLazyFrame"
            raise AttributeError(msg)
        return self

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def __native_namespace__(self) -> ModuleType:
        return get_duckdb()  # type: ignore[no-any-return]

    def __narwhals_namespace__(self) -> DuckDBNamespace:
        from narwhals._duckdb.namespace import DuckDBNamespace

        return DuckDBNamespace(version=self._version)

    def get_column(self, name: str) -> DuckDBInterchangeSeries:
        from narwhals._duckdb.series import DuckDBInterchangeSeries

        return DuckDBInterchangeSeries(self.native.select(name), version=self._version)

    def _iter_columns(self) -> Iterator[Expression]:
        for name in self.columns:
            yield col(name)

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        if backend is None or backend is Implementation.PYARROW:
            from narwhals._arrow.dataframe import ArrowDataFrame

            if self._backend_version < (1, 4):
                ret = self.native.arrow()
            else:  # pragma: no cover
                ret = self.native.fetch_arrow_table()
            return ArrowDataFrame(
                ret,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is Implementation.PANDAS:
            from narwhals._pandas_like.dataframe import PandasLikeDataFrame

            return PandasLikeDataFrame(
                self.native.df(),
                implementation=Implementation.PANDAS,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is Implementation.POLARS:
            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame(
                self.native.pl(), validate_backend_version=True, version=self._version
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise ValueError(msg)  # pragma: no cover

    def head(self, n: int) -> Self:
        return self._with_native(self.native.limit(n))

    def simple_select(self, *column_names: str) -> Self:
        return self._with_native(self.native.select(*column_names))

    def aggregate(self, *exprs: DuckDBExpr) -> Self:
        selection = [val.alias(name) for name, val in evaluate_exprs(self, *exprs)]
        try:
            return self._with_native(self.native.aggregate(selection))  # type: ignore[arg-type]
        except Exception as e:  # noqa: BLE001
            raise catch_duckdb_exception(e, self) from None

    def select(self, *exprs: DuckDBExpr) -> Self:
        selection = (val.alias(name) for name, val in evaluate_exprs(self, *exprs))
        try:
            return self._with_native(self.native.select(*selection))
        except Exception as e:  # noqa: BLE001
            raise catch_duckdb_exception(e, self) from None

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        columns_to_drop = parse_columns_to_drop(self, columns, strict=strict)
        selection = [col(name) for name in self.columns if name not in columns_to_drop]
        return self._with_native(self.native.select(*selection))

    def lazy(self, backend: None = None, **_: None) -> Self:
        # The `backend`` argument has no effect but we keep it here for
        # backwards compatibility because in `narwhals.stable.v1`
        # function `.from_native()` will return a DataFrame for DuckDB.

        if backend is not None:  # pragma: no cover
            msg = "`backend` argument is not supported for DuckDB"
            raise ValueError(msg)
        return self

    def with_columns(self, *exprs: DuckDBExpr) -> Self:
        new_columns_map = dict(evaluate_exprs(self, *exprs))
        result = [
            new_columns_map.pop(name).alias(name)
            if name in new_columns_map
            else col(name)
            for name in self.columns
        ]
        result.extend(value.alias(name) for name, value in new_columns_map.items())
        try:
            return self._with_native(self.native.select(*result))
        except Exception as e:  # noqa: BLE001
            raise catch_duckdb_exception(e, self) from None

    def filter(self, predicate: DuckDBExpr) -> Self:
        # `[0]` is safe as the predicate's expression only returns a single column
        mask = predicate(self)[0]
        try:
            return self._with_native(self.native.filter(mask))
        except Exception as e:  # noqa: BLE001
            raise catch_duckdb_exception(e, self) from None

    @property
    def schema(self) -> dict[str, DType]:
        if self._cached_native_schema is None:
            # Note: prefer `self._cached_native_schema` over `functools.cached_property`
            # due to Python3.13 failures.
            self._cached_native_schema = dict(zip(self.columns, self.native.types))

        deferred_time_zone = DeferredTimeZone(self.native)
        return {
            column_name: native_to_narwhals_dtype(
                duckdb_dtype, self._version, deferred_time_zone
            )
            for column_name, duckdb_dtype in zip_strict(
                self.native.columns, self.native.types
            )
        }

    @property
    def columns(self) -> list[str]:
        if self._cached_columns is None:
            self._cached_columns = (
                list(self.schema)
                if self._cached_native_schema is not None
                else self.native.columns
            )
        return self._cached_columns

    def to_pandas(self) -> pd.DataFrame:
        # only if version is v1, keep around for backcompat
        return self.native.df()

    def to_arrow(self) -> pa.Table:
        # only if version is v1, keep around for backcompat
        return self.lazy().collect(Implementation.PYARROW).native  # type: ignore[no-any-return]

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, version=version)

    def _with_native(self, df: duckdb.DuckDBPyRelation) -> Self:
        return self.__class__(df, version=self._version)

    def group_by(
        self, keys: Sequence[str] | Sequence[DuckDBExpr], *, drop_null_keys: bool
    ) -> DuckDBGroupBy:
        from narwhals._duckdb.group_by import DuckDBGroupBy

        return DuckDBGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def rename(self, mapping: Mapping[str, str]) -> Self:
        df = self.native
        selection = (
            col(name).alias(mapping[name]) if name in mapping else col(name)
            for name in df.columns
        )
        return self._with_native(self.native.select(*selection))

    def join(
        self,
        other: Self,
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        native_how = "outer" if how == "full" else how

        if native_how == "cross":
            if self._backend_version < (1, 1, 4):
                msg = f"'duckdb>=1.1.4' is required for cross-join, found version: {self._backend_version}"
                raise NotImplementedError(msg)
            rel = self.native.set_alias("lhs").cross(other.native.set_alias("rhs"))
        else:
            # help mypy
            assert left_on is not None  # noqa: S101
            assert right_on is not None  # noqa: S101
            it = (
                col(f'lhs."{left}"') == col(f'rhs."{right}"')
                for left, right in zip_strict(left_on, right_on)
            )
            condition: Expression = reduce(and_, it)
            rel = self.native.set_alias("lhs").join(
                other.native.set_alias("rhs"),
                # NOTE: Fixed in `--pre` https://github.com/duckdb/duckdb/pull/16933
                condition=condition,  # type: ignore[arg-type, unused-ignore]
                how=native_how,
            )

        if native_how in {"inner", "left", "cross", "outer"}:
            select = [col(f'lhs."{x}"') for x in self.columns]
            for name in other.columns:
                col_in_lhs: bool = name in self.columns
                if native_how == "outer" and not col_in_lhs:
                    select.append(col(f'rhs."{name}"'))
                elif (native_how == "outer") or (
                    col_in_lhs and (right_on is None or name not in right_on)
                ):
                    select.append(col(f'rhs."{name}"').alias(f"{name}{suffix}"))
                elif right_on is None or name not in right_on:
                    select.append(col(name))
            res = rel.select(*select).set_alias(self.native.alias)
        else:  # semi, anti
            res = rel.select("lhs.*").set_alias(self.native.alias)

        return self._with_native(res)

    def join_asof(
        self,
        other: Self,
        *,
        left_on: str,
        right_on: str,
        by_left: Sequence[str] | None,
        by_right: Sequence[str] | None,
        strategy: AsofJoinStrategy,
        suffix: str,
    ) -> Self:
        lhs = self.native
        rhs = other.native
        conditions: list[Expression] = []
        if by_left is not None and by_right is not None:
            conditions.extend(
                col(f'lhs."{left}"') == col(f'rhs."{right}"')
                for left, right in zip_strict(by_left, by_right)
            )
        else:
            by_left = by_right = []
        if strategy == "backward":
            conditions.append(col(f'lhs."{left_on}"') >= col(f'rhs."{right_on}"'))
        elif strategy == "forward":
            conditions.append(col(f'lhs."{left_on}"') <= col(f'rhs."{right_on}"'))
        else:
            msg = "Only 'backward' and 'forward' strategies are currently supported for DuckDB"
            raise NotImplementedError(msg)
        condition: Expression = reduce(and_, conditions)
        select = ["lhs.*"]
        for name in rhs.columns:
            if name in lhs.columns and (
                right_on is None or name not in {right_on, *by_right}
            ):
                select.append(f'rhs."{name}" as "{name}{suffix}"')
            elif right_on is None or name not in {right_on, *by_right}:
                select.append(str(col(name)))
        # Replace with Python API call once
        # https://github.com/duckdb/duckdb/discussions/16947 is addressed.
        query = f"""
            SELECT {",".join(select)}
            FROM lhs
            ASOF LEFT JOIN rhs
            ON {condition}
            """  # noqa: S608
        return self._with_native(duckdb.sql(query))

    def collect_schema(self) -> dict[str, DType]:
        return self.schema

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
        flags = extend_bool(True, len(order_by)) if order_by and keep == "last" else None
        if keep == "none":
            expr = window_expression(
                F("count", StarExpression()),
                subset_,
                order_by or (),
                descending=flags,
                nulls_last=flags,
            )
        else:
            expr = window_expression(
                F("row_number"),
                subset_,
                order_by or (),
                descending=flags,
                nulls_last=flags,
            )
        return self._with_native(
            self.native.select(StarExpression(), expr.alias(tmp_name)).filter(
                col(tmp_name) == lit(1)
            )
        ).drop([tmp_name], strict=False)

    def sort(self, *by: str, descending: bool | Sequence[bool], nulls_last: bool) -> Self:
        descending = extend_bool(descending, len(by))
        if nulls_last:
            it = (
                col(name).nulls_last() if not desc else col(name).desc().nulls_last()
                for name, desc in zip_strict(by, descending)
            )
        else:
            it = (
                col(name).nulls_first() if not desc else col(name).desc().nulls_first()
                for name, desc in zip_strict(by, descending)
            )
        return self._with_native(self.native.sort(*it))

    def top_k(self, k: int, *, by: Iterable[str], reverse: bool | Sequence[bool]) -> Self:
        _rel = self.native
        by = list(by)
        if isinstance(reverse, bool):
            descending = extend_bool(not reverse, len(by))
        else:
            descending = tuple(not rev for rev in reverse)
        expr = window_expression(
            F("row_number"),
            order_by=by,
            descending=descending,
            nulls_last=extend_bool(True, len(by)),
        )
        condition = expr <= lit(k)
        query = f"""
            SELECT *
            FROM _rel
            QUALIFY {condition}
        """  # noqa: S608
        return self._with_native(duckdb.sql(query))

    def drop_nulls(self, subset: Sequence[str] | None) -> Self:
        subset_ = subset if subset is not None else self.columns
        keep_condition = reduce(and_, (col(name).isnotnull() for name in subset_))
        return self._with_native(self.native.filter(keep_condition))

    def explode(self, columns: Sequence[str]) -> Self:
        dtypes = self._version.dtypes
        schema = self.collect_schema()
        for name in columns:
            dtype = schema[name]
            if dtype != dtypes.List:
                msg = (
                    f"`explode` operation not supported for dtype `{dtype}`, "
                    "expected List type"
                )
                raise InvalidOperationError(msg)

        if len(columns) != 1:
            msg = (
                "Exploding on multiple columns is not supported with DuckDB backend since "
                "we cannot guarantee that the exploded columns have matching element counts."
            )
            raise NotImplementedError(msg)

        col_to_explode = col(columns[0])
        rel = self.native
        original_columns = self.columns

        not_null_condition = col_to_explode.isnotnull() & F("len", col_to_explode) > lit(
            0
        )
        non_null_rel = rel.filter(not_null_condition).select(
            *(
                F("unnest", col_to_explode).alias(name) if name in columns else name
                for name in original_columns
            )
        )

        null_rel = rel.filter(~not_null_condition).select(
            *(
                lit(None).alias(name) if name in columns else name
                for name in original_columns
            )
        )

        return self._with_native(non_null_rel.union(null_rel))

    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        index_ = [] if index is None else index
        on_ = [c for c in self.columns if c not in index_] if on is None else on

        if variable_name == "":
            msg = "`variable_name` cannot be empty string for duckdb backend."
            raise NotImplementedError(msg)

        if value_name == "":
            msg = "`value_name` cannot be empty string for duckdb backend."
            raise NotImplementedError(msg)

        unpivot_on = join_column_names(*on_)
        _rel = self.native
        # Replace with Python API once
        # https://github.com/duckdb/duckdb/discussions/16980 is addressed.
        query = f"""
            unpivot _rel
            on {unpivot_on}
            into
                name {col(variable_name)}
                value {col(value_name)}
            """
        return self._with_native(
            duckdb.sql(query).select(*[*index_, variable_name, value_name])
        )

    @requires.backend_version((1, 3))
    def with_row_index(self, name: str, order_by: Sequence[str]) -> Self:
        if order_by is None:
            msg = "Cannot pass `order_by` to `with_row_index` for DuckDB"
            raise TypeError(msg)
        expr = (window_expression(F("row_number"), order_by=order_by) - lit(1)).alias(
            name
        )
        return self._with_native(self.native.select(expr, StarExpression()))

    def sink_parquet(self, file: str | Path | BytesIO) -> None:
        _rel = self.native
        query = f"""
            COPY (SELECT * FROM _rel)
            TO '{file}'
            (FORMAT parquet)
            """  # noqa: S608
        duckdb.sql(query)

    gather_every = not_implemented.deprecated(
        "`LazyFrame.gather_every` is deprecated and will be removed in a future version."
    )
    tail = not_implemented.deprecated(
        "`LazyFrame.tail` is deprecated and will be removed in a future version."
    )
