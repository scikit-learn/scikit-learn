from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, Any

import duckdb
from duckdb import CoalesceOperator, Expression

from narwhals._duckdb.dataframe import DuckDBLazyFrame
from narwhals._duckdb.expr import DuckDBExpr
from narwhals._duckdb.selectors import DuckDBSelectorNamespace
from narwhals._duckdb.utils import (
    DeferredTimeZone,
    F,
    concat_str,
    duckdb_dtypes,
    function,
    lit,
    narwhals_to_native_dtype,
    when,
)
from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
)
from narwhals._sql.namespace import SQLNamespace
from narwhals._sql.when_then import SQLThen, SQLWhen
from narwhals._utils import Implementation

if TYPE_CHECKING:
    from collections.abc import Iterable

    from duckdb import DuckDBPyRelation  # noqa: F401

    from narwhals._utils import Version
    from narwhals.typing import ConcatMethod, IntoDType, NonNestedLiteral

VARCHAR = duckdb_dtypes.VARCHAR


class DuckDBNamespace(
    SQLNamespace[DuckDBLazyFrame, DuckDBExpr, "DuckDBPyRelation", Expression]
):
    _implementation: Implementation = Implementation.DUCKDB

    def __init__(self, *, version: Version) -> None:
        self._version = version

    @property
    def selectors(self) -> DuckDBSelectorNamespace:
        return DuckDBSelectorNamespace.from_namespace(self)

    @property
    def _expr(self) -> type[DuckDBExpr]:
        return DuckDBExpr

    @property
    def _lazyframe(self) -> type[DuckDBLazyFrame]:
        return DuckDBLazyFrame

    def _function(self, name: str, *args: Expression) -> Expression:  # type: ignore[override]
        return function(name, *args)

    def _lit(self, value: Any) -> Expression:
        return lit(value)

    def _when(
        self,
        condition: Expression,
        value: Expression,
        otherwise: Expression | None = None,
    ) -> Expression:
        if otherwise is None:
            return when(condition, value)
        return when(condition, value).otherwise(otherwise)

    def _coalesce(self, *exprs: Expression) -> Expression:
        return CoalesceOperator(*exprs)

    def concat(
        self, items: Iterable[DuckDBLazyFrame], *, how: ConcatMethod
    ) -> DuckDBLazyFrame:
        native_items = [item._native_frame for item in items]
        items = list(items)
        first = items[0]
        schema = first.schema
        if how == "vertical" and not all(x.schema == schema for x in items[1:]):
            msg = "inputs should all have the same schema"
            raise TypeError(msg)
        if how == "diagonal":
            res = first.native
            for _item in native_items[1:]:
                # TODO(unassigned): use relational API when available https://github.com/duckdb/duckdb/discussions/16996
                res = duckdb.sql("""
                    from res select * union all by name from _item select *
                """)
            return first._with_native(res)
        res = reduce(lambda x, y: x.union(y), native_items)
        return first._with_native(res)

    def concat_str(
        self, *exprs: DuckDBExpr, separator: str, ignore_nulls: bool
    ) -> DuckDBExpr:
        def func(df: DuckDBLazyFrame) -> list[Expression]:
            cols: Iterable[Expression] = chain.from_iterable(e(df) for e in exprs)
            if ignore_nulls:
                return [concat_str(*cols, separator=separator)]
            cols = tuple(cols)
            null_mask = reduce(operator.or_, (s.isnull() for s in cols))
            cols_str = (c.cast(VARCHAR) for c in cols)
            return [when(~null_mask, concat_str(*cols_str, separator=separator))]

        return self._expr(
            call=func,
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def mean_horizontal(self, *exprs: DuckDBExpr) -> DuckDBExpr:
        def func(cols: Iterable[Expression]) -> Expression:
            cols = tuple(cols)
            total = reduce(operator.add, (CoalesceOperator(col, lit(0)) for col in cols))
            count = reduce(
                operator.add, (col.isnotnull().cast(duckdb_dtypes.BIGINT) for col in cols)
            )
            return total / count

        return self._expr._from_elementwise_horizontal_op(func, *exprs)

    def when(self, predicate: DuckDBExpr) -> DuckDBWhen:
        return DuckDBWhen.from_expr(predicate, context=self)

    def lit(self, value: NonNestedLiteral, dtype: IntoDType | None) -> DuckDBExpr:
        def func(df: DuckDBLazyFrame) -> list[Expression]:
            tz = DeferredTimeZone(df.native)
            if dtype is not None:
                target = narwhals_to_native_dtype(dtype, self._version, tz)
                return [lit(value).cast(target)]
            return [lit(value)]

        return self._expr(
            func,
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
        )

    def len(self) -> DuckDBExpr:
        def func(_df: DuckDBLazyFrame) -> list[Expression]:
            return [F("count")]

        return self._expr(
            call=func,
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            version=self._version,
        )


class DuckDBWhen(SQLWhen["DuckDBLazyFrame", Expression, DuckDBExpr]):
    @property
    def _then(self) -> type[DuckDBThen]:
        return DuckDBThen


class DuckDBThen(SQLThen["DuckDBLazyFrame", Expression, DuckDBExpr], DuckDBExpr): ...
