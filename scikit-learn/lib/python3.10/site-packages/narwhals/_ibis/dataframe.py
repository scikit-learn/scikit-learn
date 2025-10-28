from __future__ import annotations

import operator
from io import BytesIO
from typing import TYPE_CHECKING, Any, cast

import ibis
import ibis.expr.types as ir

from narwhals._ibis.expr import IbisExpr
from narwhals._ibis.utils import evaluate_exprs, lit, native_to_narwhals_dtype
from narwhals._sql.dataframe import SQLLazyFrame
from narwhals._utils import (
    Implementation,
    ValidateBackendVersion,
    Version,
    generate_temporary_column_name,
    not_implemented,
    parse_columns_to_drop,
    to_pyarrow_table,
    zip_strict,
)
from narwhals.exceptions import InvalidOperationError

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Mapping, Sequence
    from pathlib import Path
    from types import ModuleType

    import pandas as pd
    import pyarrow as pa
    from ibis.expr.operations import Binary
    from typing_extensions import Self, TypeAlias, TypeIs

    from narwhals._compliant.typing import CompliantDataFrameAny
    from narwhals._ibis.group_by import IbisGroupBy
    from narwhals._ibis.namespace import IbisNamespace
    from narwhals._ibis.series import IbisInterchangeSeries
    from narwhals._typing import _EagerAllowedImpl
    from narwhals._utils import _LimitedContext
    from narwhals.dataframe import LazyFrame
    from narwhals.dtypes import DType
    from narwhals.stable.v1 import DataFrame as DataFrameV1
    from narwhals.typing import AsofJoinStrategy, JoinStrategy, UniqueKeepStrategy

    JoinPredicates: TypeAlias = "Sequence[ir.BooleanColumn] | Sequence[str]"


class IbisLazyFrame(
    SQLLazyFrame["IbisExpr", "ir.Table", "LazyFrame[ir.Table] | DataFrameV1[ir.Table]"],
    ValidateBackendVersion,
):
    _implementation = Implementation.IBIS

    def __init__(
        self, df: ir.Table, *, version: Version, validate_backend_version: bool = False
    ) -> None:
        self._native_frame: ir.Table = df
        self._version = version
        self._cached_schema: dict[str, DType] | None = None
        self._cached_columns: list[str] | None = None
        if validate_backend_version:
            self._validate_backend_version()

    @staticmethod
    def _is_native(obj: ir.Table | Any) -> TypeIs[ir.Table]:
        return isinstance(obj, ir.Table)

    @classmethod
    def from_native(cls, data: ir.Table, /, *, context: _LimitedContext) -> Self:
        return cls(data, version=context._version)

    def to_narwhals(self) -> LazyFrame[ir.Table] | DataFrameV1[ir.Table]:
        if self._version is Version.V1:
            from narwhals.stable.v1 import DataFrame

            return DataFrame(self, level="interchange")
        return self._version.lazyframe(self, level="lazy")

    def __narwhals_dataframe__(self) -> Self:  # pragma: no cover
        # Keep around for backcompat.
        if self._version is not Version.V1:
            msg = "__narwhals_dataframe__ is not implemented for IbisLazyFrame"
            raise AttributeError(msg)
        return self

    def __narwhals_lazyframe__(self) -> Self:
        return self

    def __native_namespace__(self) -> ModuleType:
        return ibis

    def __narwhals_namespace__(self) -> IbisNamespace:
        from narwhals._ibis.namespace import IbisNamespace

        return IbisNamespace(version=self._version)

    def get_column(self, name: str) -> IbisInterchangeSeries:
        from narwhals._ibis.series import IbisInterchangeSeries

        return IbisInterchangeSeries(self.native.select(name), version=self._version)

    def _iter_columns(self) -> Iterator[ir.Expr]:
        for name in self.columns:
            yield self.native[name]

    def collect(
        self, backend: _EagerAllowedImpl | None, **kwargs: Any
    ) -> CompliantDataFrameAny:
        if backend is None or backend is Implementation.PYARROW:
            from narwhals._arrow.dataframe import ArrowDataFrame

            return ArrowDataFrame(
                to_pyarrow_table(self.native.to_pyarrow()),
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is Implementation.PANDAS:
            from narwhals._pandas_like.dataframe import PandasLikeDataFrame

            return PandasLikeDataFrame(
                self.native.to_pandas(),
                implementation=Implementation.PANDAS,
                validate_backend_version=True,
                version=self._version,
                validate_column_names=True,
            )

        if backend is Implementation.POLARS:
            from narwhals._polars.dataframe import PolarsDataFrame

            return PolarsDataFrame(
                self.native.to_polars(),
                validate_backend_version=True,
                version=self._version,
            )

        msg = f"Unsupported `backend` value: {backend}"  # pragma: no cover
        raise ValueError(msg)  # pragma: no cover

    def head(self, n: int) -> Self:
        return self._with_native(self.native.head(n))

    def simple_select(self, *column_names: str) -> Self:
        return self._with_native(self.native.select(*column_names))

    def aggregate(self, *exprs: IbisExpr) -> Self:
        selection = [
            cast("ir.Scalar", val.name(name))
            for name, val in evaluate_exprs(self, *exprs)
        ]
        return self._with_native(self.native.aggregate(selection))

    def select(self, *exprs: IbisExpr) -> Self:
        selection = [val.name(name) for name, val in evaluate_exprs(self, *exprs)]
        if not selection:
            msg = "At least one expression must be provided to `select` with the Ibis backend."
            raise ValueError(msg)

        t = self.native.select(*selection)
        return self._with_native(t)

    def drop(self, columns: Sequence[str], *, strict: bool) -> Self:
        columns_to_drop = parse_columns_to_drop(self, columns, strict=strict)
        selection = (col for col in self.columns if col not in columns_to_drop)
        return self._with_native(self.native.select(*selection))

    def lazy(self, backend: None = None, **_: None) -> Self:
        # The `backend`` argument has no effect but we keep it here for
        # backwards compatibility because in `narwhals.stable.v1`
        # function `.from_native()` will return a DataFrame for Ibis.

        if backend is not None:  # pragma: no cover
            msg = "`backend` argument is not supported for Ibis"
            raise ValueError(msg)
        return self

    def with_columns(self, *exprs: IbisExpr) -> Self:
        new_columns_map = dict(evaluate_exprs(self, *exprs))
        return self._with_native(self.native.mutate(**new_columns_map))

    def filter(self, predicate: IbisExpr) -> Self:
        # `[0]` is safe as the predicate's expression only returns a single column
        mask = cast("ir.BooleanValue", predicate(self)[0])
        return self._with_native(self.native.filter(mask))

    @property
    def schema(self) -> dict[str, DType]:
        if self._cached_schema is None:
            # Note: prefer `self._cached_schema` over `functools.cached_property`
            # due to Python3.13 failures.
            self._cached_schema = {
                name: native_to_narwhals_dtype(dtype, self._version)
                for name, dtype in self.native.schema().fields.items()
            }
        return self._cached_schema

    @property
    def columns(self) -> list[str]:
        if self._cached_columns is None:
            self._cached_columns = (
                list(self.schema)
                if self._cached_schema is not None
                else list(self.native.columns)
            )
        return self._cached_columns

    def to_pandas(self) -> pd.DataFrame:
        # only if version is v1, keep around for backcompat
        return self.native.to_pandas()

    def to_arrow(self) -> pa.Table:
        # only if version is v1, keep around for backcompat
        return self.native.to_pyarrow()

    def _with_version(self, version: Version) -> Self:
        return self.__class__(self.native, version=version)

    def _with_native(self, df: ir.Table) -> Self:
        return self.__class__(df, version=self._version)

    def group_by(
        self, keys: Sequence[str] | Sequence[IbisExpr], *, drop_null_keys: bool
    ) -> IbisGroupBy:
        from narwhals._ibis.group_by import IbisGroupBy

        return IbisGroupBy(self, keys, drop_null_keys=drop_null_keys)

    def rename(self, mapping: Mapping[str, str]) -> Self:
        def _rename(col: str) -> str:
            return mapping.get(col, col)

        return self._with_native(self.native.rename(_rename))

    @staticmethod
    def _join_drop_duplicate_columns(df: ir.Table, columns: Iterable[str], /) -> ir.Table:
        """Ibis adds a suffix to the right table col, even when it matches the left during a join."""
        duplicates = set(df.columns).intersection(columns)
        return df.drop(*duplicates) if duplicates else df

    def join(
        self,
        other: Self,
        *,
        how: JoinStrategy,
        left_on: Sequence[str] | None,
        right_on: Sequence[str] | None,
        suffix: str,
    ) -> Self:
        how_native = "outer" if how == "full" else how
        rname = "{name}" + suffix
        if other == self:
            # Ibis does not support self-references unless created as a view
            other = self._with_native(other.native.view())
        if how_native == "cross":
            joined = self.native.join(other.native, how=how_native, rname=rname)
            return self._with_native(joined)
        # help mypy
        assert left_on is not None  # noqa: S101
        assert right_on is not None  # noqa: S101
        predicates = self._convert_predicates(other, left_on, right_on)
        joined = self.native.join(other.native, predicates, how=how_native, rname=rname)
        if how_native == "left":
            right_names = (n + suffix for n in right_on)
            joined = self._join_drop_duplicate_columns(joined, right_names)
            it = (cast("Binary", p.op()) for p in predicates if not isinstance(p, str))
            to_drop = []
            for pred in it:
                right = pred.right.name
                # Mirrors how polars works.
                if right not in self.columns and pred.left.name != right:
                    to_drop.append(right)
            if to_drop:
                joined = joined.drop(*to_drop)
        return self._with_native(joined)

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
        rname = "{name}" + suffix
        strategy_op = {"backward": operator.ge, "forward": operator.le}
        predicates: JoinPredicates = []
        if op := strategy_op.get(strategy):
            on: ir.BooleanColumn = op(self.native[left_on], other.native[right_on])
        else:
            msg = "Only `backward` and `forward` strategies are currently supported for Ibis"
            raise NotImplementedError(msg)
        if by_left is not None and by_right is not None:
            predicates = self._convert_predicates(other, by_left, by_right)
        joined = self.native.asof_join(other.native, on, predicates, rname=rname)
        joined = self._join_drop_duplicate_columns(joined, [right_on + suffix])
        if by_right is not None:
            right_names = (n + suffix for n in by_right)
            joined = self._join_drop_duplicate_columns(joined, right_names)
        return self._with_native(joined)

    def _convert_predicates(
        self, other: Self, left_on: Sequence[str], right_on: Sequence[str]
    ) -> JoinPredicates:
        if left_on == right_on:
            return left_on
        return [
            cast("ir.BooleanColumn", (self.native[left] == other.native[right]))
            for left, right in zip_strict(left_on, right_on)
        ]

    def collect_schema(self) -> dict[str, DType]:
        return {
            name: native_to_narwhals_dtype(dtype, self._version)
            for name, dtype in self.native.schema().fields.items()
        }

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
        if order_by and keep == "last":
            order_by_ = IbisExpr._sort(*order_by, descending=True, nulls_last=True)
        elif order_by:
            order_by_ = IbisExpr._sort(*order_by, descending=False, nulls_last=False)
        else:
            order_by_ = lit(1)
        window = ibis.window(group_by=subset_, order_by=order_by_)
        if keep == "none":
            expr = self.native.count().over(window)
        else:
            expr = ibis.row_number().over(window) + lit(1)
        df = (
            self.native.mutate(**{tmp_name: expr})
            .filter(ibis._[tmp_name] == lit(1))
            .drop(tmp_name)
        )
        return self._with_native(df)

    def sort(self, *by: str, descending: bool | Sequence[bool], nulls_last: bool) -> Self:
        from narwhals._ibis.expr import IbisExpr

        cols = IbisExpr._sort(*by, descending=descending, nulls_last=nulls_last)
        return self._with_native(self.native.order_by(*cols))

    def top_k(self, k: int, *, by: Iterable[str], reverse: bool | Sequence[bool]) -> Self:
        from narwhals._ibis.expr import IbisExpr

        desc = not reverse if isinstance(reverse, bool) else [not el for el in reverse]
        cols = IbisExpr._sort(*by, descending=desc, nulls_last=True)
        return self._with_native(self.native.order_by(*cols).head(k))

    def drop_nulls(self, subset: Sequence[str] | None) -> Self:
        subset_ = subset if subset is not None else self.columns
        return self._with_native(self.native.drop_null(subset_))

    def explode(self, columns: Sequence[str]) -> Self:
        dtypes = self._version.dtypes
        schema = self.collect_schema()
        for col in columns:
            dtype = schema[col]

            if dtype != dtypes.List:
                msg = (
                    f"`explode` operation not supported for dtype `{dtype}`, "
                    "expected List type"
                )
                raise InvalidOperationError(msg)

        if len(columns) != 1:
            msg = (
                "Exploding on multiple columns is not supported with Ibis backend since "
                "we cannot guarantee that the exploded columns have matching element counts."
            )
            raise NotImplementedError(msg)

        return self._with_native(self.native.unnest(columns[0], keep_empty=True))

    def unpivot(
        self,
        on: Sequence[str] | None,
        index: Sequence[str] | None,
        variable_name: str,
        value_name: str,
    ) -> Self:
        import ibis.selectors as s

        index_: Sequence[str] = [] if index is None else index
        on_: Sequence[str] = (
            [c for c in self.columns if c not in index_] if on is None else on
        )

        # Discard columns not in the index
        final_columns = list(dict.fromkeys([*index_, variable_name, value_name]))

        unpivoted = self.native.pivot_longer(
            s.cols(*on_), names_to=variable_name, values_to=value_name
        )
        return self._with_native(unpivoted.select(*final_columns))

    def with_row_index(self, name: str, order_by: Sequence[str]) -> Self:
        to_select = [
            ibis.row_number().over(ibis.window(order_by=order_by)).name(name),
            ibis.selectors.all(),
        ]
        return self._with_native(self.native.select(*to_select))

    def sink_parquet(self, file: str | Path | BytesIO) -> None:
        if isinstance(file, BytesIO):  # pragma: no cover
            msg = "Writing to BytesIO is not supported for Ibis backend."
            raise NotImplementedError(msg)
        self.native.to_parquet(file)

    gather_every = not_implemented.deprecated(
        "`LazyFrame.gather_every` is deprecated and will be removed in a future version."
    )
    tail = not_implemented.deprecated(
        "`LazyFrame.tail` is deprecated and will be removed in a future version."
    )

    # Intentionally not implemented, as Ibis does its own expression rewriting.
    _evaluate_window_expr = not_implemented()
