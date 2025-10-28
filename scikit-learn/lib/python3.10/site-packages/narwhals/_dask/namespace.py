from __future__ import annotations

import operator
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, cast

import dask.dataframe as dd
import pandas as pd

from narwhals._compliant import (
    CompliantThen,
    CompliantWhen,
    DepthTrackingNamespace,
    LazyNamespace,
)
from narwhals._dask.dataframe import DaskLazyFrame
from narwhals._dask.expr import DaskExpr
from narwhals._dask.selectors import DaskSelectorNamespace
from narwhals._dask.utils import (
    align_series_full_broadcast,
    narwhals_to_native_dtype,
    validate_comparand,
)
from narwhals._expression_parsing import (
    ExprKind,
    combine_alias_output_names,
    combine_evaluate_output_names,
)
from narwhals._utils import Implementation, zip_strict

if TYPE_CHECKING:
    from collections.abc import Iterable, Iterator, Sequence

    import dask.dataframe.dask_expr as dx

    from narwhals._compliant.typing import ScalarKwargs
    from narwhals._utils import Version
    from narwhals.typing import ConcatMethod, IntoDType, NonNestedLiteral


class DaskNamespace(
    LazyNamespace[DaskLazyFrame, DaskExpr, dd.DataFrame],
    DepthTrackingNamespace[DaskLazyFrame, DaskExpr],
):
    _implementation: Implementation = Implementation.DASK

    @property
    def selectors(self) -> DaskSelectorNamespace:
        return DaskSelectorNamespace.from_namespace(self)

    @property
    def _expr(self) -> type[DaskExpr]:
        return DaskExpr

    @property
    def _lazyframe(self) -> type[DaskLazyFrame]:
        return DaskLazyFrame

    def __init__(self, *, version: Version) -> None:
        self._version = version

    def lit(self, value: NonNestedLiteral, dtype: IntoDType | None) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            if dtype is not None:
                native_dtype = narwhals_to_native_dtype(dtype, self._version)
                native_pd_series = pd.Series([value], dtype=native_dtype, name="literal")
            else:
                native_pd_series = pd.Series([value], name="literal")
            npartitions = df._native_frame.npartitions
            dask_series = dd.from_pandas(native_pd_series, npartitions=npartitions)
            return [dask_series[0].to_series()]

        return self._expr(
            func,
            depth=0,
            function_name="lit",
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            version=self._version,
        )

    def len(self) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            # We don't allow dataframes with 0 columns, so `[0]` is safe.
            return [df._native_frame[df.columns[0]].size.to_series()]

        return self._expr(
            func,
            depth=0,
            function_name="len",
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            version=self._version,
        )

    def all_horizontal(self, *exprs: DaskExpr, ignore_nulls: bool) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            series: Iterator[dx.Series] = chain.from_iterable(e(df) for e in exprs)
            # Note on `ignore_nulls`: Dask doesn't support storing arbitrary Python
            # objects in `object` dtype, so we don't need the same check we have for pandas-like.
            if ignore_nulls:
                # NumPy-backed 'bool' dtype can't contain nulls so doesn't need filling.
                series = (s if s.dtype == "bool" else s.fillna(True) for s in series)
            return [reduce(operator.and_, align_series_full_broadcast(df, *series))]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="all_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def any_horizontal(self, *exprs: DaskExpr, ignore_nulls: bool) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            series: Iterator[dx.Series] = chain.from_iterable(e(df) for e in exprs)
            if ignore_nulls:
                series = (s if s.dtype == "bool" else s.fillna(False) for s in series)
            return [reduce(operator.or_, align_series_full_broadcast(df, *series))]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="any_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def sum_horizontal(self, *exprs: DaskExpr) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            series = align_series_full_broadcast(
                df, *(s for _expr in exprs for s in _expr(df))
            )
            return [dd.concat(series, axis=1).sum(axis=1)]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="sum_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def concat(
        self, items: Iterable[DaskLazyFrame], *, how: ConcatMethod
    ) -> DaskLazyFrame:
        if not items:
            msg = "No items to concatenate"  # pragma: no cover
            raise AssertionError(msg)
        dfs = [i._native_frame for i in items]
        cols_0 = dfs[0].columns
        if how == "vertical":
            for i, df in enumerate(dfs[1:], start=1):
                cols_current = df.columns
                if not (
                    (len(cols_current) == len(cols_0)) and (cols_current == cols_0).all()
                ):
                    msg = (
                        "unable to vstack, column names don't match:\n"
                        f"   - dataframe 0: {cols_0.to_list()}\n"
                        f"   - dataframe {i}: {cols_current.to_list()}\n"
                    )
                    raise TypeError(msg)
            return DaskLazyFrame(
                dd.concat(dfs, axis=0, join="inner"), version=self._version
            )
        if how == "diagonal":
            return DaskLazyFrame(
                dd.concat(dfs, axis=0, join="outer"), version=self._version
            )

        raise NotImplementedError

    def mean_horizontal(self, *exprs: DaskExpr) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            expr_results = [s for _expr in exprs for s in _expr(df)]
            series = align_series_full_broadcast(df, *(s.fillna(0) for s in expr_results))
            non_na = align_series_full_broadcast(
                df, *(1 - s.isna() for s in expr_results)
            )
            num = reduce(lambda x, y: x + y, series)  # pyright: ignore[reportOperatorIssue]
            den = reduce(lambda x, y: x + y, non_na)  # pyright: ignore[reportOperatorIssue]
            return [cast("dx.Series", num / den)]  # pyright: ignore[reportOperatorIssue]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="mean_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def min_horizontal(self, *exprs: DaskExpr) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            series = align_series_full_broadcast(
                df, *(s for _expr in exprs for s in _expr(df))
            )

            return [dd.concat(series, axis=1).min(axis=1)]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="min_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def max_horizontal(self, *exprs: DaskExpr) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            series = align_series_full_broadcast(
                df, *(s for _expr in exprs for s in _expr(df))
            )

            return [dd.concat(series, axis=1).max(axis=1)]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="max_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )

    def when(self, predicate: DaskExpr) -> DaskWhen:
        return DaskWhen.from_expr(predicate, context=self)

    def concat_str(
        self, *exprs: DaskExpr, separator: str, ignore_nulls: bool
    ) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            expr_results = [s for _expr in exprs for s in _expr(df)]
            series = (
                s.astype(str) for s in align_series_full_broadcast(df, *expr_results)
            )
            null_mask = [s.isna() for s in align_series_full_broadcast(df, *expr_results)]

            if not ignore_nulls:
                null_mask_result = reduce(operator.or_, null_mask)
                result = reduce(lambda x, y: x + separator + y, series).where(
                    ~null_mask_result, None
                )
            else:
                init_value, *values = [
                    s.where(~nm, "") for s, nm in zip_strict(series, null_mask)
                ]

                separators = (
                    nm.map({True: "", False: separator}, meta=str)
                    for nm in null_mask[:-1]
                )
                result = reduce(
                    operator.add,
                    (s + v for s, v in zip_strict(separators, values)),
                    init_value,
                )

            return [result]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="concat_str",
            evaluate_output_names=getattr(
                exprs[0], "_evaluate_output_names", lambda _df: ["literal"]
            ),
            alias_output_names=getattr(exprs[0], "_alias_output_names", None),
            version=self._version,
        )

    def coalesce(self, *exprs: DaskExpr) -> DaskExpr:
        def func(df: DaskLazyFrame) -> list[dx.Series]:
            series = align_series_full_broadcast(
                df, *(s for _expr in exprs for s in _expr(df))
            )
            return [reduce(lambda x, y: x.fillna(y), series)]

        return self._expr(
            call=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="coalesce",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            version=self._version,
        )


class DaskWhen(CompliantWhen[DaskLazyFrame, "dx.Series", DaskExpr]):  # pyright: ignore[reportInvalidTypeArguments]
    @property
    def _then(self) -> type[DaskThen]:
        return DaskThen

    def __call__(self, df: DaskLazyFrame) -> Sequence[dx.Series]:
        then_value = (
            self._then_value(df)[0]
            if isinstance(self._then_value, DaskExpr)
            else self._then_value
        )
        otherwise_value = (
            self._otherwise_value(df)[0]
            if isinstance(self._otherwise_value, DaskExpr)
            else self._otherwise_value
        )

        condition = self._condition(df)[0]
        # re-evaluate DataFrame if the condition aggregates to force
        #   then/otherwise to be evaluated against the aggregated frame
        assert self._condition._metadata is not None  # noqa: S101
        if self._condition._metadata.is_scalar_like:
            new_df = df._with_native(condition.to_frame())
            condition = self._condition.broadcast(ExprKind.AGGREGATION)(df)[0]
            df = new_df

        if self._otherwise_value is None:
            (condition, then_series) = align_series_full_broadcast(
                df, condition, then_value
            )
            validate_comparand(condition, then_series)
            return [then_series.where(condition)]  # pyright: ignore[reportArgumentType]
        (condition, then_series, otherwise_series) = align_series_full_broadcast(
            df, condition, then_value, otherwise_value
        )
        validate_comparand(condition, then_series)
        validate_comparand(condition, otherwise_series)
        return [then_series.where(condition, otherwise_series)]  # pyright: ignore[reportArgumentType]


class DaskThen(CompliantThen[DaskLazyFrame, "dx.Series", DaskExpr, DaskWhen], DaskExpr):  # pyright: ignore[reportInvalidTypeArguments]
    _depth: int = 0
    _scalar_kwargs: ScalarKwargs = {}  # noqa: RUF012
    _function_name: str = "whenthen"
