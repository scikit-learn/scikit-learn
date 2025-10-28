from __future__ import annotations

import operator
import warnings
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, Any, Literal, Protocol, overload

from narwhals._compliant import CompliantThen, EagerNamespace, EagerWhen
from narwhals._expression_parsing import (
    combine_alias_output_names,
    combine_evaluate_output_names,
)
from narwhals._pandas_like.dataframe import PandasLikeDataFrame
from narwhals._pandas_like.expr import PandasLikeExpr
from narwhals._pandas_like.selectors import PandasSelectorNamespace
from narwhals._pandas_like.series import PandasLikeSeries
from narwhals._pandas_like.typing import NativeDataFrameT, NativeSeriesT
from narwhals._pandas_like.utils import is_non_nullable_boolean
from narwhals._utils import zip_strict

if TYPE_CHECKING:
    from collections.abc import Iterable, Sequence

    from typing_extensions import TypeAlias

    from narwhals._compliant.typing import ScalarKwargs
    from narwhals._utils import Implementation, Version
    from narwhals.typing import IntoDType, NonNestedLiteral


Incomplete: TypeAlias = Any
"""Escape hatch, but leaving a trace that this isn't ideal."""


_Vertical: TypeAlias = Literal[0]
_Horizontal: TypeAlias = Literal[1]
Axis: TypeAlias = Literal[_Vertical, _Horizontal]

VERTICAL: _Vertical = 0
HORIZONTAL: _Horizontal = 1


class PandasLikeNamespace(
    EagerNamespace[
        PandasLikeDataFrame,
        PandasLikeSeries,
        PandasLikeExpr,
        NativeDataFrameT,
        NativeSeriesT,
    ]
):
    @property
    def _dataframe(self) -> type[PandasLikeDataFrame]:
        return PandasLikeDataFrame

    @property
    def _expr(self) -> type[PandasLikeExpr]:
        return PandasLikeExpr

    @property
    def _series(self) -> type[PandasLikeSeries]:
        return PandasLikeSeries

    @property
    def selectors(self) -> PandasSelectorNamespace:
        return PandasSelectorNamespace.from_namespace(self)

    def __init__(self, implementation: Implementation, version: Version) -> None:
        self._implementation = implementation
        self._version = version

    def coalesce(self, *exprs: PandasLikeExpr) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            series = (s for _expr in exprs for s in _expr(df))
            return [
                reduce(lambda x, y: x.fill_null(y, strategy=None, limit=None), series)
            ]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="coalesce",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def lit(self, value: NonNestedLiteral, dtype: IntoDType | None) -> PandasLikeExpr:
        def _lit_pandas_series(df: PandasLikeDataFrame) -> PandasLikeSeries:
            pandas_series = self._series.from_iterable(
                data=[value],
                name="literal",
                index=df._native_frame.index[0:1],
                context=self,
            )
            if dtype:
                return pandas_series.cast(dtype)
            return pandas_series

        return PandasLikeExpr(
            lambda df: [_lit_pandas_series(df)],
            depth=0,
            function_name="lit",
            evaluate_output_names=lambda _df: ["literal"],
            alias_output_names=None,
            implementation=self._implementation,
            version=self._version,
        )

    def len(self) -> PandasLikeExpr:
        return PandasLikeExpr(
            lambda df: [
                self._series.from_iterable(
                    [len(df._native_frame)], name="len", index=[0], context=self
                )
            ],
            depth=0,
            function_name="len",
            evaluate_output_names=lambda _df: ["len"],
            alias_output_names=None,
            implementation=self._implementation,
            version=self._version,
        )

    # --- horizontal ---
    def sum_horizontal(self, *exprs: PandasLikeExpr) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            it = chain.from_iterable(expr(df) for expr in exprs)
            native_series = (s.fill_null(0, None, None) for s in it)
            return [reduce(operator.add, native_series)]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="sum_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def all_horizontal(
        self, *exprs: PandasLikeExpr, ignore_nulls: bool
    ) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            series = [s for _expr in exprs for s in _expr(df)]
            if not ignore_nulls and any(
                s.native.dtype == "object" and s.is_null().any() for s in series
            ):
                # classical NumPy boolean columns don't support missing values, so
                # only do the full scan with `is_null` if we have `object` dtype.
                msg = "Cannot use `ignore_nulls=False` in `all_horizontal` for non-nullable NumPy-backed pandas Series when nulls are present."
                raise ValueError(msg)
            it = (
                (
                    # NumPy-backed 'bool' dtype can't contain nulls so doesn't need filling.
                    s if is_non_nullable_boolean(s) else s.fill_null(True, None, None)
                    for s in series
                )
                if ignore_nulls
                else iter(series)
            )
            return [reduce(operator.and_, it)]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="all_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def any_horizontal(
        self, *exprs: PandasLikeExpr, ignore_nulls: bool
    ) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            series = [s for _expr in exprs for s in _expr(df)]
            if not ignore_nulls and any(
                s.native.dtype == "object" and s.is_null().any() for s in series
            ):
                # classical NumPy boolean columns don't support missing values, so
                # only do the full scan with `is_null` if we have `object` dtype.
                msg = "Cannot use `ignore_nulls=False` in `any_horizontal` for non-nullable NumPy-backed pandas Series when nulls are present."
                raise ValueError(msg)
            it = (
                (
                    # NumPy-backed 'bool' dtype can't contain nulls so doesn't need filling.
                    s if is_non_nullable_boolean(s) else s.fill_null(False, None, None)
                    for s in series
                )
                if ignore_nulls
                else iter(series)
            )
            return [reduce(operator.or_, it)]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="any_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def mean_horizontal(self, *exprs: PandasLikeExpr) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            expr_results = [s for _expr in exprs for s in _expr(df)]
            series = (s.fill_null(0, strategy=None, limit=None) for s in expr_results)
            non_na = (1 - s.is_null() for s in expr_results)
            return [reduce(operator.add, series) / reduce(operator.add, non_na)]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="mean_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def min_horizontal(self, *exprs: PandasLikeExpr) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            series = list(chain.from_iterable(expr(df) for expr in exprs))
            return [
                PandasLikeSeries(
                    self.concat(
                        (s.to_frame() for s in series), how="horizontal"
                    )._native_frame.min(axis=1),
                    implementation=self._implementation,
                    version=self._version,
                ).alias(series[0].name)
            ]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="min_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    def max_horizontal(self, *exprs: PandasLikeExpr) -> PandasLikeExpr:
        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            series = list(chain.from_iterable(expr(df) for expr in exprs))
            return [
                PandasLikeSeries(
                    self.concat(
                        (s.to_frame() for s in series), how="horizontal"
                    ).native.max(axis=1),
                    implementation=self._implementation,
                    version=self._version,
                ).alias(series[0].name)
            ]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="max_horizontal",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )

    @property
    def _concat(self) -> _NativeConcat[NativeDataFrameT, NativeSeriesT]:
        """Concatenate pandas objects along a particular axis.

        Return the **native** equivalent of `pd.concat`.
        """
        return self._implementation.to_native_namespace().concat

    def _concat_diagonal(self, dfs: Sequence[NativeDataFrameT], /) -> NativeDataFrameT:
        if self._implementation.is_pandas() and self._backend_version < (3,):
            return self._concat(dfs, axis=VERTICAL, copy=False)
        return self._concat(dfs, axis=VERTICAL)

    def _concat_horizontal(
        self, dfs: Sequence[NativeDataFrameT | NativeSeriesT], /
    ) -> NativeDataFrameT:
        if self._implementation.is_cudf():
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    "ignore",
                    message="The behavior of array concatenation with empty entries is deprecated",
                    category=FutureWarning,
                )
                return self._concat(dfs, axis=HORIZONTAL)
        elif self._implementation.is_pandas() and self._backend_version < (3,):
            return self._concat(dfs, axis=HORIZONTAL, copy=False)
        return self._concat(dfs, axis=HORIZONTAL)

    def _concat_vertical(self, dfs: Sequence[NativeDataFrameT], /) -> NativeDataFrameT:
        cols_0 = dfs[0].columns
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
        if self._implementation.is_pandas() and self._backend_version < (3,):
            return self._concat(dfs, axis=VERTICAL, copy=False)
        return self._concat(dfs, axis=VERTICAL)

    def when(self, predicate: PandasLikeExpr) -> PandasWhen[NativeSeriesT]:
        return PandasWhen[NativeSeriesT].from_expr(predicate, context=self)

    def concat_str(
        self, *exprs: PandasLikeExpr, separator: str, ignore_nulls: bool
    ) -> PandasLikeExpr:
        string = self._version.dtypes.String()

        def func(df: PandasLikeDataFrame) -> list[PandasLikeSeries]:
            expr_results = [s for _expr in exprs for s in _expr(df)]
            series = [s.cast(string) for s in expr_results]
            null_mask = [s.is_null() for s in expr_results]

            if not ignore_nulls:
                null_mask_result = reduce(operator.or_, null_mask)
                result = reduce(lambda x, y: x + separator + y, series).zip_with(
                    ~null_mask_result, None
                )
            else:
                # NOTE: Trying to help `mypy` later
                # error: Cannot determine type of "values"  [has-type]
                values: list[PandasLikeSeries]
                init_value, *values = (
                    s.zip_with(~nm, "") for s, nm in zip_strict(series, null_mask)
                )
                sep_array = init_value._with_native(
                    init_value.__native_namespace__().Series(
                        separator,
                        name="sep",
                        index=init_value.native.index,
                        dtype=init_value.native.dtype,
                    )
                )
                separators = (sep_array.zip_with(~nm, "") for nm in null_mask[:-1])
                result = reduce(
                    operator.add,
                    (s + v for s, v in zip_strict(separators, values)),
                    init_value,
                )

            return [result]

        return self._expr._from_callable(
            func=func,
            depth=max(x._depth for x in exprs) + 1,
            function_name="concat_str",
            evaluate_output_names=combine_evaluate_output_names(*exprs),
            alias_output_names=combine_alias_output_names(*exprs),
            context=self,
        )


class _NativeConcat(Protocol[NativeDataFrameT, NativeSeriesT]):
    @overload
    def __call__(
        self,
        objs: Iterable[NativeDataFrameT],
        *,
        axis: _Vertical,
        copy: bool | None = ...,
    ) -> NativeDataFrameT: ...
    @overload
    def __call__(
        self, objs: Iterable[NativeSeriesT], *, axis: _Vertical, copy: bool | None = ...
    ) -> NativeSeriesT: ...
    @overload
    def __call__(
        self,
        objs: Iterable[NativeDataFrameT | NativeSeriesT],
        *,
        axis: _Horizontal,
        copy: bool | None = ...,
    ) -> NativeDataFrameT: ...
    @overload
    def __call__(
        self,
        objs: Iterable[NativeDataFrameT | NativeSeriesT],
        *,
        axis: Axis,
        copy: bool | None = ...,
    ) -> NativeDataFrameT | NativeSeriesT: ...

    def __call__(
        self,
        objs: Iterable[NativeDataFrameT | NativeSeriesT],
        *,
        axis: Axis,
        copy: bool | None = None,
    ) -> NativeDataFrameT | NativeSeriesT: ...


class PandasWhen(
    EagerWhen[PandasLikeDataFrame, PandasLikeSeries, PandasLikeExpr, NativeSeriesT]
):
    @property
    # Signature of "_then" incompatible with supertype "CompliantWhen"
    # ArrowWhen seems to follow the same pattern, but no mypy complaint there?
    def _then(self) -> type[PandasThen]:  # type: ignore[override]
        return PandasThen

    def _if_then_else(
        self,
        when: NativeSeriesT,
        then: NativeSeriesT,
        otherwise: NativeSeriesT | NonNestedLiteral,
    ) -> NativeSeriesT:
        where: Incomplete = then.where
        return where(when) if otherwise is None else where(when, otherwise)


class PandasThen(
    CompliantThen[PandasLikeDataFrame, PandasLikeSeries, PandasLikeExpr, PandasWhen],
    PandasLikeExpr,
):
    _depth: int = 0
    _scalar_kwargs: ScalarKwargs = {}  # noqa: RUF012
    _function_name: str = "whenthen"
