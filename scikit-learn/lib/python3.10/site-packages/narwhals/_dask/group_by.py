from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any, Callable, ClassVar

import dask.dataframe as dd

from narwhals._compliant import DepthTrackingGroupBy
from narwhals._expression_parsing import evaluate_output_names_and_aliases
from narwhals._utils import zip_strict

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    import pandas as pd
    from dask.dataframe.api import GroupBy as _DaskGroupBy
    from pandas.core.groupby import SeriesGroupBy as _PandasSeriesGroupBy
    from typing_extensions import TypeAlias

    from narwhals._compliant.typing import NarwhalsAggregation
    from narwhals._dask.dataframe import DaskLazyFrame
    from narwhals._dask.expr import DaskExpr

    PandasSeriesGroupBy: TypeAlias = _PandasSeriesGroupBy[Any, Any]
    _AggFn: TypeAlias = Callable[..., Any]

else:
    try:
        import dask.dataframe.dask_expr as dx
    except ModuleNotFoundError:
        import dask_expr as dx
    _DaskGroupBy = dx._groupby.GroupBy

Aggregation: TypeAlias = "str | _AggFn"
"""The name of an aggregation function, or the function itself."""


def n_unique() -> dd.Aggregation:
    def chunk(s: PandasSeriesGroupBy) -> pd.Series[Any]:
        return s.nunique(dropna=False)

    def agg(s0: PandasSeriesGroupBy) -> pd.Series[Any]:
        return s0.sum()

    return dd.Aggregation(name="nunique", chunk=chunk, agg=agg)


def _all() -> dd.Aggregation:
    def chunk(s: PandasSeriesGroupBy) -> pd.Series[Any]:
        return s.all(skipna=True)

    def agg(s0: PandasSeriesGroupBy) -> pd.Series[Any]:
        return s0.all(skipna=True)

    return dd.Aggregation(name="all", chunk=chunk, agg=agg)


def _any() -> dd.Aggregation:
    def chunk(s: PandasSeriesGroupBy) -> pd.Series[Any]:
        return s.any(skipna=True)

    def agg(s0: PandasSeriesGroupBy) -> pd.Series[Any]:
        return s0.any(skipna=True)

    return dd.Aggregation(name="any", chunk=chunk, agg=agg)


def var(ddof: int) -> _AggFn:
    return partial(_DaskGroupBy.var, ddof=ddof)


def std(ddof: int) -> _AggFn:
    return partial(_DaskGroupBy.std, ddof=ddof)


class DaskLazyGroupBy(DepthTrackingGroupBy["DaskLazyFrame", "DaskExpr", Aggregation]):
    _REMAP_AGGS: ClassVar[Mapping[NarwhalsAggregation, Aggregation]] = {
        "sum": "sum",
        "mean": "mean",
        "median": "median",
        "max": "max",
        "min": "min",
        "std": std,
        "var": var,
        "len": "size",
        "n_unique": n_unique,
        "count": "count",
        "quantile": "quantile",
        "all": _all,
        "any": _any,
    }

    def __init__(
        self,
        df: DaskLazyFrame,
        keys: Sequence[DaskExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._compliant_frame, self._keys, self._output_key_names = self._parse_keys(
            df, keys=keys
        )
        self._grouped = self.compliant.native.groupby(
            self._keys, dropna=drop_null_keys, observed=True
        )

    def agg(self, *exprs: DaskExpr) -> DaskLazyFrame:
        from narwhals._dask.dataframe import DaskLazyFrame

        if not exprs:
            # No aggregation provided
            return (
                self.compliant.simple_select(*self._keys)
                .unique(self._keys, keep="any", order_by=None)
                .rename(dict(zip(self._keys, self._output_key_names)))
            )

        self._ensure_all_simple(exprs)
        # This should be the fastpath, but cuDF is too far behind to use it.
        # - https://github.com/rapidsai/cudf/issues/15118
        # - https://github.com/rapidsai/cudf/issues/15084
        simple_aggregations: dict[str, tuple[str, Aggregation]] = {}
        exclude = (*self._keys, *self._output_key_names)
        for expr in exprs:
            output_names, aliases = evaluate_output_names_and_aliases(
                expr, self.compliant, exclude
            )
            if expr._depth == 0:
                # e.g. `agg(nw.len())`
                column = self._keys[0]
                agg_fn = self._remap_expr_name(expr._function_name)
                simple_aggregations.update(dict.fromkeys(aliases, (column, agg_fn)))
                continue

            # e.g. `agg(nw.mean('a'))`
            agg_fn = self._remap_expr_name(self._leaf_name(expr))
            # deal with n_unique case in a "lazy" mode to not depend on dask globally
            agg_fn = agg_fn(**expr._scalar_kwargs) if callable(agg_fn) else agg_fn
            simple_aggregations.update(
                (alias, (output_name, agg_fn))
                for alias, output_name in zip_strict(aliases, output_names)
            )
        return DaskLazyFrame(
            self._grouped.agg(**simple_aggregations).reset_index(),
            version=self.compliant._version,
        ).rename(dict(zip(self._keys, self._output_key_names)))
