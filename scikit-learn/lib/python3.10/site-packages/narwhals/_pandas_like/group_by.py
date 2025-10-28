from __future__ import annotations

import warnings
from functools import lru_cache
from itertools import chain
from operator import methodcaller
from typing import TYPE_CHECKING, Any, ClassVar, Literal

from narwhals._compliant import EagerGroupBy
from narwhals._exceptions import issue_warning
from narwhals._expression_parsing import evaluate_output_names_and_aliases
from narwhals._utils import zip_strict
from narwhals.dependencies import is_pandas_like_dataframe

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence

    import pandas as pd
    from pandas.api.typing import DataFrameGroupBy as _NativeGroupBy
    from typing_extensions import TypeAlias, Unpack

    from narwhals._compliant.typing import NarwhalsAggregation, ScalarKwargs
    from narwhals._pandas_like.dataframe import PandasLikeDataFrame
    from narwhals._pandas_like.expr import PandasLikeExpr

    NativeGroupBy: TypeAlias = "_NativeGroupBy[tuple[str, ...], Literal[True]]"

NativeApply: TypeAlias = "Callable[[pd.DataFrame], pd.Series[Any]]"
InefficientNativeAggregation: TypeAlias = Literal["cov", "skew"]
NativeAggregation: TypeAlias = Literal[
    "any",
    "all",
    "count",
    "idxmax",
    "idxmin",
    "max",
    "mean",
    "median",
    "min",
    "mode",
    "nth",
    "nunique",
    "prod",
    "quantile",
    "sem",
    "size",
    "std",
    "sum",
    "var",
    InefficientNativeAggregation,
]
"""https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html#built-in-aggregation-methods"""

_NativeAgg: TypeAlias = "Callable[[Any], pd.DataFrame | pd.Series[Any]]"
"""Equivalent to a partial method call on `DataFrameGroupBy`."""


NonStrHashable: TypeAlias = Any
"""Because `pandas` allows *"names"* like that 😭"""

_REMAP_ORDERED_INDEX: Mapping[NarwhalsAggregation, Literal[0, -1]] = {
    "first": 0,
    "last": -1,
}


@lru_cache(maxsize=32)
def _native_agg(name: NativeAggregation, /, **kwds: Unpack[ScalarKwargs]) -> _NativeAgg:
    if name == "nunique":
        return methodcaller(name, dropna=False)
    if not kwds or kwds.get("ddof") == 1:
        return methodcaller(name)
    return methodcaller(name, **kwds)


class AggExpr:
    """Wrapper storing the intermediate state per-`PandasLikeExpr`.

    There's a lot of edge cases to handle, so aim to evaluate as little
    as possible - and store anything that's needed twice.

    Warning:
        While a `PandasLikeExpr` can be reused - this wrapper is valid **only**
        in a single `.agg(...)` operation.
    """

    expr: PandasLikeExpr
    output_names: Sequence[str]
    aliases: Sequence[str]

    def __init__(self, expr: PandasLikeExpr) -> None:
        self.expr = expr
        self.output_names = ()
        self.aliases = ()
        self._leaf_name: NarwhalsAggregation | Any = ""

    def with_expand_names(self, group_by: PandasLikeGroupBy, /) -> AggExpr:
        """**Mutating operation**.

        Stores the results of `evaluate_output_names_and_aliases`.
        """
        df = group_by.compliant
        exclude = group_by.exclude
        self.output_names, self.aliases = evaluate_output_names_and_aliases(
            self.expr, df, exclude
        )
        return self

    def _getitem_aggs(
        self, group_by: PandasLikeGroupBy, /
    ) -> pd.DataFrame | pd.Series[Any]:
        """Evaluate the wrapped expression as a group_by operation."""
        result: pd.DataFrame | pd.Series[Any]
        names = self.output_names
        if self.is_len() and self.is_top_level_function():
            result = group_by._grouped.size()
        elif self.is_len():
            result_single = group_by._grouped.size()
            ns = group_by.compliant.__narwhals_namespace__()
            result = ns._concat_horizontal(
                [ns.from_native(result_single).alias(name).native for name in names]
            )
        elif self.is_mode():
            compliant = group_by.compliant
            if (keep := self.kwargs.get("keep")) != "any":  # pragma: no cover
                msg = (
                    f"`Expr.mode(keep='{keep}')` is not implemented in group by context for "
                    f"backend {compliant._implementation}\n\n"
                    "Hint: Use `nw.col(...).mode(keep='any')` instead."
                )
                raise NotImplementedError(msg)

            cols = list(names)
            native = compliant.native
            keys, kwargs = group_by._keys, group_by._kwargs

            # Implementation based on the following suggestion:
            # https://github.com/pandas-dev/pandas/issues/19254#issuecomment-778661578
            ns = compliant.__narwhals_namespace__()
            result = ns._concat_horizontal(
                [
                    native.groupby([*keys, col], **kwargs)
                    .size()
                    .sort_values(ascending=False)
                    .reset_index(col)
                    .groupby(keys, **kwargs)[col]
                    .head(1)
                    .sort_index()
                    for col in cols
                ]
            )
        elif self.is_last() or self.is_first():
            result = self.native_agg()(group_by._grouped[[*group_by._keys, *names]])
            result.set_index(group_by._keys, inplace=True)  # noqa: PD002
        else:
            select = names[0] if len(names) == 1 else list(names)
            result = self.native_agg()(group_by._grouped[select])
        if is_pandas_like_dataframe(result):
            result.columns = list(self.aliases)
        else:
            result.name = self.aliases[0]
        return result

    def is_len(self) -> bool:
        return self.leaf_name == "len"

    def is_last(self) -> bool:
        return self.leaf_name == "last"

    def is_first(self) -> bool:
        return self.leaf_name == "first"

    def is_mode(self) -> bool:
        return self.leaf_name == "mode"

    def is_top_level_function(self) -> bool:
        # e.g. `nw.len()`.
        return self.expr._depth == 0

    @property
    def kwargs(self) -> ScalarKwargs:
        return self.expr._scalar_kwargs

    @property
    def leaf_name(self) -> NarwhalsAggregation | Any:
        if name := self._leaf_name:
            return name
        self._leaf_name = PandasLikeGroupBy._leaf_name(self.expr)
        return self._leaf_name

    def native_agg(self) -> _NativeAgg:
        """Return a partial `DataFrameGroupBy` method, missing only `self`."""
        native_name = PandasLikeGroupBy._remap_expr_name(self.leaf_name)
        if self.leaf_name in _REMAP_ORDERED_INDEX:
            return methodcaller("nth", n=_REMAP_ORDERED_INDEX[self.leaf_name])
        return _native_agg(native_name, **self.kwargs)


class PandasLikeGroupBy(
    EagerGroupBy["PandasLikeDataFrame", "PandasLikeExpr", NativeAggregation]
):
    _REMAP_AGGS: ClassVar[Mapping[NarwhalsAggregation, NativeAggregation]] = {
        "sum": "sum",
        "mean": "mean",
        "median": "median",
        "max": "max",
        "min": "min",
        "mode": "mode",
        "std": "std",
        "var": "var",
        "len": "size",
        "n_unique": "nunique",
        "count": "count",
        "quantile": "quantile",
        "all": "all",
        "any": "any",
        "first": "nth",
        "last": "nth",
    }
    _original_columns: tuple[str, ...]
    """Column names *prior* to any aliasing in `ParseKeysGroupBy`."""

    _keys: list[str]
    """Stores the **aliased** version of group keys from `ParseKeysGroupBy`."""

    _output_key_names: list[str]
    """Stores the **original** version of group keys."""

    _kwargs: Mapping[str, bool]
    """Stores keyword arguments for `DataFrame.groupby` other than `by`."""

    @property
    def exclude(self) -> tuple[str, ...]:
        """Group keys to ignore when expanding multi-output aggregations."""
        return self._exclude

    def __init__(
        self,
        df: PandasLikeDataFrame,
        keys: Sequence[PandasLikeExpr] | Sequence[str],
        /,
        *,
        drop_null_keys: bool,
    ) -> None:
        self._original_columns = tuple(df.columns)
        self._drop_null_keys = drop_null_keys
        self._compliant_frame, self._keys, self._output_key_names = self._parse_keys(
            df, keys
        )
        self._exclude: tuple[str, ...] = (*self._keys, *self._output_key_names)
        # Drop index to avoid potential collisions:
        # https://github.com/narwhals-dev/narwhals/issues/1907.
        native = self.compliant.native
        if set(native.index.names).intersection(self.compliant.columns):
            native = native.reset_index(drop=True)

        self._kwargs = {
            "sort": False,
            "as_index": True,
            "dropna": drop_null_keys,
            "observed": True,
        }
        self._grouped: NativeGroupBy = native.groupby(self._keys.copy(), **self._kwargs)

    def agg(self, *exprs: PandasLikeExpr) -> PandasLikeDataFrame:
        all_aggs_are_simple = True
        agg_exprs: list[AggExpr] = []
        for expr in exprs:
            agg_exprs.append(AggExpr(expr).with_expand_names(self))
            if not self._is_simple(expr):
                all_aggs_are_simple = False

        if all_aggs_are_simple:
            result: pd.DataFrame
            if agg_exprs:
                ns = self.compliant.__narwhals_namespace__()
                result = ns._concat_horizontal(self._getitem_aggs(agg_exprs))
            else:
                result = self.compliant.__native_namespace__().DataFrame(
                    list(self._grouped.groups), columns=self._keys
                )
        elif self.compliant.native.empty:
            raise empty_results_error()
        else:
            result = self._apply_aggs(exprs)
        # NOTE: Keep `inplace=True` to avoid making a redundant copy.
        # This may need updating, depending on https://github.com/pandas-dev/pandas/pull/51466/files
        result.reset_index(inplace=True)  # noqa: PD002
        return self._select_results(result, agg_exprs)

    def _select_results(
        self, df: pd.DataFrame, /, agg_exprs: Sequence[AggExpr]
    ) -> PandasLikeDataFrame:
        """Responsible for remapping temp column names back to original.

        See `ParseKeysGroupBy`.
        """
        new_names = chain.from_iterable(e.aliases for e in agg_exprs)
        return (
            self.compliant._with_native(df, validate_column_names=False)
            .simple_select(*self._keys, *new_names)
            .rename(dict(zip(self._keys, self._output_key_names)))
        )

    def _getitem_aggs(
        self, exprs: Iterable[AggExpr], /
    ) -> list[pd.DataFrame | pd.Series[Any]]:
        return [e._getitem_aggs(self) for e in exprs]

    def _apply_aggs(self, exprs: Iterable[PandasLikeExpr]) -> pd.DataFrame:
        """Stub issue for `include_groups` [pandas-dev/pandas-stubs#1270].

        - [User guide] mentions `include_groups` 4 times without deprecation.
        - [`DataFrameGroupBy.apply`] doc says the default value of `True` is deprecated since `2.2.0`.
        - `False` is explicitly the only *non-deprecated* option, but entirely omitted since [pandas-dev/pandas-stubs#1268].

        [pandas-dev/pandas-stubs#1270]: https://github.com/pandas-dev/pandas-stubs/issues/1270
        [User guide]: https://pandas.pydata.org/pandas-docs/stable/user_guide/groupby.html
        [`DataFrameGroupBy.apply`]: https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.core.groupby.DataFrameGroupBy.apply.html
        [pandas-dev/pandas-stubs#1268]: https://github.com/pandas-dev/pandas-stubs/pull/1268
        """
        warn_complex_group_by()
        impl = self.compliant._implementation
        func = self._apply_exprs_function(exprs)
        apply = self._grouped.apply
        if impl.is_pandas() and impl._backend_version() >= (2, 2):
            return apply(func, include_groups=False)  # type: ignore[call-overload]
        return apply(func)  # pragma: no cover

    def _apply_exprs_function(self, exprs: Iterable[PandasLikeExpr]) -> NativeApply:
        ns = self.compliant.__narwhals_namespace__()
        into_series = ns._series.from_iterable

        def fn(df: pd.DataFrame) -> pd.Series[Any]:
            compliant = self.compliant._with_native(df)
            results = (
                (keys.native.iloc[0], keys.name)
                for expr in exprs
                for keys in expr(compliant)
            )
            out_group, out_names = zip_strict(*results) if results else ([], [])
            return into_series(out_group, index=out_names, context=ns).native

        return fn

    def __iter__(self) -> Iterator[tuple[Any, PandasLikeDataFrame]]:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=".*a length 1 tuple will be returned",
                category=FutureWarning,
            )
            with_native = self.compliant._with_native
            for key, group in self._grouped:
                yield (key, with_native(group).simple_select(*self._original_columns))


def empty_results_error() -> ValueError:
    """Don't even attempt this, it's way too inconsistent across pandas versions."""
    msg = (
        "No results for group-by aggregation.\n\n"
        "Hint: you were probably trying to apply a non-elementary aggregation with a "
        "pandas-like API.\n"
        "Please rewrite your query such that group-by aggregations "
        "are elementary. For example, instead of:\n\n"
        "    df.group_by('a').agg(nw.col('b').round(2).mean())\n\n"
        "use:\n\n"
        "    df.with_columns(nw.col('b').round(2)).group_by('a').agg(nw.col('b').mean())\n\n"
    )
    return ValueError(msg)


def warn_complex_group_by() -> None:
    issue_warning(
        "Found complex group-by expression, which can't be expressed efficiently with the "
        "pandas API. If you can, please rewrite your query such that group-by aggregations "
        "are simple (e.g. mean, std, min, max, ...). \n\n"
        "Please see: "
        "https://narwhals-dev.github.io/narwhals/concepts/improve_group_by_operation/",
        UserWarning,
    )
