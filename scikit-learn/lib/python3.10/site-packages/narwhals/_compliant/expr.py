from __future__ import annotations

from collections.abc import Mapping
from functools import partial
from operator import methodcaller
from typing import TYPE_CHECKING, Any, Callable, Generic, Literal, Protocol

from narwhals._compliant.any_namespace import (
    CatNamespace,
    DateTimeNamespace,
    ListNamespace,
    NameNamespace,
    StringNamespace,
    StructNamespace,
)
from narwhals._compliant.column import CompliantColumn
from narwhals._compliant.namespace import CompliantNamespace
from narwhals._compliant.typing import (
    AliasName,
    AliasNames,
    CompliantExprT_co,
    CompliantFrameT,
    CompliantLazyFrameT,
    CompliantSeriesOrNativeExprT_co,
    EagerDataFrameT,
    EagerExprT,
    EagerSeriesT,
    LazyExprT,
    NativeExprT,
)
from narwhals._utils import (
    _StoresCompliant,
    not_implemented,
    qualified_type_name,
    zip_strict,
)
from narwhals.dependencies import is_numpy_array, is_numpy_scalar

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from typing_extensions import Self, TypeIs

    from narwhals._compliant.namespace import CompliantNamespace, EagerNamespace
    from narwhals._compliant.series import CompliantSeries
    from narwhals._compliant.typing import AliasNames, EvalNames, EvalSeries, ScalarKwargs
    from narwhals._expression_parsing import ExprKind, ExprMetadata
    from narwhals._utils import Implementation, Version, _LimitedContext
    from narwhals.typing import (
        ClosedInterval,
        FillNullStrategy,
        IntoDType,
        ModeKeepStrategy,
        NonNestedLiteral,
        NumericLiteral,
        RankMethod,
        RollingInterpolationMethod,
        TemporalLiteral,
        TimeUnit,
    )

__all__ = ["CompliantExpr", "DepthTrackingExpr", "EagerExpr", "LazyExpr", "NativeExpr"]


class NativeExpr(Protocol):
    """An `Expr`-like object from a package with [Lazy-only support](https://narwhals-dev.github.io/narwhals/extending/#levels-of-support).

    Protocol members are chosen *purely* for matching statically - as they
    are common to all currently supported packages.
    """

    def between(self, *args: Any, **kwds: Any) -> Any: ...

    # NOTE: None of these are annotated for `dx.Series`, but are added imperatively
    # Probably better to define a sub-protocol for `NativeSQLExpr`
    # - match `dx.Series` to `NativeExpr`
    # - match the others to `NativeSQLExpr`
    def __gt__(self, value: Any, /) -> Self: ...
    def __lt__(self, value: Any, /) -> Self: ...
    def __ge__(self, value: Any, /) -> Self: ...
    def __le__(self, value: Any, /) -> Self: ...
    def __eq__(self, value: Any, /) -> Self: ...  # type: ignore[override]
    def __ne__(self, value: Any, /) -> Self: ...  # type: ignore[override]


class CompliantExpr(
    CompliantColumn, Protocol[CompliantFrameT, CompliantSeriesOrNativeExprT_co]
):
    # NOTE: `narwhals`
    _implementation: Implementation
    _evaluate_output_names: EvalNames[CompliantFrameT]
    _alias_output_names: AliasNames | None
    _metadata: ExprMetadata | None

    def __call__(
        self, df: CompliantFrameT
    ) -> Sequence[CompliantSeriesOrNativeExprT_co]: ...
    def __narwhals_expr__(self) -> Self:  # pragma: no cover
        return self

    def __narwhals_namespace__(self) -> CompliantNamespace[CompliantFrameT, Self]: ...
    @classmethod
    def from_column_indices(
        cls, *column_indices: int, context: _LimitedContext
    ) -> Self: ...
    @classmethod
    def from_column_names(
        cls,
        evaluate_column_names: EvalNames[CompliantFrameT],
        /,
        *,
        context: _LimitedContext,
    ) -> Self: ...
    def broadcast(
        self, kind: Literal[ExprKind.AGGREGATION, ExprKind.LITERAL]
    ) -> Self: ...

    # NOTE: `polars`
    def all(self) -> Self: ...
    def any(self) -> Self: ...
    def count(self) -> Self: ...
    def min(self) -> Self: ...
    def max(self) -> Self: ...
    def mean(self) -> Self: ...
    def sum(self) -> Self: ...
    def median(self) -> Self: ...
    def first(self) -> Self: ...
    def last(self) -> Self: ...
    def skew(self) -> Self: ...
    def kurtosis(self) -> Self: ...
    def std(self, *, ddof: int) -> Self: ...
    def var(self, *, ddof: int) -> Self: ...
    def n_unique(self) -> Self: ...
    def null_count(self) -> Self: ...
    def len(self) -> Self: ...
    def over(self, partition_by: Sequence[str], order_by: Sequence[str]) -> Self: ...
    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> Self: ...
    def map_batches(
        self,
        function: Callable[[CompliantSeries[Any]], CompliantExpr[Any, Any]],
        return_dtype: IntoDType | None,
        *,
        returns_scalar: bool,
    ) -> Self: ...
    @property
    def name(self) -> NameNamespace[Self]: ...


class ImplExpr(
    CompliantExpr[CompliantFrameT, CompliantSeriesOrNativeExprT_co],
    Protocol[CompliantFrameT, CompliantSeriesOrNativeExprT_co],
):
    @staticmethod
    def _eval_names_indices(indices: Sequence[int], /) -> EvalNames[CompliantFrameT]:
        def fn(df: CompliantFrameT) -> Sequence[str]:
            column_names = df.columns
            return [column_names[i] for i in indices]

        return fn

    def _evaluate_aliases(self, frame: CompliantFrameT, /) -> Sequence[str]:
        # NOTE: Ignore intermittent [False Negative]
        # Argument of type "CompliantFrameT@ImplExpr" cannot be assigned to parameter of type "CompliantFrameT@ImplExpr"
        #  Type "CompliantFrameT@ImplExpr" is not assignable to type "CompliantFrameT@ImplExpr"
        names = self._evaluate_output_names(frame)  # pyright: ignore[reportArgumentType]
        return alias(names) if (alias := self._alias_output_names) else names


class DepthTrackingExpr(
    ImplExpr[CompliantFrameT, CompliantSeriesOrNativeExprT_co],
    Protocol[CompliantFrameT, CompliantSeriesOrNativeExprT_co],
):
    _depth: int
    _function_name: str

    # NOTE: pyright bug?
    # Method "from_column_names" overrides class "CompliantExpr" in an incompatible manner
    # Parameter 2 type mismatch: base parameter is type "EvalNames[CompliantFrameT@DepthTrackingExpr]", override parameter is type "EvalNames[CompliantFrameT@DepthTrackingExpr]"
    #   Type "EvalNames[CompliantFrameT@DepthTrackingExpr]" is not assignable to type "EvalNames[CompliantFrameT@DepthTrackingExpr]"
    #     Parameter 1: type "CompliantFrameT@DepthTrackingExpr" is incompatible with type "CompliantFrameT@DepthTrackingExpr"
    #       Type "CompliantFrameT@DepthTrackingExpr" is not assignable to type "CompliantFrameT@DepthTrackingExpr"
    @classmethod
    def from_column_names(  # pyright: ignore[reportIncompatibleMethodOverride]
        cls: type[Self],
        evaluate_column_names: EvalNames[CompliantFrameT],
        /,
        *,
        context: _LimitedContext,
        function_name: str = "",
    ) -> Self: ...

    def _is_elementary(self) -> bool:
        """Check if expr is elementary.

        Examples:
            - nw.col('a').mean()  # depth 1
            - nw.mean('a')  # depth 1
            - nw.len()  # depth 0

        as opposed to, say

            - nw.col('a').filter(nw.col('b')>nw.col('c')).max()

        Elementary expressions are the only ones supported properly in
        pandas, PyArrow, and Dask.
        """
        return self._depth < 2

    def __repr__(self) -> str:  # pragma: no cover
        return f"{type(self).__name__}(depth={self._depth}, function_name={self._function_name})"


class EagerExpr(
    DepthTrackingExpr[EagerDataFrameT, EagerSeriesT],
    Protocol[EagerDataFrameT, EagerSeriesT],
):
    _call: EvalSeries[EagerDataFrameT, EagerSeriesT]
    _scalar_kwargs: ScalarKwargs

    def __init__(
        self,
        call: EvalSeries[EagerDataFrameT, EagerSeriesT],
        *,
        depth: int,
        function_name: str,
        evaluate_output_names: EvalNames[EagerDataFrameT],
        alias_output_names: AliasNames | None,
        implementation: Implementation,
        version: Version,
        scalar_kwargs: ScalarKwargs | None = None,
    ) -> None: ...

    def __call__(self, df: EagerDataFrameT) -> Sequence[EagerSeriesT]:
        return self._call(df)

    def __narwhals_namespace__(
        self,
    ) -> EagerNamespace[EagerDataFrameT, EagerSeriesT, Self, Any, Any]: ...
    @classmethod
    def _from_callable(
        cls,
        func: EvalSeries[EagerDataFrameT, EagerSeriesT],
        *,
        depth: int,
        function_name: str,
        evaluate_output_names: EvalNames[EagerDataFrameT],
        alias_output_names: AliasNames | None,
        context: _LimitedContext,
        scalar_kwargs: ScalarKwargs | None = None,
    ) -> Self:
        return cls(
            func,
            depth=depth,
            function_name=function_name,
            evaluate_output_names=evaluate_output_names,
            alias_output_names=alias_output_names,
            implementation=context._implementation,
            version=context._version,
            scalar_kwargs=scalar_kwargs,
        )

    @classmethod
    def _from_series(cls, series: EagerSeriesT) -> Self:
        return cls(
            lambda _df: [series],
            depth=0,
            function_name="series",
            evaluate_output_names=lambda _df: [series.name],
            alias_output_names=None,
            implementation=series._implementation,
            version=series._version,
        )

    def _with_alias_output_names(self, alias_name: AliasName | None, /) -> Self:
        current_alias_output_names = self._alias_output_names
        alias_output_names: AliasNames | None = (
            None
            if alias_name is None
            else (
                lambda output_names: [
                    alias_name(x) for x in current_alias_output_names(output_names)
                ]
            )
            if current_alias_output_names is not None
            else (lambda output_names: [alias_name(x) for x in output_names])
        )

        def func(df: EagerDataFrameT) -> list[EagerSeriesT]:
            if alias_output_names:
                return [
                    series.alias(name)
                    for series, name in zip_strict(
                        self(df), alias_output_names(self._evaluate_output_names(df))
                    )
                ]
            return [
                series.alias(name)
                for series, name in zip_strict(self(df), self._evaluate_output_names(df))
            ]

        return self.__class__(
            func,
            depth=self._depth,
            function_name=self._function_name,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=alias_output_names,
            implementation=self._implementation,
            version=self._version,
            scalar_kwargs=self._scalar_kwargs,
        )

    def _reuse_series(
        self,
        method_name: str,
        *,
        returns_scalar: bool = False,
        scalar_kwargs: ScalarKwargs | None = None,
        **expressifiable_args: Any,
    ) -> Self:
        """Reuse Series implementation for expression.

        If Series.foo is already defined, and we'd like Expr.foo to be the same, we can
        leverage this method to do that for us.

        Arguments:
            method_name: name of method.
            returns_scalar: whether the Series version returns a scalar. In this case,
                the expression version should return a 1-row Series.
            scalar_kwargs: non-expressifiable args which we may need to reuse in `agg` or `over`,
                such as `ddof` for `std` and `var`.
            expressifiable_args: keyword arguments to pass to function, which may
                be expressifiable (e.g. `nw.col('a').is_between(3, nw.col('b')))`).
        """
        func = partial(
            self._reuse_series_inner,
            method_name=method_name,
            returns_scalar=returns_scalar,
            scalar_kwargs=scalar_kwargs or {},
            expressifiable_args=expressifiable_args,
        )
        return self._from_callable(
            func,
            depth=self._depth + 1,
            function_name=f"{self._function_name}->{method_name}",
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            scalar_kwargs=scalar_kwargs,
            context=self,
        )

    # For PyArrow.Series, we return Python Scalars (like Polars does) instead of PyArrow Scalars.
    # However, when working with expressions, we keep everything PyArrow-native.
    def _reuse_series_extra_kwargs(
        self, *, returns_scalar: bool = False
    ) -> dict[str, Any]:
        return {}

    @classmethod
    def _is_expr(cls, obj: Self | Any) -> TypeIs[Self]:
        return hasattr(obj, "__narwhals_expr__")

    def _reuse_series_inner(
        self,
        df: EagerDataFrameT,
        *,
        method_name: str,
        returns_scalar: bool,
        scalar_kwargs: ScalarKwargs,
        expressifiable_args: dict[str, Any],
    ) -> Sequence[EagerSeriesT]:
        kwargs = {
            **scalar_kwargs,
            **{
                name: df._evaluate_expr(value) if self._is_expr(value) else value
                for name, value in expressifiable_args.items()
            },
        }
        method = methodcaller(
            method_name,
            **self._reuse_series_extra_kwargs(returns_scalar=returns_scalar),
            **kwargs,
        )
        out: Sequence[EagerSeriesT] = [
            series._from_scalar(method(series)) if returns_scalar else method(series)
            for series in self(df)
        ]
        aliases, names = self._evaluate_aliases(df), (s.name for s in out)
        if any(
            alias != name for alias, name in zip_strict(aliases, names)
        ):  # pragma: no cover
            msg = (
                f"Safety assertion failed, please report a bug to https://github.com/narwhals-dev/narwhals/issues\n"
                f"Expression aliases: {aliases}\n"
            )
            raise AssertionError(msg)
        return out

    def _reuse_series_namespace(
        self,
        series_namespace: Literal["cat", "dt", "list", "name", "str", "struct"],
        method_name: str,
        **expressifiable_args: Any,
    ) -> Self:
        """Reuse Series implementation for expression.

        Just like `_reuse_series`, but for e.g. `Expr.dt.foo` instead
        of `Expr.foo`.

        Arguments:
            series_namespace: The Series namespace.
            method_name: name of method, within `series_namespace`.
            expressifiable_args: keyword arguments to pass to function, which may
                be expressifiable (e.g. `nw.col('a').str.replace('abc', nw.col('b')))`).
        """

        def inner(df: EagerDataFrameT) -> list[EagerSeriesT]:
            kwargs = {
                name: df._evaluate_expr(value) if self._is_expr(value) else value
                for name, value in expressifiable_args.items()
            }
            return [
                getattr(getattr(series, series_namespace), method_name)(**kwargs)
                for series in self(df)
            ]

        return self._from_callable(
            inner,
            depth=self._depth + 1,
            function_name=f"{self._function_name}->{series_namespace}.{method_name}",
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            scalar_kwargs=self._scalar_kwargs,
            context=self,
        )

    def broadcast(self, kind: Literal[ExprKind.AGGREGATION, ExprKind.LITERAL]) -> Self:
        # Mark the resulting Series with `_broadcast = True`.
        # Then, when extracting native objects, `extract_native` will
        # know what to do.
        def func(df: EagerDataFrameT) -> list[EagerSeriesT]:
            results = []
            for result in self(df):
                result._broadcast = True
                results.append(result)
            return results

        return type(self)(
            func,
            depth=self._depth,
            function_name=self._function_name,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            implementation=self._implementation,
            version=self._version,
            scalar_kwargs=self._scalar_kwargs,
        )

    def cast(self, dtype: IntoDType) -> Self:
        return self._reuse_series("cast", dtype=dtype)

    def _with_binary(self, operator: str, other: Self | Any, /) -> Self:
        return self._reuse_series(operator, other=other)

    def _with_binary_right(self, operator: str, other: Self | Any, /) -> Self:
        return self.alias("literal")._reuse_series(operator, other=other)

    def __eq__(self, other: Self | Any) -> Self:  # type: ignore[override]
        return self._with_binary("__eq__", other)

    def __ne__(self, other: Self | Any) -> Self:  # type: ignore[override]
        return self._with_binary("__ne__", other)

    def __ge__(self, other: Self | Any) -> Self:
        return self._with_binary("__ge__", other)

    def __gt__(self, other: Self | Any) -> Self:
        return self._with_binary("__gt__", other)

    def __le__(self, other: Self | Any) -> Self:
        return self._with_binary("__le__", other)

    def __lt__(self, other: Self | Any) -> Self:
        return self._with_binary("__lt__", other)

    def __and__(self, other: Self | bool | Any) -> Self:
        return self._with_binary("__and__", other)

    def __or__(self, other: Self | bool | Any) -> Self:
        return self._with_binary("__or__", other)

    def __add__(self, other: Self | Any) -> Self:
        return self._with_binary("__add__", other)

    def __sub__(self, other: Self | Any) -> Self:
        return self._with_binary("__sub__", other)

    def __rsub__(self, other: Self | Any) -> Self:
        return self._with_binary_right("__rsub__", other)

    def __mul__(self, other: Self | Any) -> Self:
        return self._with_binary("__mul__", other)

    def __truediv__(self, other: Self | Any) -> Self:
        return self._with_binary("__truediv__", other)

    def __rtruediv__(self, other: Self | Any) -> Self:
        return self._with_binary_right("__rtruediv__", other)

    def __floordiv__(self, other: Self | Any) -> Self:
        return self._with_binary("__floordiv__", other)

    def __rfloordiv__(self, other: Self | Any) -> Self:
        return self._with_binary_right("__rfloordiv__", other)

    def __pow__(self, other: Self | Any) -> Self:
        return self._with_binary("__pow__", other)

    def __rpow__(self, other: Self | Any) -> Self:
        return self._with_binary_right("__rpow__", other)

    def __mod__(self, other: Self | Any) -> Self:
        return self._with_binary("__mod__", other)

    def __rmod__(self, other: Self | Any) -> Self:
        return self._with_binary_right("__rmod__", other)

    # Unary
    def __invert__(self) -> Self:
        return self._reuse_series("__invert__")

    # Reductions
    def null_count(self) -> Self:
        return self._reuse_series("null_count", returns_scalar=True)

    def n_unique(self) -> Self:
        return self._reuse_series("n_unique", returns_scalar=True)

    def sum(self) -> Self:
        return self._reuse_series("sum", returns_scalar=True)

    def count(self) -> Self:
        return self._reuse_series("count", returns_scalar=True)

    def mean(self) -> Self:
        return self._reuse_series("mean", returns_scalar=True)

    def median(self) -> Self:
        return self._reuse_series("median", returns_scalar=True)

    def std(self, *, ddof: int) -> Self:
        return self._reuse_series(
            "std", returns_scalar=True, scalar_kwargs={"ddof": ddof}
        )

    def var(self, *, ddof: int) -> Self:
        return self._reuse_series(
            "var", returns_scalar=True, scalar_kwargs={"ddof": ddof}
        )

    def skew(self) -> Self:
        return self._reuse_series("skew", returns_scalar=True)

    def kurtosis(self) -> Self:
        return self._reuse_series("kurtosis", returns_scalar=True)

    def any(self) -> Self:
        return self._reuse_series("any", returns_scalar=True)

    def all(self) -> Self:
        return self._reuse_series("all", returns_scalar=True)

    def max(self) -> Self:
        return self._reuse_series("max", returns_scalar=True)

    def min(self) -> Self:
        return self._reuse_series("min", returns_scalar=True)

    def arg_min(self) -> Self:
        return self._reuse_series("arg_min", returns_scalar=True)

    def arg_max(self) -> Self:
        return self._reuse_series("arg_max", returns_scalar=True)

    # Other

    def clip(
        self,
        lower_bound: Self | NumericLiteral | TemporalLiteral | None,
        upper_bound: Self | NumericLiteral | TemporalLiteral | None,
    ) -> Self:
        return self._reuse_series(
            "clip", lower_bound=lower_bound, upper_bound=upper_bound
        )

    def is_null(self) -> Self:
        return self._reuse_series("is_null")

    def is_nan(self) -> Self:
        return self._reuse_series("is_nan")

    def fill_nan(self, value: float | None) -> Self:
        return self._reuse_series("fill_nan", value=value)

    def fill_null(
        self,
        value: Self | NonNestedLiteral,
        strategy: FillNullStrategy | None,
        limit: int | None,
    ) -> Self:
        return self._reuse_series(
            "fill_null", value=value, scalar_kwargs={"strategy": strategy, "limit": limit}
        )

    def is_in(self, other: Any) -> Self:
        return self._reuse_series("is_in", other=other)

    def arg_true(self) -> Self:
        return self._reuse_series("arg_true")

    def filter(self, *predicates: Self) -> Self:
        plx = self.__narwhals_namespace__()
        predicate = plx.all_horizontal(*predicates, ignore_nulls=False)
        return self._reuse_series("filter", predicate=predicate)

    def drop_nulls(self) -> Self:
        return self._reuse_series("drop_nulls")

    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> Self:
        return self._reuse_series(
            "replace_strict", old=old, new=new, return_dtype=return_dtype
        )

    def sort(self, *, descending: bool, nulls_last: bool) -> Self:
        return self._reuse_series("sort", descending=descending, nulls_last=nulls_last)

    def abs(self) -> Self:
        return self._reuse_series("abs")

    def unique(self) -> Self:
        return self._reuse_series("unique", maintain_order=False)

    def diff(self) -> Self:
        return self._reuse_series("diff")

    def sample(
        self,
        n: int | None,
        *,
        fraction: float | None,
        with_replacement: bool,
        seed: int | None,
    ) -> Self:
        return self._reuse_series(
            "sample", n=n, fraction=fraction, with_replacement=with_replacement, seed=seed
        )

    def alias(self, name: str) -> Self:
        def alias_output_names(names: Sequence[str]) -> Sequence[str]:
            if len(names) != 1:
                msg = f"Expected function with single output, found output names: {names}"
                raise ValueError(msg)
            return [name]

        # Define this one manually, so that we can
        # override `output_names` and not increase depth
        return type(self)(
            lambda df: [series.alias(name) for series in self(df)],
            depth=self._depth,
            function_name=self._function_name,
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=alias_output_names,
            implementation=self._implementation,
            version=self._version,
            scalar_kwargs=self._scalar_kwargs,
        )

    def is_unique(self) -> Self:
        return self._reuse_series("is_unique")

    def is_first_distinct(self) -> Self:
        return self._reuse_series("is_first_distinct")

    def is_last_distinct(self) -> Self:
        return self._reuse_series("is_last_distinct")

    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> Self:
        return self._reuse_series(
            "quantile",
            returns_scalar=True,
            scalar_kwargs={"quantile": quantile, "interpolation": interpolation},
        )

    def head(self, n: int) -> Self:
        return self._reuse_series("head", scalar_kwargs={"n": n})

    def tail(self, n: int) -> Self:
        return self._reuse_series("tail", scalar_kwargs={"n": n})

    def round(self, decimals: int) -> Self:
        return self._reuse_series("round", decimals=decimals)

    def floor(self) -> Self:
        return self._reuse_series("floor")

    def ceil(self) -> Self:
        return self._reuse_series("ceil")

    def len(self) -> Self:
        return self._reuse_series("len", returns_scalar=True)

    def gather_every(self, n: int, offset: int) -> Self:
        return self._reuse_series("gather_every", n=n, offset=offset)

    def mode(self, *, keep: ModeKeepStrategy) -> Self:
        return self._reuse_series("mode", scalar_kwargs={"keep": keep})

    def is_finite(self) -> Self:
        return self._reuse_series("is_finite")

    def rolling_mean(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        return self._reuse_series(
            "rolling_mean",
            scalar_kwargs={
                "window_size": window_size,
                "min_samples": min_samples,
                "center": center,
            },
        )

    def rolling_std(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        return self._reuse_series(
            "rolling_std",
            scalar_kwargs={
                "window_size": window_size,
                "min_samples": min_samples,
                "center": center,
                "ddof": ddof,
            },
        )

    def rolling_sum(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        return self._reuse_series(
            "rolling_sum",
            scalar_kwargs={
                "window_size": window_size,
                "min_samples": min_samples,
                "center": center,
            },
        )

    def rolling_var(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        return self._reuse_series(
            "rolling_var",
            scalar_kwargs={
                "window_size": window_size,
                "min_samples": min_samples,
                "center": center,
                "ddof": ddof,
            },
        )

    def map_batches(
        self,
        function: Callable[[Any], Any],
        return_dtype: IntoDType | None,
        *,
        returns_scalar: bool,
    ) -> Self:
        def func(df: EagerDataFrameT) -> Sequence[EagerSeriesT]:
            udf_series_in = self(df)
            output_names = (input_series.name for input_series in udf_series_in)
            udf_series_out = tuple(function(series) for series in udf_series_in)
            _first_in, _first_out = udf_series_in[0], udf_series_out[0]

            result: Sequence[EagerSeriesT]
            it = zip_strict(udf_series_out, output_names)
            if is_numpy_array(_first_out) or is_numpy_scalar(_first_out):
                from_numpy = partial(_first_in.from_numpy, context=self)
                result = tuple(from_numpy(arr).alias(out_name) for arr, out_name in it)
            elif isinstance(_first_out, _first_in.__class__):  # compliant series
                result = tuple(series.alias(out_name) for series, out_name in it)
            else:  # If everything else fails, assume scalar case
                from_scalar = _first_in._from_scalar
                result = tuple(from_scalar(val).alias(out_name) for val, out_name in it)

            if return_dtype is not None:
                result = tuple(series.cast(return_dtype) for series in result)

            is_scalar_result = tuple(len(r) == 1 for r in result)
            if (not returns_scalar) and any(is_scalar_result) and (len(df) > 1):
                _idx = is_scalar_result.index(True)  # Index of first result with length 1
                _type = type(udf_series_out[_idx])
                msg = (
                    "`map_batches` with `returns_scalar=False` must return a Series; "
                    f"found '{qualified_type_name(_type)}'.\n\nIf `returns_scalar` "
                    "is set to `True`, a returned value can be a scalar value."
                )
                raise TypeError(msg)
            return result

        return self._from_callable(
            func,
            depth=self._depth + 1,
            function_name=self._function_name + "->map_batches",
            evaluate_output_names=self._evaluate_output_names,
            alias_output_names=self._alias_output_names,
            context=self,
        )

    def shift(self, n: int) -> Self:
        return self._reuse_series("shift", scalar_kwargs={"n": n})

    def cum_sum(self, *, reverse: bool) -> Self:
        return self._reuse_series("cum_sum", scalar_kwargs={"reverse": reverse})

    def cum_count(self, *, reverse: bool) -> Self:
        return self._reuse_series("cum_count", scalar_kwargs={"reverse": reverse})

    def cum_min(self, *, reverse: bool) -> Self:
        return self._reuse_series("cum_min", scalar_kwargs={"reverse": reverse})

    def cum_max(self, *, reverse: bool) -> Self:
        return self._reuse_series("cum_max", scalar_kwargs={"reverse": reverse})

    def cum_prod(self, *, reverse: bool) -> Self:
        return self._reuse_series("cum_prod", scalar_kwargs={"reverse": reverse})

    def rank(self, method: RankMethod, *, descending: bool) -> Self:
        return self._reuse_series(
            "rank", scalar_kwargs={"method": method, "descending": descending}
        )

    def log(self, base: float) -> Self:
        return self._reuse_series("log", base=base)

    def exp(self) -> Self:
        return self._reuse_series("exp")

    def sqrt(self) -> Self:
        return self._reuse_series("sqrt")

    def is_between(
        self, lower_bound: Any, upper_bound: Any, closed: ClosedInterval
    ) -> Self:
        return self._reuse_series(
            "is_between", lower_bound=lower_bound, upper_bound=upper_bound, closed=closed
        )

    def is_close(
        self,
        other: Self | NumericLiteral,
        *,
        abs_tol: float,
        rel_tol: float,
        nans_equal: bool,
    ) -> Self:
        return self._reuse_series(
            "is_close",
            other=other,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
            nans_equal=nans_equal,
        )

    def first(self) -> Self:
        return self._reuse_series("first", returns_scalar=True)

    def last(self) -> Self:
        return self._reuse_series("last", returns_scalar=True)

    @property
    def cat(self) -> EagerExprCatNamespace[Self]:
        return EagerExprCatNamespace(self)

    @property
    def dt(self) -> EagerExprDateTimeNamespace[Self]:
        return EagerExprDateTimeNamespace(self)

    @property
    def list(self) -> EagerExprListNamespace[Self]:
        return EagerExprListNamespace(self)

    @property
    def name(self) -> EagerExprNameNamespace[Self]:
        return EagerExprNameNamespace(self)

    @property
    def str(self) -> EagerExprStringNamespace[Self]:
        return EagerExprStringNamespace(self)

    @property
    def struct(self) -> EagerExprStructNamespace[Self]:
        return EagerExprStructNamespace(self)


# mypy thinks `NativeExprT` should be covariant, pyright thinks it should be invariant
class LazyExpr(  # type: ignore[misc]
    ImplExpr[CompliantLazyFrameT, NativeExprT], Protocol[CompliantLazyFrameT, NativeExprT]
):
    def _with_alias_output_names(self, func: AliasNames | None, /) -> Self: ...
    def alias(self, name: str) -> Self:
        def fn(names: Sequence[str]) -> Sequence[str]:
            if len(names) != 1:
                msg = f"Expected function with single output, found output names: {names}"
                raise ValueError(msg)
            return [name]

        return self._with_alias_output_names(fn)

    @property
    def name(self) -> LazyExprNameNamespace[Self]:
        return LazyExprNameNamespace(self)

    ewm_mean = not_implemented()  # type: ignore[misc]
    map_batches = not_implemented()  # type: ignore[misc]
    replace_strict = not_implemented()  # type: ignore[misc]

    cat: not_implemented = not_implemented()  # type: ignore[assignment]


class _ExprNamespace(  # type: ignore[misc]
    _StoresCompliant[CompliantExprT_co], Protocol[CompliantExprT_co]
):
    _compliant_expr: CompliantExprT_co

    @property
    def compliant(self) -> CompliantExprT_co:
        return self._compliant_expr


class EagerExprNamespace(_ExprNamespace[EagerExprT], Generic[EagerExprT]):
    def __init__(self, expr: EagerExprT, /) -> None:
        self._compliant_expr = expr


class LazyExprNamespace(_ExprNamespace[LazyExprT], Generic[LazyExprT]):
    def __init__(self, expr: LazyExprT, /) -> None:
        self._compliant_expr = expr


class EagerExprCatNamespace(
    EagerExprNamespace[EagerExprT], CatNamespace[EagerExprT], Generic[EagerExprT]
):
    def get_categories(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("cat", "get_categories")


class EagerExprDateTimeNamespace(
    EagerExprNamespace[EagerExprT], DateTimeNamespace[EagerExprT], Generic[EagerExprT]
):
    def to_string(self, format: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "to_string", format=format)

    def replace_time_zone(self, time_zone: str | None) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "dt", "replace_time_zone", time_zone=time_zone
        )

    def convert_time_zone(self, time_zone: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "dt", "convert_time_zone", time_zone=time_zone
        )

    def timestamp(self, time_unit: TimeUnit) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "dt", "timestamp", time_unit=time_unit
        )

    def date(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "date")

    def year(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "year")

    def month(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "month")

    def day(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "day")

    def hour(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "hour")

    def minute(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "minute")

    def second(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "second")

    def millisecond(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "millisecond")

    def microsecond(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "microsecond")

    def nanosecond(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "nanosecond")

    def ordinal_day(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "ordinal_day")

    def weekday(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "weekday")

    def total_minutes(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "total_minutes")

    def total_seconds(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "total_seconds")

    def total_milliseconds(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "total_milliseconds")

    def total_microseconds(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "total_microseconds")

    def total_nanoseconds(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "total_nanoseconds")

    def truncate(self, every: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "truncate", every=every)

    def offset_by(self, by: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("dt", "offset_by", by=by)


class EagerExprListNamespace(
    EagerExprNamespace[EagerExprT], ListNamespace[EagerExprT], Generic[EagerExprT]
):
    def len(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("list", "len")

    def unique(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("list", "unique")

    def contains(self, item: NonNestedLiteral) -> EagerExprT:
        return self.compliant._reuse_series_namespace("list", "contains", item=item)

    def get(self, index: int) -> EagerExprT:
        return self.compliant._reuse_series_namespace("list", "get", index=index)


class CompliantExprNameNamespace(  # type: ignore[misc]
    _ExprNamespace[CompliantExprT_co],
    NameNamespace[CompliantExprT_co],
    Protocol[CompliantExprT_co],
):
    def keep(self) -> CompliantExprT_co:
        return self._from_callable(None)

    def map(self, function: AliasName) -> CompliantExprT_co:
        return self._from_callable(function)

    def prefix(self, prefix: str) -> CompliantExprT_co:
        return self._from_callable(lambda name: f"{prefix}{name}")

    def suffix(self, suffix: str) -> CompliantExprT_co:
        return self._from_callable(lambda name: f"{name}{suffix}")

    def to_lowercase(self) -> CompliantExprT_co:
        return self._from_callable(str.lower)

    def to_uppercase(self) -> CompliantExprT_co:
        return self._from_callable(str.upper)

    @staticmethod
    def _alias_output_names(func: AliasName, /) -> AliasNames:
        def fn(output_names: Sequence[str], /) -> Sequence[str]:
            return [func(name) for name in output_names]

        return fn

    def _from_callable(self, func: AliasName | None, /) -> CompliantExprT_co: ...


class EagerExprNameNamespace(
    EagerExprNamespace[EagerExprT],
    CompliantExprNameNamespace[EagerExprT],
    Generic[EagerExprT],
):
    def _from_callable(self, func: AliasName | None) -> EagerExprT:
        expr = self.compliant
        return expr._with_alias_output_names(func)


class LazyExprNameNamespace(
    LazyExprNamespace[LazyExprT],
    CompliantExprNameNamespace[LazyExprT],
    Generic[LazyExprT],
):
    def _from_callable(self, func: AliasName | None) -> LazyExprT:
        expr = self.compliant
        output_names = self._alias_output_names(func) if func else None
        return expr._with_alias_output_names(output_names)


class EagerExprStringNamespace(
    EagerExprNamespace[EagerExprT], StringNamespace[EagerExprT], Generic[EagerExprT]
):
    def len_chars(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "len_chars")

    def replace(self, pattern: str, value: str, *, literal: bool, n: int) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "str", "replace", pattern=pattern, value=value, literal=literal, n=n
        )

    def replace_all(self, pattern: str, value: str, *, literal: bool) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "str", "replace_all", pattern=pattern, value=value, literal=literal
        )

    def strip_chars(self, characters: str | None) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "str", "strip_chars", characters=characters
        )

    def starts_with(self, prefix: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "starts_with", prefix=prefix)

    def ends_with(self, suffix: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "ends_with", suffix=suffix)

    def contains(self, pattern: str, *, literal: bool) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "str", "contains", pattern=pattern, literal=literal
        )

    def slice(self, offset: int, length: int | None) -> EagerExprT:
        return self.compliant._reuse_series_namespace(
            "str", "slice", offset=offset, length=length
        )

    def split(self, by: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "split", by=by)

    def to_datetime(self, format: str | None) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "to_datetime", format=format)

    def to_date(self, format: str | None) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "to_date", format=format)

    def to_lowercase(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "to_lowercase")

    def to_uppercase(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "to_uppercase")

    def zfill(self, width: int) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "zfill", width=width)

    def to_titlecase(self) -> EagerExprT:
        return self.compliant._reuse_series_namespace("str", "to_titlecase")


class EagerExprStructNamespace(
    EagerExprNamespace[EagerExprT], StructNamespace[EagerExprT], Generic[EagerExprT]
):
    def field(self, name: str) -> EagerExprT:
        return self.compliant._reuse_series_namespace("struct", "field", name=name).alias(
            name
        )
