from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any, Callable, Literal, TypedDict, TypeVar

if TYPE_CHECKING:
    from typing_extensions import TypeAlias

    from narwhals._compliant.dataframe import (
        CompliantDataFrame,
        CompliantLazyFrame,
        EagerDataFrame,
    )
    from narwhals._compliant.expr import (
        CompliantExpr,
        DepthTrackingExpr,
        EagerExpr,
        LazyExpr,
        NativeExpr,
    )
    from narwhals._compliant.namespace import CompliantNamespace, EagerNamespace
    from narwhals._compliant.series import CompliantSeries, EagerSeries
    from narwhals._compliant.window import WindowInputs
    from narwhals.typing import (
        FillNullStrategy,
        NativeFrame,
        NativeSeries,
        RankMethod,
        RollingInterpolationMethod,
    )

    class ScalarKwargs(TypedDict, total=False):
        """Non-expressifiable args which we may need to reuse in `agg` or `over`."""

        adjust: bool
        alpha: float | None
        center: int
        com: float | None
        ddof: int
        descending: bool
        half_life: float | None
        ignore_nulls: bool
        interpolation: RollingInterpolationMethod
        limit: int | None
        method: RankMethod
        min_samples: int
        n: int
        quantile: float
        reverse: bool
        span: float | None
        strategy: FillNullStrategy | None
        window_size: int


__all__ = [
    "AliasName",
    "AliasNames",
    "CompliantDataFrameT",
    "CompliantFrameT",
    "CompliantLazyFrameT",
    "CompliantSeriesT",
    "EvalNames",
    "EvalSeries",
    "IntoCompliantExpr",
    "NarwhalsAggregation",
    "NativeFrameT_co",
    "NativeSeriesT_co",
]
CompliantExprAny: TypeAlias = "CompliantExpr[Any, Any]"
CompliantSeriesAny: TypeAlias = "CompliantSeries[Any]"
CompliantSeriesOrNativeExprAny: TypeAlias = "CompliantSeriesAny | NativeExpr"
CompliantDataFrameAny: TypeAlias = "CompliantDataFrame[Any, Any, Any, Any]"
CompliantLazyFrameAny: TypeAlias = "CompliantLazyFrame[Any, Any, Any]"
CompliantFrameAny: TypeAlias = "CompliantDataFrameAny | CompliantLazyFrameAny"
CompliantNamespaceAny: TypeAlias = "CompliantNamespace[Any, Any]"

DepthTrackingExprAny: TypeAlias = "DepthTrackingExpr[Any, Any]"

EagerDataFrameAny: TypeAlias = "EagerDataFrame[Any, Any, Any, Any]"
EagerSeriesAny: TypeAlias = "EagerSeries[Any]"
EagerExprAny: TypeAlias = "EagerExpr[Any, Any]"
EagerNamespaceAny: TypeAlias = "EagerNamespace[EagerDataFrameAny, EagerSeriesAny, EagerExprAny, NativeFrame, NativeSeries]"

LazyExprAny: TypeAlias = "LazyExpr[Any, Any]"

NativeExprT = TypeVar("NativeExprT", bound="NativeExpr")
NativeExprT_co = TypeVar("NativeExprT_co", bound="NativeExpr", covariant=True)
NativeSeriesT = TypeVar("NativeSeriesT", bound="NativeSeries")
NativeSeriesT_co = TypeVar("NativeSeriesT_co", bound="NativeSeries", covariant=True)
NativeSeriesT_contra = TypeVar(
    "NativeSeriesT_contra", bound="NativeSeries", contravariant=True
)
NativeFrameT = TypeVar("NativeFrameT", bound="NativeFrame")
NativeFrameT_co = TypeVar("NativeFrameT_co", bound="NativeFrame", covariant=True)
NativeFrameT_contra = TypeVar(
    "NativeFrameT_contra", bound="NativeFrame", contravariant=True
)

CompliantExprT = TypeVar("CompliantExprT", bound=CompliantExprAny)
CompliantExprT_co = TypeVar("CompliantExprT_co", bound=CompliantExprAny, covariant=True)
CompliantExprT_contra = TypeVar(
    "CompliantExprT_contra", bound=CompliantExprAny, contravariant=True
)
CompliantSeriesT = TypeVar("CompliantSeriesT", bound=CompliantSeriesAny)
CompliantSeriesT_co = TypeVar(
    "CompliantSeriesT_co", bound=CompliantSeriesAny, covariant=True
)
CompliantSeriesOrNativeExprT = TypeVar(
    "CompliantSeriesOrNativeExprT", bound=CompliantSeriesOrNativeExprAny
)
CompliantSeriesOrNativeExprT_co = TypeVar(
    "CompliantSeriesOrNativeExprT_co",
    bound=CompliantSeriesOrNativeExprAny,
    covariant=True,
)
CompliantFrameT = TypeVar("CompliantFrameT", bound=CompliantFrameAny)
CompliantFrameT_co = TypeVar(
    "CompliantFrameT_co", bound=CompliantFrameAny, covariant=True
)
CompliantDataFrameT = TypeVar("CompliantDataFrameT", bound=CompliantDataFrameAny)
CompliantDataFrameT_co = TypeVar(
    "CompliantDataFrameT_co", bound=CompliantDataFrameAny, covariant=True
)
CompliantLazyFrameT = TypeVar("CompliantLazyFrameT", bound=CompliantLazyFrameAny)
CompliantLazyFrameT_co = TypeVar(
    "CompliantLazyFrameT_co", bound=CompliantLazyFrameAny, covariant=True
)
CompliantNamespaceT = TypeVar("CompliantNamespaceT", bound=CompliantNamespaceAny)
CompliantNamespaceT_co = TypeVar(
    "CompliantNamespaceT_co", bound=CompliantNamespaceAny, covariant=True
)

IntoCompliantExpr: TypeAlias = "CompliantExpr[CompliantFrameT, CompliantSeriesOrNativeExprT_co] | CompliantSeriesOrNativeExprT_co"

DepthTrackingExprT = TypeVar("DepthTrackingExprT", bound=DepthTrackingExprAny)
DepthTrackingExprT_contra = TypeVar(
    "DepthTrackingExprT_contra", bound=DepthTrackingExprAny, contravariant=True
)

EagerExprT = TypeVar("EagerExprT", bound=EagerExprAny)
EagerExprT_contra = TypeVar("EagerExprT_contra", bound=EagerExprAny, contravariant=True)
EagerSeriesT = TypeVar("EagerSeriesT", bound=EagerSeriesAny)
EagerSeriesT_co = TypeVar("EagerSeriesT_co", bound=EagerSeriesAny, covariant=True)

# NOTE: `pyright` gives false (8) positives if this uses `EagerDataFrameAny`?
EagerDataFrameT = TypeVar("EagerDataFrameT", bound="EagerDataFrame[Any, Any, Any, Any]")

LazyExprT = TypeVar("LazyExprT", bound=LazyExprAny)
LazyExprT_contra = TypeVar("LazyExprT_contra", bound=LazyExprAny, contravariant=True)

AliasNames: TypeAlias = Callable[[Sequence[str]], Sequence[str]]
"""A function aliasing a *sequence* of column names."""

AliasName: TypeAlias = Callable[[str], str]
"""A function aliasing a *single* column name."""

EvalSeries: TypeAlias = Callable[
    [CompliantFrameT], Sequence[CompliantSeriesOrNativeExprT]
]
"""A function from a `Frame` to a sequence of `Series`*.

See [underwater unicorn magic](https://narwhals-dev.github.io/narwhals/how_it_works/).
"""

EvalNames: TypeAlias = Callable[[CompliantFrameT], Sequence[str]]
"""A function from a `Frame` to a sequence of columns names *before* any aliasing takes place."""

WindowFunction: TypeAlias = (
    "Callable[[CompliantFrameT, WindowInputs[NativeExprT]], Sequence[NativeExprT]]"
)
"""A function evaluated with `over(partition_by=..., order_by=...)`."""

NarwhalsAggregation: TypeAlias = Literal[
    "sum",
    "mean",
    "median",
    "max",
    "min",
    "std",
    "var",
    "len",
    "n_unique",
    "count",
    "quantile",
    "all",
    "any",
]
"""`Expr` methods we aim to support in `DepthTrackingGroupBy`.

Be sure to update me if you're working on one of these:
- https://github.com/narwhals-dev/narwhals/issues/981
- https://github.com/narwhals-dev/narwhals/issues/2385
- https://github.com/narwhals-dev/narwhals/issues/2484
- https://github.com/narwhals-dev/narwhals/issues/2526
- https://github.com/narwhals-dev/narwhals/issues/2660
"""
