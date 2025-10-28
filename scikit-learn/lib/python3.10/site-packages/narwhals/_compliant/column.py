from __future__ import annotations

from typing import TYPE_CHECKING, Any, Protocol

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from typing_extensions import Self

    from narwhals._compliant.any_namespace import (
        CatNamespace,
        DateTimeNamespace,
        ListNamespace,
        StringNamespace,
        StructNamespace,
    )
    from narwhals._compliant.namespace import CompliantNamespace
    from narwhals._utils import Version
    from narwhals.typing import (
        ClosedInterval,
        FillNullStrategy,
        IntoDType,
        ModeKeepStrategy,
        NonNestedLiteral,
        NumericLiteral,
        RankMethod,
        TemporalLiteral,
    )

__all__ = ["CompliantColumn"]


class CompliantColumn(Protocol):
    """Common parts of `Expr`, `Series`."""

    _version: Version

    def __add__(self, other: Any) -> Self: ...
    def __and__(self, other: Any) -> Self: ...
    def __eq__(self, other: object) -> Self: ...  # type: ignore[override]
    def __floordiv__(self, other: Any) -> Self: ...
    def __ge__(self, other: Any) -> Self: ...
    def __gt__(self, other: Any) -> Self: ...
    def __invert__(self) -> Self: ...
    def __le__(self, other: Any) -> Self: ...
    def __lt__(self, other: Any) -> Self: ...
    def __mod__(self, other: Any) -> Self: ...
    def __mul__(self, other: Any) -> Self: ...
    def __ne__(self, other: object) -> Self: ...  # type: ignore[override]
    def __or__(self, other: Any) -> Self: ...
    def __pow__(self, other: Any) -> Self: ...
    def __rfloordiv__(self, other: Any) -> Self: ...
    def __rmod__(self, other: Any) -> Self: ...
    def __rpow__(self, other: Any) -> Self: ...
    def __rsub__(self, other: Any) -> Self: ...
    def __rtruediv__(self, other: Any) -> Self: ...
    def __sub__(self, other: Any) -> Self: ...
    def __truediv__(self, other: Any) -> Self: ...

    def __narwhals_namespace__(self) -> CompliantNamespace[Any, Any]: ...

    def abs(self) -> Self: ...
    def alias(self, name: str) -> Self: ...
    def cast(self, dtype: IntoDType) -> Self: ...
    def clip(
        self,
        lower_bound: Self | NumericLiteral | TemporalLiteral | None,
        upper_bound: Self | NumericLiteral | TemporalLiteral | None,
    ) -> Self: ...
    def cum_count(self, *, reverse: bool) -> Self: ...
    def cum_max(self, *, reverse: bool) -> Self: ...
    def cum_min(self, *, reverse: bool) -> Self: ...
    def cum_prod(self, *, reverse: bool) -> Self: ...
    def cum_sum(self, *, reverse: bool) -> Self: ...
    def diff(self) -> Self: ...
    def drop_nulls(self) -> Self: ...
    def ewm_mean(
        self,
        *,
        com: float | None,
        span: float | None,
        half_life: float | None,
        alpha: float | None,
        adjust: bool,
        min_samples: int,
        ignore_nulls: bool,
    ) -> Self: ...
    def exp(self) -> Self: ...
    def sqrt(self) -> Self: ...
    def fill_nan(self, value: float | None) -> Self: ...
    def fill_null(
        self,
        value: Self | NonNestedLiteral,
        strategy: FillNullStrategy | None,
        limit: int | None,
    ) -> Self: ...
    def is_between(
        self, lower_bound: Self, upper_bound: Self, closed: ClosedInterval
    ) -> Self:
        if closed == "left":
            return (self >= lower_bound) & (self < upper_bound)
        if closed == "right":
            return (self > lower_bound) & (self <= upper_bound)
        if closed == "none":
            return (self > lower_bound) & (self < upper_bound)
        return (self >= lower_bound) & (self <= upper_bound)

    def is_close(
        self,
        other: Self | NumericLiteral,
        *,
        abs_tol: float,
        rel_tol: float,
        nans_equal: bool,
    ) -> Self:
        from decimal import Decimal

        other_abs: Self | NumericLiteral
        other_is_nan: Self | bool
        other_is_inf: Self | bool
        other_is_not_inf: Self | bool

        if isinstance(other, (float, int, Decimal)):
            from math import isinf, isnan

            # NOTE: See https://discuss.python.org/t/inferred-type-of-function-that-calls-dunder-abs-abs/101447
            other_abs = other.__abs__()
            other_is_nan = isnan(other)
            other_is_inf = isinf(other)

            # Define the other_is_not_inf variable to prevent triggering the following warning:
            # > DeprecationWarning: Bitwise inversion '~' on bool is deprecated and will be
            # >     removed in Python 3.16.
            other_is_not_inf = not other_is_inf

        else:
            other_abs, other_is_nan = other.abs(), other.is_nan()
            other_is_not_inf = other.is_finite() | other_is_nan
            other_is_inf = ~other_is_not_inf

        rel_threshold = self.abs().clip(lower_bound=other_abs, upper_bound=None) * rel_tol
        tolerance = rel_threshold.clip(lower_bound=abs_tol, upper_bound=None)

        self_is_nan = self.is_nan()
        self_is_not_inf = self.is_finite() | self_is_nan

        # Values are close if abs_diff <= tolerance, and both finite
        is_close = (
            ((self - other).abs() <= tolerance) & self_is_not_inf & other_is_not_inf
        )

        # Handle infinity cases: infinities are close/equal if they have the same sign
        self_sign, other_sign = self > 0, other > 0
        is_same_inf = (~self_is_not_inf) & other_is_inf & (self_sign == other_sign)

        # Handle nan cases:
        #   * If any value is NaN, then False (via `& ~either_nan`)
        #   * However, if `nans_equals = True` and if _both_ values are NaN, then True
        either_nan = self_is_nan | other_is_nan
        result = (is_close | is_same_inf) & ~either_nan

        if nans_equal:
            both_nan = self_is_nan & other_is_nan
            result = result | both_nan

        return result

    def is_duplicated(self) -> Self:
        return ~self.is_unique()

    def is_finite(self) -> Self: ...
    def is_first_distinct(self) -> Self: ...
    def is_in(self, other: Any) -> Self: ...
    def is_last_distinct(self) -> Self: ...
    def is_nan(self) -> Self: ...
    def is_null(self) -> Self: ...
    def is_unique(self) -> Self: ...
    def log(self, base: float) -> Self: ...
    def mode(self, *, keep: ModeKeepStrategy) -> Self: ...
    def rank(self, method: RankMethod, *, descending: bool) -> Self: ...
    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> Self: ...
    def rolling_mean(
        self, window_size: int, *, min_samples: int, center: bool
    ) -> Self: ...
    def rolling_std(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self: ...
    def rolling_sum(
        self, window_size: int, *, min_samples: int, center: bool
    ) -> Self: ...
    def rolling_var(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self: ...
    def round(self, decimals: int) -> Self: ...
    def floor(self) -> Self: ...
    def ceil(self) -> Self: ...
    def shift(self, n: int) -> Self: ...
    def unique(self) -> Self: ...

    @property
    def str(self) -> StringNamespace[Self]: ...
    @property
    def dt(self) -> DateTimeNamespace[Self]: ...
    @property
    def cat(self) -> CatNamespace[Self]: ...
    @property
    def list(self) -> ListNamespace[Self]: ...
    @property
    def struct(self) -> StructNamespace[Self]: ...
