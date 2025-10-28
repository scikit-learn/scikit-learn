from __future__ import annotations

from typing import TYPE_CHECKING, Any, Callable, ClassVar, Literal

import polars as pl

from narwhals._polars.utils import (
    BACKEND_VERSION,
    PolarsAnyNamespace,
    PolarsCatNamespace,
    PolarsDateTimeNamespace,
    PolarsListNamespace,
    PolarsStringNamespace,
    PolarsStructNamespace,
    extract_args_kwargs,
    extract_native,
    narwhals_to_native_dtype,
)
from narwhals._utils import Implementation, requires

if TYPE_CHECKING:
    from collections.abc import Mapping, Sequence

    from typing_extensions import Self

    from narwhals._compliant.typing import Accessor
    from narwhals._expression_parsing import ExprKind, ExprMetadata
    from narwhals._polars.dataframe import Method
    from narwhals._polars.namespace import PolarsNamespace
    from narwhals._utils import Version
    from narwhals.typing import IntoDType, ModeKeepStrategy, NumericLiteral


class PolarsExpr:
    # CompliantExpr
    _implementation: Implementation = Implementation.POLARS
    _version: Version
    _native_expr: pl.Expr
    _metadata: ExprMetadata | None = None
    _evaluate_output_names: Any
    _alias_output_names: Any
    __call__: Any

    # CompliantExpr + builtin descriptor
    # TODO @dangotbanned: Remove in #2713
    @classmethod
    def from_column_names(cls, *_: Any, **__: Any) -> Self:
        raise NotImplementedError

    @classmethod
    def from_column_indices(cls, *_: Any, **__: Any) -> Self:
        raise NotImplementedError

    def __narwhals_expr__(self) -> Self:  # pragma: no cover
        return self

    def __narwhals_namespace__(self) -> PolarsNamespace:  # pragma: no cover
        from narwhals._polars.namespace import PolarsNamespace

        return PolarsNamespace(version=self._version)

    def __init__(self, expr: pl.Expr, version: Version) -> None:
        self._native_expr = expr
        self._version = version

    @property
    def _backend_version(self) -> tuple[int, ...]:
        return self._implementation._backend_version()

    @property
    def native(self) -> pl.Expr:
        return self._native_expr

    def __repr__(self) -> str:  # pragma: no cover
        return "PolarsExpr"

    def _with_native(self, expr: pl.Expr) -> Self:
        return self.__class__(expr, self._version)

    def broadcast(self, kind: Literal[ExprKind.AGGREGATION, ExprKind.LITERAL]) -> Self:
        # Let Polars do its thing.
        return self

    def __getattr__(self, attr: str) -> Any:
        def func(*args: Any, **kwargs: Any) -> Any:
            pos, kwds = extract_args_kwargs(args, kwargs)
            return self._with_native(getattr(self.native, attr)(*pos, **kwds))

        return func

    def _renamed_min_periods(self, min_samples: int, /) -> dict[str, Any]:
        name = "min_periods" if self._backend_version < (1, 21, 0) else "min_samples"
        return {name: min_samples}

    def cast(self, dtype: IntoDType) -> Self:
        dtype_pl = narwhals_to_native_dtype(dtype, self._version)
        return self._with_native(self.native.cast(dtype_pl))

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
    ) -> Self:
        native = self.native.ewm_mean(
            com=com,
            span=span,
            half_life=half_life,
            alpha=alpha,
            adjust=adjust,
            ignore_nulls=ignore_nulls,
            **self._renamed_min_periods(min_samples),
        )
        if self._backend_version < (1,):  # pragma: no cover
            native = pl.when(~self.native.is_null()).then(native).otherwise(None)
        return self._with_native(native)

    def is_nan(self) -> Self:
        if self._backend_version >= (1, 18):
            native = self.native.is_nan()
        else:  # pragma: no cover
            native = pl.when(self.native.is_not_null()).then(self.native.is_nan())
        return self._with_native(native)

    def over(self, partition_by: Sequence[str], order_by: Sequence[str]) -> Self:
        # Use `pl.repeat(1, pl.len())` instead of `pl.lit(1)` to avoid issues for
        # non-numeric types: https://github.com/pola-rs/polars/issues/24756.
        pl_partition_by = partition_by or pl.repeat(1, pl.len())
        if self._backend_version < (1, 9):
            if order_by:
                msg = "`order_by` in Polars requires version 1.10 or greater"
                raise NotImplementedError(msg)
            native = self.native.over(pl_partition_by)
        else:
            native = self.native.over(pl_partition_by, order_by=order_by or None)
        return self._with_native(native)

    @requires.backend_version((1,))
    def rolling_var(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        kwds = self._renamed_min_periods(min_samples)
        native = self.native.rolling_var(
            window_size=window_size, center=center, ddof=ddof, **kwds
        )
        return self._with_native(native)

    @requires.backend_version((1,))
    def rolling_std(
        self, window_size: int, *, min_samples: int, center: bool, ddof: int
    ) -> Self:
        kwds = self._renamed_min_periods(min_samples)
        native = self.native.rolling_std(
            window_size=window_size, center=center, ddof=ddof, **kwds
        )
        return self._with_native(native)

    def rolling_sum(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        kwds = self._renamed_min_periods(min_samples)
        native = self.native.rolling_sum(window_size=window_size, center=center, **kwds)
        return self._with_native(native)

    def rolling_mean(self, window_size: int, *, min_samples: int, center: bool) -> Self:
        kwds = self._renamed_min_periods(min_samples)
        native = self.native.rolling_mean(window_size=window_size, center=center, **kwds)
        return self._with_native(native)

    def map_batches(
        self,
        function: Callable[[Any], Any],
        return_dtype: IntoDType | None,
        *,
        returns_scalar: bool,
    ) -> Self:
        pl_version = self._backend_version
        return_dtype_pl = (
            narwhals_to_native_dtype(return_dtype, self._version)
            if return_dtype is not None
            else None
            if pl_version < (1, 32)
            else pl.self_dtype()
        )
        kwargs = {} if pl_version < (0, 20, 31) else {"returns_scalar": returns_scalar}
        native = self.native.map_batches(function, return_dtype_pl, **kwargs)
        return self._with_native(native)

    @requires.backend_version((1,))
    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any],
        *,
        return_dtype: IntoDType | None,
    ) -> Self:
        return_dtype_pl = (
            narwhals_to_native_dtype(return_dtype, self._version)
            if return_dtype
            else None
        )
        native = self.native.replace_strict(old, new, return_dtype=return_dtype_pl)
        return self._with_native(native)

    def __eq__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_native(self.native.__eq__(extract_native(other)))  # type: ignore[operator]

    def __ne__(self, other: object) -> Self:  # type: ignore[override]
        return self._with_native(self.native.__ne__(extract_native(other)))  # type: ignore[operator]

    def __ge__(self, other: Any) -> Self:
        return self._with_native(self.native.__ge__(extract_native(other)))

    def __gt__(self, other: Any) -> Self:
        return self._with_native(self.native.__gt__(extract_native(other)))

    def __le__(self, other: Any) -> Self:
        return self._with_native(self.native.__le__(extract_native(other)))

    def __lt__(self, other: Any) -> Self:
        return self._with_native(self.native.__lt__(extract_native(other)))

    def __and__(self, other: PolarsExpr | bool | Any) -> Self:
        return self._with_native(self.native.__and__(extract_native(other)))  # type: ignore[operator]

    def __or__(self, other: PolarsExpr | bool | Any) -> Self:
        return self._with_native(self.native.__or__(extract_native(other)))  # type: ignore[operator]

    def __add__(self, other: Any) -> Self:
        return self._with_native(self.native.__add__(extract_native(other)))

    def __sub__(self, other: Any) -> Self:
        return self._with_native(self.native.__sub__(extract_native(other)))

    def __mul__(self, other: Any) -> Self:
        return self._with_native(self.native.__mul__(extract_native(other)))

    def __pow__(self, other: Any) -> Self:
        return self._with_native(self.native.__pow__(extract_native(other)))

    def __truediv__(self, other: Any) -> Self:
        return self._with_native(self.native.__truediv__(extract_native(other)))

    def __floordiv__(self, other: Any) -> Self:
        return self._with_native(self.native.__floordiv__(extract_native(other)))

    def __rfloordiv__(self, other: Any) -> Self:
        native = self.native
        result = native.__rfloordiv__(extract_native(other))
        if self._backend_version < (1, 10, 0):
            # Polars 1.9.0 and earlier returns 0 for division by 0 in rfloordiv.
            result = pl.when(native != 0).then(result).otherwise(None)
        return self._with_native(result)

    def __mod__(self, other: Any) -> Self:
        return self._with_native(self.native.__mod__(extract_native(other)))

    def __invert__(self) -> Self:
        return self._with_native(self.native.__invert__())

    def cum_count(self, *, reverse: bool) -> Self:
        return self._with_native(self.native.cum_count(reverse=reverse))

    def is_close(
        self,
        other: Self | NumericLiteral,
        *,
        abs_tol: float,
        rel_tol: float,
        nans_equal: bool,
    ) -> Self:
        left = self.native
        right = other.native if isinstance(other, PolarsExpr) else pl.lit(other)

        if self._backend_version < (1, 32, 0):
            lower_bound = right.abs()
            tolerance = (left.abs().clip(lower_bound) * rel_tol).clip(abs_tol)

            # Values are close if abs_diff <= tolerance, and both finite
            abs_diff = (left - right).abs()
            all_ = pl.all_horizontal
            is_close = all_((abs_diff <= tolerance), left.is_finite(), right.is_finite())

            # Handle infinity cases: infinities are "close" only if they have the same sign
            is_same_inf = all_(
                left.is_infinite(), right.is_infinite(), (left.sign() == right.sign())
            )

            # Handle nan cases:
            #   * nans_equals = True => if both values are NaN, then True
            #   * nans_equals = False => if any value is NaN, then False
            left_is_nan, right_is_nan = left.is_nan(), right.is_nan()
            either_nan = left_is_nan | right_is_nan
            result = (is_close | is_same_inf) & either_nan.not_()

            if nans_equal:
                result = result | (left_is_nan & right_is_nan)
        else:
            result = left.is_close(
                right, abs_tol=abs_tol, rel_tol=rel_tol, nans_equal=nans_equal
            )
        return self._with_native(result)

    def mode(self, *, keep: ModeKeepStrategy) -> Self:
        result = self.native.mode()
        return self._with_native(result.first() if keep == "any" else result)

    @property
    def dt(self) -> PolarsExprDateTimeNamespace:
        return PolarsExprDateTimeNamespace(self)

    @property
    def str(self) -> PolarsExprStringNamespace:
        return PolarsExprStringNamespace(self)

    @property
    def cat(self) -> PolarsExprCatNamespace:
        return PolarsExprCatNamespace(self)

    @property
    def name(self) -> PolarsExprNameNamespace:
        return PolarsExprNameNamespace(self)

    @property
    def list(self) -> PolarsExprListNamespace:
        return PolarsExprListNamespace(self)

    @property
    def struct(self) -> PolarsExprStructNamespace:
        return PolarsExprStructNamespace(self)

    # Polars
    abs: Method[Self]
    all: Method[Self]
    any: Method[Self]
    alias: Method[Self]
    arg_max: Method[Self]
    arg_min: Method[Self]
    arg_true: Method[Self]
    ceil: Method[Self]
    clip: Method[Self]
    count: Method[Self]
    cum_max: Method[Self]
    cum_min: Method[Self]
    cum_prod: Method[Self]
    cum_sum: Method[Self]
    diff: Method[Self]
    drop_nulls: Method[Self]
    exp: Method[Self]
    fill_null: Method[Self]
    fill_nan: Method[Self]
    first: Method[Self]
    floor: Method[Self]
    last: Method[Self]
    gather_every: Method[Self]
    head: Method[Self]
    is_between: Method[Self]
    is_duplicated: Method[Self]
    is_finite: Method[Self]
    is_first_distinct: Method[Self]
    is_in: Method[Self]
    is_last_distinct: Method[Self]
    is_null: Method[Self]
    is_unique: Method[Self]
    kurtosis: Method[Self]
    len: Method[Self]
    log: Method[Self]
    max: Method[Self]
    mean: Method[Self]
    median: Method[Self]
    min: Method[Self]
    n_unique: Method[Self]
    null_count: Method[Self]
    quantile: Method[Self]
    rank: Method[Self]
    round: Method[Self]
    sample: Method[Self]
    shift: Method[Self]
    skew: Method[Self]
    sqrt: Method[Self]
    std: Method[Self]
    sum: Method[Self]
    sort: Method[Self]
    tail: Method[Self]
    unique: Method[Self]
    var: Method[Self]
    __rsub__: Method[Self]
    __rmod__: Method[Self]
    __rpow__: Method[Self]
    __rtruediv__: Method[Self]


class PolarsExprNamespace(PolarsAnyNamespace[PolarsExpr, pl.Expr]):
    def __init__(self, expr: PolarsExpr) -> None:
        self._expr = expr

    @property
    def compliant(self) -> PolarsExpr:
        return self._expr

    @property
    def native(self) -> pl.Expr:
        return self._expr.native


class PolarsExprDateTimeNamespace(
    PolarsExprNamespace, PolarsDateTimeNamespace[PolarsExpr, pl.Expr]
): ...


class PolarsExprStringNamespace(
    PolarsExprNamespace, PolarsStringNamespace[PolarsExpr, pl.Expr]
):
    def to_titlecase(self) -> PolarsExpr:
        native_expr = self.native

        if BACKEND_VERSION < (1, 5):
            native_result = (
                native_expr.str.to_lowercase()
                .str.extract_all(r"[a-z0-9]*[^a-z0-9]*")
                .list.eval(pl.element().str.to_titlecase())
                .list.join("")
            )
        else:
            native_result = native_expr.str.to_titlecase()

        return self.compliant._with_native(native_result)

    @requires.backend_version((0, 20, 5))
    def zfill(self, width: int) -> PolarsExpr:
        backend_version = self.compliant._backend_version
        native_result = self.native.str.zfill(width)

        if backend_version <= (1, 30, 0):
            length = self.native.str.len_chars()
            less_than_width = length < width
            plus = "+"
            starts_with_plus = self.native.str.starts_with(plus)
            native_result = (
                pl.when(starts_with_plus & less_than_width)
                .then(
                    self.native.str.slice(1, length)
                    .str.zfill(width - 1)
                    .str.pad_start(width, plus)
                )
                .otherwise(native_result)
            )

        return self.compliant._with_native(native_result)


class PolarsExprCatNamespace(
    PolarsExprNamespace, PolarsCatNamespace[PolarsExpr, pl.Expr]
): ...


class PolarsExprNameNamespace(PolarsExprNamespace):
    _accessor: ClassVar[Accessor] = "name"
    keep: Method[PolarsExpr]
    map: Method[PolarsExpr]
    prefix: Method[PolarsExpr]
    suffix: Method[PolarsExpr]
    to_lowercase: Method[PolarsExpr]
    to_uppercase: Method[PolarsExpr]


class PolarsExprListNamespace(
    PolarsExprNamespace, PolarsListNamespace[PolarsExpr, pl.Expr]
):
    def len(self) -> PolarsExpr:
        native_expr = self.native
        native_result = native_expr.list.len()

        if self.compliant._backend_version < (1, 16):  # pragma: no cover
            native_result = (
                pl.when(~native_expr.is_null()).then(native_result).cast(pl.UInt32())
            )
        elif self.compliant._backend_version < (1, 17):  # pragma: no cover
            native_result = native_result.cast(pl.UInt32())

        return self.compliant._with_native(native_result)

    def contains(self, item: Any) -> PolarsExpr:
        if self.compliant._backend_version < (1, 28):
            result: pl.Expr = pl.when(self.native.is_not_null()).then(
                self.native.list.contains(item)
            )
        else:
            result = self.native.list.contains(item)
        return self.compliant._with_native(result)


class PolarsExprStructNamespace(
    PolarsExprNamespace, PolarsStructNamespace[PolarsExpr, pl.Expr]
): ...
