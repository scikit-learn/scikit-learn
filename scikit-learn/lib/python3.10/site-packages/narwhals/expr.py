from __future__ import annotations

import math
import operator as op
from collections.abc import Iterable, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Callable

from narwhals._expression_parsing import (
    ExprMetadata,
    apply_n_ary_operation,
    combine_metadata,
)
from narwhals._utils import _validate_rolling_arguments, ensure_type, flatten
from narwhals.dtypes import _validate_dtype
from narwhals.exceptions import ComputeError, InvalidOperationError
from narwhals.expr_cat import ExprCatNamespace
from narwhals.expr_dt import ExprDateTimeNamespace
from narwhals.expr_list import ExprListNamespace
from narwhals.expr_name import ExprNameNamespace
from narwhals.expr_str import ExprStringNamespace
from narwhals.expr_struct import ExprStructNamespace
from narwhals.translate import to_native

if TYPE_CHECKING:
    from typing import NoReturn, TypeVar

    from typing_extensions import Concatenate, ParamSpec, Self, TypeAlias

    from narwhals._compliant import CompliantExpr, CompliantNamespace
    from narwhals.dtypes import DType
    from narwhals.typing import (
        ClosedInterval,
        FillNullStrategy,
        IntoDType,
        IntoExpr,
        ModeKeepStrategy,
        NonNestedLiteral,
        NumericLiteral,
        RankMethod,
        RollingInterpolationMethod,
        TemporalLiteral,
        _1DArray,
    )

    PS = ParamSpec("PS")
    R = TypeVar("R")
    _ToCompliant: TypeAlias = Callable[
        [CompliantNamespace[Any, Any]], CompliantExpr[Any, Any]
    ]


class Expr:
    def __init__(self, to_compliant_expr: _ToCompliant, metadata: ExprMetadata) -> None:
        # callable from CompliantNamespace to CompliantExpr
        def func(plx: CompliantNamespace[Any, Any]) -> CompliantExpr[Any, Any]:
            result = to_compliant_expr(plx)
            result._metadata = self._metadata
            return result

        self._to_compliant_expr: _ToCompliant = func
        self._metadata = metadata

    def _with_elementwise(self, to_compliant_expr: Callable[[Any], Any]) -> Self:
        return self.__class__(to_compliant_expr, self._metadata.with_elementwise_op())

    def _with_aggregation(self, to_compliant_expr: Callable[[Any], Any]) -> Self:
        return self.__class__(to_compliant_expr, self._metadata.with_aggregation())

    def _with_orderable_aggregation(
        self, to_compliant_expr: Callable[[Any], Any]
    ) -> Self:
        return self.__class__(
            to_compliant_expr, self._metadata.with_orderable_aggregation()
        )

    def _with_orderable_window(self, to_compliant_expr: Callable[[Any], Any]) -> Self:
        return self.__class__(to_compliant_expr, self._metadata.with_orderable_window())

    def _with_window(self, to_compliant_expr: Callable[[Any], Any]) -> Self:
        return self.__class__(to_compliant_expr, self._metadata.with_window())

    def _with_filtration(self, to_compliant_expr: Callable[[Any], Any]) -> Self:
        return self.__class__(to_compliant_expr, self._metadata.with_filtration())

    def _with_orderable_filtration(self, to_compliant_expr: Callable[[Any], Any]) -> Self:
        return self.__class__(
            to_compliant_expr, self._metadata.with_orderable_filtration()
        )

    def _with_nary(
        self,
        n_ary_function: Callable[..., Any],
        *args: IntoExpr | NonNestedLiteral | _1DArray,
    ) -> Self:
        return self.__class__(
            lambda plx: apply_n_ary_operation(
                plx, n_ary_function, self, *args, str_as_lit=False
            ),
            combine_metadata(
                self,
                *args,
                str_as_lit=False,
                allow_multi_output=False,
                to_single_output=False,
            ),
        )

    def __repr__(self) -> str:
        return f"Narwhals Expr\nmetadata: {self._metadata}\n"

    def __bool__(self) -> NoReturn:
        msg = (
            f"the truth value of {type(self)} is ambiguous"
            "\n\n"
            "You probably got here by using a Python standard library function instead "
            "of the native expressions API.\n"
            "Here are some things you might want to try:\n"
            "- instead of `nw.col('a') and nw.col('b')`, use `nw.col('a') & nw.col('b')`\n"
            "- instead of `nw.col('a') in [y, z]`, use `nw.col('a').is_in([y, z])`\n"
            "- instead of `max(nw.col('a'), nw.col('b'))`, use `nw.max_horizontal(nw.col('a'), nw.col('b'))`\n"
        )
        raise TypeError(msg)

    def _taxicab_norm(self) -> Self:
        # This is just used to test out the stable api feature in a realistic-ish way.
        # It's not intended to be used.
        return self._with_aggregation(
            lambda plx: self._to_compliant_expr(plx).abs().sum()
        )

    # --- convert ---
    def alias(self, name: str) -> Self:
        """Rename the expression.

        Arguments:
            name: The new name.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2], "b": [4, 5]})
            >>> df = nw.from_native(df_native)
            >>> df.select((nw.col("b") + 10).alias("c"))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |          c       |
            |      0  14       |
            |      1  15       |
            └──────────────────┘
        """
        # Don't use `_with_elementwise` so that `_metadata.last_node` is preserved.
        return self.__class__(
            lambda plx: self._to_compliant_expr(plx).alias(name), self._metadata
        )

    def pipe(
        self,
        function: Callable[Concatenate[Self, PS], R],
        *args: PS.args,
        **kwargs: PS.kwargs,
    ) -> R:
        """Pipe function call.

        Arguments:
            function: Function to apply.
            args: Positional arguments to pass to function.
            kwargs: Keyword arguments to pass to function.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 4]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_piped=nw.col("a").pipe(lambda x: x + 1))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |     a  a_piped   |
            |  0  1        2   |
            |  1  2        3   |
            |  2  3        4   |
            |  3  4        5   |
            └──────────────────┘
        """
        return function(self, *args, **kwargs)

    def cast(self, dtype: IntoDType) -> Self:
        """Redefine an object's data type.

        Arguments:
            dtype: Data type that the object will be cast into.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"foo": [1, 2, 3], "bar": [6.0, 7.0, 8.0]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("foo").cast(nw.Float32), nw.col("bar").cast(nw.UInt8))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |      foo  bar    |
            |   0  1.0    6    |
            |   1  2.0    7    |
            |   2  3.0    8    |
            └──────────────────┘
        """
        _validate_dtype(dtype)
        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).cast(dtype)
        )

    # --- binary ---
    def _with_binary(
        self,
        function: Callable[[Any, Any], Any],
        other: Self | Any,
        *,
        str_as_lit: bool = True,
    ) -> Self:
        return self.__class__(
            lambda plx: apply_n_ary_operation(
                plx, function, self, other, str_as_lit=str_as_lit
            ),
            ExprMetadata.from_binary_op(self, other),
        )

    def __eq__(self, other: Self | Any) -> Self:  # type: ignore[override]
        return self._with_binary(op.eq, other)

    def __ne__(self, other: Self | Any) -> Self:  # type: ignore[override]
        return self._with_binary(op.ne, other)

    def __and__(self, other: Any) -> Self:
        return self._with_binary(op.and_, other)

    def __rand__(self, other: Any) -> Self:
        return (self & other).alias("literal")  # type: ignore[no-any-return]

    def __or__(self, other: Any) -> Self:
        return self._with_binary(op.or_, other)

    def __ror__(self, other: Any) -> Self:
        return (self | other).alias("literal")  # type: ignore[no-any-return]

    def __add__(self, other: Any) -> Self:
        return self._with_binary(op.add, other)

    def __radd__(self, other: Any) -> Self:
        return (self + other).alias("literal")  # type: ignore[no-any-return]

    def __sub__(self, other: Any) -> Self:
        return self._with_binary(op.sub, other)

    def __rsub__(self, other: Any) -> Self:
        return self._with_binary(lambda x, y: x.__rsub__(y), other)

    def __truediv__(self, other: Any) -> Self:
        return self._with_binary(op.truediv, other)

    def __rtruediv__(self, other: Any) -> Self:
        return self._with_binary(lambda x, y: x.__rtruediv__(y), other)

    def __mul__(self, other: Any) -> Self:
        return self._with_binary(op.mul, other)

    def __rmul__(self, other: Any) -> Self:
        return (self * other).alias("literal")  # type: ignore[no-any-return]

    def __le__(self, other: Any) -> Self:
        return self._with_binary(op.le, other)

    def __lt__(self, other: Any) -> Self:
        return self._with_binary(op.lt, other)

    def __gt__(self, other: Any) -> Self:
        return self._with_binary(op.gt, other)

    def __ge__(self, other: Any) -> Self:
        return self._with_binary(op.ge, other)

    def __pow__(self, other: Any) -> Self:
        return self._with_binary(op.pow, other)

    def __rpow__(self, other: Any) -> Self:
        return self._with_binary(lambda x, y: x.__rpow__(y), other)

    def __floordiv__(self, other: Any) -> Self:
        return self._with_binary(op.floordiv, other)

    def __rfloordiv__(self, other: Any) -> Self:
        return self._with_binary(lambda x, y: x.__rfloordiv__(y), other)

    def __mod__(self, other: Any) -> Self:
        return self._with_binary(op.mod, other)

    def __rmod__(self, other: Any) -> Self:
        return self._with_binary(lambda x, y: x.__rmod__(y), other)

    # --- unary ---
    def __invert__(self) -> Self:
        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).__invert__()
        )

    def any(self) -> Self:
        """Return whether any of the values in the column are `True`.

        If there are no non-null elements, the result is `False`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [True, False], "b": [True, True]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").any())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a     b   |
            |  0  True  True   |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).any())

    def all(self) -> Self:
        """Return whether all values in the column are `True`.

        If there are no non-null elements, the result is `True`.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [True, False], "b": [True, True]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").all())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |         a     b  |
            |  0  False  True  |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).all())

    def ewm_mean(
        self,
        *,
        com: float | None = None,
        span: float | None = None,
        half_life: float | None = None,
        alpha: float | None = None,
        adjust: bool = True,
        min_samples: int = 1,
        ignore_nulls: bool = False,
    ) -> Self:
        r"""Compute exponentially-weighted moving average.

        Arguments:
            com: Specify decay in terms of center of mass, $\gamma$, with <br> $\alpha = \frac{1}{1+\gamma}\forall\gamma\geq0$
            span: Specify decay in terms of span, $\theta$, with <br> $\alpha = \frac{2}{\theta + 1} \forall \theta \geq 1$
            half_life: Specify decay in terms of half-life, $\tau$, with <br> $\alpha = 1 - \exp \left\{ \frac{ -\ln(2) }{ \tau } \right\} \forall \tau > 0$
            alpha: Specify smoothing factor alpha directly, $0 < \alpha \leq 1$.
            adjust: Divide by decaying adjustment factor in beginning periods to account for imbalance in relative weightings

                - When `adjust=True` (the default) the EW function is calculated
                  using weights $w_i = (1 - \alpha)^i$
                - When `adjust=False` the EW function is calculated recursively by
                  $$
                  y_0=x_0
                  $$
                  $$
                  y_t = (1 - \alpha)y_{t - 1} + \alpha x_t
                  $$
            min_samples: Minimum number of observations in window required to have a value, (otherwise result is null).
            ignore_nulls: Ignore missing values when calculating weights.

                - When `ignore_nulls=False` (default), weights are based on absolute
                  positions.
                  For example, the weights of $x_0$ and $x_2$ used in
                  calculating the final weighted average of $[x_0, None, x_2]$ are
                  $(1-\alpha)^2$ and $1$ if `adjust=True`, and
                  $(1-\alpha)^2$ and $\alpha$ if `adjust=False`.
                - When `ignore_nulls=True`, weights are based
                  on relative positions. For example, the weights of
                  $x_0$ and $x_2$ used in calculating the final weighted
                  average of $[x_0, None, x_2]$ are
                  $1-\alpha$ and $1$ if `adjust=True`,
                  and $1-\alpha$ and $\alpha$ if `adjust=False`.

        Examples:
            >>> import pandas as pd
            >>> import polars as pl
            >>> import narwhals as nw
            >>> from narwhals.typing import IntoFrameT
            >>>
            >>> data = {"a": [1, 2, 3]}
            >>> df_pd = pd.DataFrame(data)
            >>> df_pl = pl.DataFrame(data)

            We define a library agnostic function:

            >>> def agnostic_ewm_mean(df_native: IntoFrameT) -> IntoFrameT:
            ...     df = nw.from_native(df_native)
            ...     return df.select(
            ...         nw.col("a").ewm_mean(com=1, ignore_nulls=False)
            ...     ).to_native()

            We can then pass either pandas or Polars to `agnostic_ewm_mean`:

            >>> agnostic_ewm_mean(df_pd)
                      a
            0  1.000000
            1  1.666667
            2  2.428571

            >>> agnostic_ewm_mean(df_pl)  # doctest: +NORMALIZE_WHITESPACE
            shape: (3, 1)
            ┌──────────┐
            │ a        │
            │ ---      │
            │ f64      │
            ╞══════════╡
            │ 1.0      │
            │ 1.666667 │
            │ 2.428571 │
            └──────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).ewm_mean(
                com=com,
                span=span,
                half_life=half_life,
                alpha=alpha,
                adjust=adjust,
                min_samples=min_samples,
                ignore_nulls=ignore_nulls,
            )
        )

    def mean(self) -> Self:
        """Get mean value.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [-1, 0, 1], "b": [2, 4, 6]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").mean())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a    b    |
            |   0  0.0  4.0    |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).mean())

    def median(self) -> Self:
        """Get median value.

        Notes:
            Results might slightly differ across backends due to differences in the underlying algorithms used to compute the median.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 8, 3], "b": [4, 5, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").median())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a    b    |
            |   0  3.0  4.0    |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).median())

    def std(self, *, ddof: int = 1) -> Self:
        """Get standard deviation.

        Arguments:
            ddof: "Delta Degrees of Freedom": the divisor used in the calculation is N - ddof,
                where N represents the number of elements. By default ddof is 1.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [20, 25, 60], "b": [1.5, 1, -1.4]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").std(ddof=0))
            ┌─────────────────────┐
            | Narwhals DataFrame  |
            |---------------------|
            |          a         b|
            |0  17.79513  1.265789|
            └─────────────────────┘
        """
        return self._with_aggregation(
            lambda plx: self._to_compliant_expr(plx).std(ddof=ddof)
        )

    def var(self, *, ddof: int = 1) -> Self:
        """Get variance.

        Arguments:
            ddof: "Delta Degrees of Freedom": the divisor used in the calculation is N - ddof,
                     where N represents the number of elements. By default ddof is 1.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [20, 25, 60], "b": [1.5, 1, -1.4]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").var(ddof=0))
            ┌───────────────────────┐
            |  Narwhals DataFrame   |
            |-----------------------|
            |            a         b|
            |0  316.666667  1.602222|
            └───────────────────────┘
        """
        return self._with_aggregation(
            lambda plx: self._to_compliant_expr(plx).var(ddof=ddof)
        )

    def map_batches(
        self,
        function: Callable[[Any], CompliantExpr[Any, Any]],
        return_dtype: DType | None = None,
        *,
        returns_scalar: bool = False,
    ) -> Self:
        """Apply a custom python function to a whole Series or sequence of Series.

        The output of this custom function is presumed to be either a Series,
        or a NumPy array (in which case it will be automatically converted into
        a Series).

        Arguments:
            function: Function to apply to Series.
            return_dtype: Dtype of the output Series.
                If not set, the dtype will be inferred based on the first non-null value
                that is returned by the function.
            returns_scalar: If the function returns a scalar, by default it will be wrapped
                in a list in the output, since the assumption is that the function always
                returns something Series-like.
                If you want to keep the result as a scalar, set this argument to True.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a", "b")
            ...     .map_batches(lambda s: s.to_numpy() + 1, return_dtype=nw.Float64)
            ...     .name.suffix("_mapped")
            ... )
            ┌───────────────────────────┐
            |    Narwhals DataFrame     |
            |---------------------------|
            |   a  b  a_mapped  b_mapped|
            |0  1  4       2.0       5.0|
            |1  2  5       3.0       6.0|
            |2  3  6       4.0       7.0|
            └───────────────────────────┘
        """

        def compliant_expr(plx: Any) -> Any:
            return self._to_compliant_expr(plx).map_batches(
                function=function,
                return_dtype=return_dtype,
                returns_scalar=returns_scalar,
            )

        if returns_scalar:
            return self._with_orderable_aggregation(compliant_expr)
        # safest assumptions
        return self._with_orderable_filtration(compliant_expr)

    def skew(self) -> Self:
        """Calculate the sample skewness of a column.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 4, 5], "b": [1, 1, 2, 10, 100]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").skew())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |      a         b |
            | 0  0.0  1.472427 |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).skew())

    def kurtosis(self) -> Self:
        """Compute the kurtosis (Fisher's definition) without bias correction.

        Kurtosis is the fourth central moment divided by the square of the variance.
        The Fisher's definition is used where 3.0 is subtracted from the result to give 0.0 for a normal distribution.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 4, 5], "b": [1, 1, 2, 10, 100]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").kurtosis())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |      a         b |
            | 0 -1.3  0.210657 |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).kurtosis())

    def sum(self) -> Expr:
        """Return the sum value.

        If there are no non-null elements, the result is zero.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql("SELECT * FROM VALUES (5, 50), (10, 100) df(a, b)")
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").sum())
            ┌───────────────────┐
            |Narwhals LazyFrame |
            |-------------------|
            |┌────────┬────────┐|
            |│   a    │   b    │|
            |│ int128 │ int128 │|
            |├────────┼────────┤|
            |│     15 │    150 │|
            |└────────┴────────┘|
            └───────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).sum())

    def min(self) -> Self:
        """Returns the minimum value(s) from a column(s).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2], "b": [4, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.min("a", "b"))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a  b      |
            |     0  1  3      |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).min())

    def max(self) -> Self:
        """Returns the maximum value(s) from a column(s).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [10, 20], "b": [50, 100]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.max("a", "b"))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a    b    |
            |    0  20  100    |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).max())

    def count(self) -> Self:
        """Returns the number of non-null elements in the column.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3], "b": [None, 4, 4]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.all().count())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a  b      |
            |     0  3  2      |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).count())

    def n_unique(self) -> Self:
        """Returns count of unique values.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 4, 5], "b": [1, 1, 3, 3, 5]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").n_unique())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a  b      |
            |     0  5  3      |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).n_unique())

    def unique(self) -> Self:
        """Return unique values of this expression.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 1, 3, 5, 5], "b": [2, 4, 4, 6, 6]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").unique().sum())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a   b     |
            |     0  9  12     |
            └──────────────────┘
        """
        return self._with_filtration(lambda plx: self._to_compliant_expr(plx).unique())

    def abs(self) -> Self:
        """Return absolute value of each element.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, -2], "b": [-3, 4]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a", "b").abs().name.suffix("_abs"))
            ┌─────────────────────┐
            | Narwhals DataFrame  |
            |---------------------|
            |   a  b  a_abs  b_abs|
            |0  1 -3      1      3|
            |1 -2  4      2      4|
            └─────────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).abs())

    def cum_sum(self, *, reverse: bool = False) -> Self:
        """Return cumulative sum.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 1, 3, 5, 5], "b": [2, 4, 4, 6, 6]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_cum_sum=nw.col("a").cum_sum())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |   a  b  a_cum_sum|
            |0  1  2          1|
            |1  1  4          2|
            |2  3  4          5|
            |3  5  6         10|
            |4  5  6         15|
            └──────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).cum_sum(reverse=reverse)
        )

    def diff(self) -> Self:
        """Returns the difference between each element and the previous one.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Notes:
            pandas may change the dtype here, for example when introducing missing
            values in an integer column. To ensure, that the dtype doesn't change,
            you may want to use `fill_null` and `cast`. For example, to calculate
            the diff and fill missing values with `0` in a Int64 column, you could
            do:

                nw.col("a").diff().fill_null(0).cast(nw.Int64)

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 1, 3, 5, 5]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_diff=nw.col("a").diff())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            | shape: (5, 2)    |
            | ┌─────┬────────┐ |
            | │ a   ┆ a_diff │ |
            | │ --- ┆ ---    │ |
            | │ i64 ┆ i64    │ |
            | ╞═════╪════════╡ |
            | │ 1   ┆ null   │ |
            | │ 1   ┆ 0      │ |
            | │ 3   ┆ 2      │ |
            | │ 5   ┆ 2      │ |
            | │ 5   ┆ 0      │ |
            | └─────┴────────┘ |
            └──────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).diff()
        )

    def shift(self, n: int) -> Self:
        """Shift values by `n` positions.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            n: Number of positions to shift values by.

        Notes:
            pandas may change the dtype here, for example when introducing missing
            values in an integer column. To ensure, that the dtype doesn't change,
            you may want to use `fill_null` and `cast`. For example, to shift
            and fill missing values with `0` in a Int64 column, you could
            do:

                nw.col("a").shift(1).fill_null(0).cast(nw.Int64)

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [1, 1, 3, 5, 5]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_shift=nw.col("a").shift(n=1))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |shape: (5, 2)     |
            |┌─────┬─────────┐ |
            |│ a   ┆ a_shift │ |
            |│ --- ┆ ---     │ |
            |│ i64 ┆ i64     │ |
            |╞═════╪═════════╡ |
            |│ 1   ┆ null    │ |
            |│ 1   ┆ 1       │ |
            |│ 3   ┆ 1       │ |
            |│ 5   ┆ 3       │ |
            |│ 5   ┆ 5       │ |
            |└─────┴─────────┘ |
            └──────────────────┘
        """
        ensure_type(n, int, param_name="n")

        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).shift(n)
        )

    def replace_strict(
        self,
        old: Sequence[Any] | Mapping[Any, Any],
        new: Sequence[Any] | None = None,
        *,
        return_dtype: IntoDType | None = None,
    ) -> Self:
        """Replace all values by different values.

        This function must replace all non-null input values (else it raises an error).

        Arguments:
            old: Sequence of values to replace. It also accepts a mapping of values to
                their replacement as syntactic sugar for
                `replace_strict(old=list(mapping.keys()), new=list(mapping.values()))`.
            new: Sequence of values to replace by. Length must match the length of `old`.
            return_dtype: The data type of the resulting expression. If set to `None`
                (default), the data type is determined automatically based on the other
                inputs.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [3, 0, 1, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     b=nw.col("a").replace_strict(
            ...         [0, 1, 2, 3],
            ...         ["zero", "one", "two", "three"],
            ...         return_dtype=nw.String,
            ...     )
            ... )
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |      a      b    |
            |   0  3  three    |
            |   1  0   zero    |
            |   2  1    one    |
            |   3  2    two    |
            └──────────────────┘
        """
        if new is None:
            if not isinstance(old, Mapping):
                msg = "`new` argument is required if `old` argument is not a Mapping type"
                raise TypeError(msg)

            new = list(old.values())
            old = list(old.keys())

        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).replace_strict(
                old, new, return_dtype=return_dtype
            )
        )

    # --- transform ---
    def is_between(
        self,
        lower_bound: Any | IntoExpr,
        upper_bound: Any | IntoExpr,
        closed: ClosedInterval = "both",
    ) -> Self:
        """Check if this expression is between the given lower and upper bounds.

        Arguments:
            lower_bound: Lower bound value. String literals are interpreted as column names.
            upper_bound: Upper bound value. String literals are interpreted as column names.
            closed: Define which sides of the interval are closed (inclusive).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 4, 5]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(b=nw.col("a").is_between(2, 4, "right"))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |      a      b    |
            |   0  1  False    |
            |   1  2  False    |
            |   2  3   True    |
            |   3  4   True    |
            |   4  5  False    |
            └──────────────────┘
        """
        return self._with_nary(
            lambda expr, lb, ub: expr.is_between(lb, ub, closed=closed),
            lower_bound,
            upper_bound,
        )

    def is_in(self, other: Any) -> Self:
        """Check if elements of this expression are present in the other iterable.

        Arguments:
            other: iterable

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 9, 10]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(b=nw.col("a").is_in([1, 2]))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |       a      b   |
            |   0   1   True   |
            |   1   2   True   |
            |   2   9  False   |
            |   3  10  False   |
            └──────────────────┘
        """
        if isinstance(other, Iterable) and not isinstance(other, (str, bytes)):
            return self._with_elementwise(
                lambda plx: self._to_compliant_expr(plx).is_in(
                    to_native(other, pass_through=True)
                )
            )
        msg = "Narwhals `is_in` doesn't accept expressions as an argument, as opposed to Polars. You should provide an iterable instead."
        raise NotImplementedError(msg)

    def filter(self, *predicates: Any) -> Self:
        """Filters elements based on a condition, returning a new expression.

        Arguments:
            predicates: Conditions to filter by (which get AND-ed together).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": [2, 3, 4, 5, 6, 7], "b": [10, 11, 12, 13, 14, 15]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(
            ...     nw.col("a").filter(nw.col("a") > 4),
            ...     nw.col("b").filter(nw.col("b") < 13),
            ... )
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a   b     |
            |     3  5  10     |
            |     4  6  11     |
            |     5  7  12     |
            └──────────────────┘
        """
        flat_predicates = flatten(predicates)
        metadata = combine_metadata(
            self,
            *flat_predicates,
            str_as_lit=False,
            allow_multi_output=True,
            to_single_output=False,
        ).with_filtration()
        return self.__class__(
            lambda plx: apply_n_ary_operation(
                plx,
                lambda *exprs: exprs[0].filter(*exprs[1:]),
                self,
                *flat_predicates,
                str_as_lit=False,
            ),
            metadata,
        )

    def is_null(self) -> Self:
        """Returns a boolean Series indicating which values are null.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES (null, CAST('NaN' AS DOUBLE)), (2, 2.) df(a, b)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_is_null=nw.col("a").is_null(), b_is_null=nw.col("b").is_null()
            ... )
            ┌──────────────────────────────────────────┐
            |            Narwhals LazyFrame            |
            |------------------------------------------|
            |┌───────┬────────┬───────────┬───────────┐|
            |│   a   │   b    │ a_is_null │ b_is_null │|
            |│ int32 │ double │  boolean  │  boolean  │|
            |├───────┼────────┼───────────┼───────────┤|
            |│  NULL │    nan │ true      │ false     │|
            |│     2 │    2.0 │ false     │ false     │|
            |└───────┴────────┴───────────┴───────────┘|
            └──────────────────────────────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).is_null())

    def is_nan(self) -> Self:
        """Indicate which values are NaN.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES (null, CAST('NaN' AS DOUBLE)), (2, 2.) df(a, b)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_is_nan=nw.col("a").is_nan(), b_is_nan=nw.col("b").is_nan()
            ... )
            ┌────────────────────────────────────────┐
            |           Narwhals LazyFrame           |
            |----------------------------------------|
            |┌───────┬────────┬──────────┬──────────┐|
            |│   a   │   b    │ a_is_nan │ b_is_nan │|
            |│ int32 │ double │ boolean  │ boolean  │|
            |├───────┼────────┼──────────┼──────────┤|
            |│  NULL │    nan │ NULL     │ true     │|
            |│     2 │    2.0 │ false    │ false    │|
            |└───────┴────────┴──────────┴──────────┘|
            └────────────────────────────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).is_nan())

    def fill_null(
        self,
        value: Expr | NonNestedLiteral = None,
        strategy: FillNullStrategy | None = None,
        limit: int | None = None,
    ) -> Self:
        """Fill null values with given value.

        Arguments:
            value: Value or expression used to fill null values.
            strategy: Strategy used to fill null values.
            limit: Number of consecutive null values to fill when using the 'forward' or 'backward' strategy.

        Notes:
            - pandas handles null values differently from other libraries.
              See [null_handling](../concepts/null_handling.md/)
              for reference.
            - For pandas Series of `object` dtype, `fill_null` will not automatically change the
              Series' dtype as pandas used to do. Explicitly call `cast` if you want the dtype to change.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame(
            ...     {
            ...         "a": [2, None, None, 3],
            ...         "b": [2.0, float("nan"), float("nan"), 3.0],
            ...         "c": [1, 2, 3, 4],
            ...     }
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a", "b").fill_null(0).name.suffix("_filled"),
            ...     nw.col("a").fill_null(nw.col("c")).name.suffix("_filled_with_c"),
            ... )
            ┌────────────────────────────────────────────────────────────┐
            |                     Narwhals DataFrame                     |
            |------------------------------------------------------------|
            |shape: (4, 6)                                               |
            |┌──────┬─────┬─────┬──────────┬──────────┬─────────────────┐|
            |│ a    ┆ b   ┆ c   ┆ a_filled ┆ b_filled ┆ a_filled_with_c │|
            |│ ---  ┆ --- ┆ --- ┆ ---      ┆ ---      ┆ ---             │|
            |│ i64  ┆ f64 ┆ i64 ┆ i64      ┆ f64      ┆ i64             │|
            |╞══════╪═════╪═════╪══════════╪══════════╪═════════════════╡|
            |│ 2    ┆ 2.0 ┆ 1   ┆ 2        ┆ 2.0      ┆ 2               │|
            |│ null ┆ NaN ┆ 2   ┆ 0        ┆ NaN      ┆ 2               │|
            |│ null ┆ NaN ┆ 3   ┆ 0        ┆ NaN      ┆ 3               │|
            |│ 3    ┆ 3.0 ┆ 4   ┆ 3        ┆ 3.0      ┆ 3               │|
            |└──────┴─────┴─────┴──────────┴──────────┴─────────────────┘|
            └────────────────────────────────────────────────────────────┘

            Using a strategy:

            >>> df.select(
            ...     nw.col("a", "b"),
            ...     nw.col("a", "b")
            ...     .fill_null(strategy="forward", limit=1)
            ...     .name.suffix("_nulls_forward_filled"),
            ... )
            ┌────────────────────────────────────────────────────────────────┐
            |                       Narwhals DataFrame                       |
            |----------------------------------------------------------------|
            |shape: (4, 4)                                                   |
            |┌──────┬─────┬────────────────────────┬────────────────────────┐|
            |│ a    ┆ b   ┆ a_nulls_forward_filled ┆ b_nulls_forward_filled │|
            |│ ---  ┆ --- ┆ ---                    ┆ ---                    │|
            |│ i64  ┆ f64 ┆ i64                    ┆ f64                    │|
            |╞══════╪═════╪════════════════════════╪════════════════════════╡|
            |│ 2    ┆ 2.0 ┆ 2                      ┆ 2.0                    │|
            |│ null ┆ NaN ┆ 2                      ┆ NaN                    │|
            |│ null ┆ NaN ┆ null                   ┆ NaN                    │|
            |│ 3    ┆ 3.0 ┆ 3                      ┆ 3.0                    │|
            |└──────┴─────┴────────────────────────┴────────────────────────┘|
            └────────────────────────────────────────────────────────────────┘
        """
        if value is not None and strategy is not None:
            msg = "cannot specify both `value` and `strategy`"
            raise ValueError(msg)
        if value is None and strategy is None:
            msg = "must specify either a fill `value` or `strategy`"
            raise ValueError(msg)
        if strategy is not None and strategy not in {"forward", "backward"}:
            msg = f"strategy not supported: {strategy}"
            raise ValueError(msg)

        return self.__class__(
            lambda plx: apply_n_ary_operation(
                plx,
                lambda *exprs: exprs[0].fill_null(
                    exprs[1], strategy=strategy, limit=limit
                ),
                self,
                value,
                str_as_lit=True,
            ),
            self._metadata.with_orderable_window()
            if strategy is not None
            else self._metadata,
        )

    def fill_nan(self, value: float | None) -> Self:
        """Fill floating point NaN values with given value.

        Arguments:
            value: Value used to fill NaN values.

        Notes:
            This function only fills `'NaN'` values, not null ones, except for pandas
            which doesn't distinguish between them.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import duckdb
            >>> import narwhals as nw
            >>> df_native = duckdb.sql(
            ...     "SELECT * FROM VALUES (5.::DOUBLE, 50.::DOUBLE), ('NaN', null) df(a, b)"
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.col("a", "b").fill_nan(0).name.suffix("_nans_filled"))
            ┌───────────────────────────────────────────────────┐
            |                Narwhals LazyFrame                 |
            |---------------------------------------------------|
            |┌────────┬────────┬───────────────┬───────────────┐|
            |│   a    │   b    │ a_nans_filled │ b_nans_filled │|
            |│ double │ double │    double     │    double     │|
            |├────────┼────────┼───────────────┼───────────────┤|
            |│    5.0 │   50.0 │           5.0 │          50.0 │|
            |│    nan │   NULL │           0.0 │          NULL │|
            |└────────┴────────┴───────────────┴───────────────┘|
            └───────────────────────────────────────────────────┘
        """
        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).fill_nan(value)
        )

    # --- partial reduction ---
    def drop_nulls(self) -> Self:
        """Drop null values.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [2.0, 4.0, float("nan"), 3.0, None, 5.0]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a").drop_nulls())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  shape: (5, 1)   |
            |  ┌─────┐         |
            |  │ a   │         |
            |  │ --- │         |
            |  │ f64 │         |
            |  ╞═════╡         |
            |  │ 2.0 │         |
            |  │ 4.0 │         |
            |  │ NaN │         |
            |  │ 3.0 │         |
            |  │ 5.0 │         |
            |  └─────┘         |
            └──────────────────┘
        """
        return self._with_filtration(
            lambda plx: self._to_compliant_expr(plx).drop_nulls()
        )

    def over(
        self,
        *partition_by: str | Sequence[str],
        order_by: str | Sequence[str] | None = None,
    ) -> Self:
        """Compute expressions over the given groups (optionally with given order).

        Arguments:
            partition_by: Names of columns to compute window expression over.
                Must be names of columns, as opposed to expressions -
                so, this is a bit less flexible than Polars' `Expr.over`.
            order_by: Column(s) to order window functions by.
                For lazy backends, this argument is required when `over` is applied
                to order-dependent functions, see [order-dependence](../concepts/order_dependence.md).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 4], "b": ["x", "x", "y"]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_min_per_group=nw.col("a").min().over("b"))
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |   a  b  a_min_per_group|
            |0  1  x                1|
            |1  2  x                1|
            |2  4  y                4|
            └────────────────────────┘

            Cumulative operations are also supported, but (currently) only for
            pandas and Polars:

            >>> df.with_columns(a_cum_sum_per_group=nw.col("a").cum_sum().over("b"))
            ┌────────────────────────────┐
            |     Narwhals DataFrame     |
            |----------------------------|
            |   a  b  a_cum_sum_per_group|
            |0  1  x                    1|
            |1  2  x                    3|
            |2  4  y                    4|
            └────────────────────────────┘
        """
        flat_partition_by = flatten(partition_by)
        flat_order_by = [order_by] if isinstance(order_by, str) else (order_by or [])
        if not flat_partition_by and not flat_order_by:  # pragma: no cover
            msg = "At least one of `partition_by` or `order_by` must be specified."
            raise ValueError(msg)

        current_meta = self._metadata
        if flat_order_by:
            next_meta = current_meta.with_ordered_over()
        elif not flat_partition_by:  # pragma: no cover
            msg = "At least one of `partition_by` or `order_by` must be specified."
            raise InvalidOperationError(msg)
        else:
            next_meta = current_meta.with_partitioned_over()

        return self.__class__(
            lambda plx: self._to_compliant_expr(plx).over(
                flat_partition_by, flat_order_by
            ),
            next_meta,
        )

    def is_duplicated(self) -> Self:
        r"""Return a boolean mask indicating duplicated values.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 1], "b": ["a", "a", "b", "c"]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.all().is_duplicated().name.suffix("_is_duplicated"))
            ┌─────────────────────────────────────────┐
            |           Narwhals DataFrame            |
            |-----------------------------------------|
            |   a  b  a_is_duplicated  b_is_duplicated|
            |0  1  a             True             True|
            |1  2  a            False             True|
            |2  3  b            False            False|
            |3  1  c             True            False|
            └─────────────────────────────────────────┘
        """
        return self._with_window(lambda plx: self._to_compliant_expr(plx).is_duplicated())

    def is_unique(self) -> Self:
        r"""Return a boolean mask indicating unique values.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 1], "b": ["a", "a", "b", "c"]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(nw.all().is_unique().name.suffix("_is_unique"))
            ┌─────────────────────────────────┐
            |       Narwhals DataFrame        |
            |---------------------------------|
            |   a  b  a_is_unique  b_is_unique|
            |0  1  a        False        False|
            |1  2  a         True        False|
            |2  3  b         True         True|
            |3  1  c        False         True|
            └─────────────────────────────────┘
        """
        return self._with_window(lambda plx: self._to_compliant_expr(plx).is_unique())

    def null_count(self) -> Self:
        r"""Count null values.

        Notes:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": [1, 2, None, 1], "b": ["a", None, "b", None]}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.all().null_count())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a  b      |
            |     0  1  2      |
            └──────────────────┘
        """
        return self._with_aggregation(
            lambda plx: self._to_compliant_expr(plx).null_count()
        )

    def is_first_distinct(self) -> Self:
        r"""Return a boolean mask indicating the first occurrence of each distinct value.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 1], "b": ["a", "a", "b", "c"]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.all().is_first_distinct().name.suffix("_is_first_distinct")
            ... )
            ┌─────────────────────────────────────────────────┐
            |               Narwhals DataFrame                |
            |-------------------------------------------------|
            |   a  b  a_is_first_distinct  b_is_first_distinct|
            |0  1  a                 True                 True|
            |1  2  a                 True                False|
            |2  3  b                 True                 True|
            |3  1  c                False                 True|
            └─────────────────────────────────────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).is_first_distinct()
        )

    def is_last_distinct(self) -> Self:
        r"""Return a boolean mask indicating the last occurrence of each distinct value.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3, 1], "b": ["a", "a", "b", "c"]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.all().is_last_distinct().name.suffix("_is_last_distinct")
            ... )
            ┌───────────────────────────────────────────────┐
            |              Narwhals DataFrame               |
            |-----------------------------------------------|
            |   a  b  a_is_last_distinct  b_is_last_distinct|
            |0  1  a               False               False|
            |1  2  a                True                True|
            |2  3  b                True                True|
            |3  1  c                True                True|
            └───────────────────────────────────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).is_last_distinct()
        )

    def quantile(
        self, quantile: float, interpolation: RollingInterpolationMethod
    ) -> Self:
        r"""Get quantile value.

        Arguments:
            quantile: Quantile between 0.0 and 1.0.
            interpolation: Interpolation method.

        Note:
            - pandas and Polars may have implementation differences for a given interpolation method.
            - [dask](https://docs.dask.org/en/stable/generated/dask.dataframe.Series.quantile.html) has
                its own method to approximate quantile and it doesn't implement 'nearest', 'higher',
                'lower', 'midpoint' as interpolation method - use 'linear' which is closest to the
                native 'dask' - method.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame(
            ...     {"a": list(range(50)), "b": list(range(50, 100))}
            ... )
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a", "b").quantile(0.5, interpolation="linear"))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |        a     b   |
            |  0  24.5  74.5   |
            └──────────────────┘
        """
        return self._with_aggregation(
            lambda plx: self._to_compliant_expr(plx).quantile(quantile, interpolation)
        )

    def round(self, decimals: int = 0) -> Self:
        r"""Round underlying floating point data by `decimals` digits.

        Arguments:
            decimals: Number of decimals to round by.


        Notes:
            For values exactly halfway between rounded decimal values pandas behaves differently than Polars and Arrow.

            pandas rounds to the nearest even value (e.g. -0.5 and 0.5 round to 0.0, 1.5 and 2.5 round to 2.0, 3.5 and
            4.5 to 4.0, etc..).

            Polars and Arrow round away from 0 (e.g. -0.5 to -1.0, 0.5 to 1.0, 1.5 to 2.0, 2.5 to 3.0, etc..).

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1.12345, 2.56789, 3.901234]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_rounded=nw.col("a").round(1))
            ┌──────────────────────┐
            |  Narwhals DataFrame  |
            |----------------------|
            |          a  a_rounded|
            |0  1.123450        1.1|
            |1  2.567890        2.6|
            |2  3.901234        3.9|
            └──────────────────────┘
        """
        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).round(decimals)
        )

    def floor(self) -> Self:
        r"""Compute the numerical floor.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"values": [1.1, 4.3, -1.3]})
            >>> df = nw.from_native(df_native)
            >>> result = df.with_columns(floor=nw.col("values").floor())
            >>> result
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |pyarrow.Table           |
            |values: double          |
            |floor: double           |
            |----                    |
            |values: [[1.1,4.3,-1.3]]|
            |floor: [[1,4,-2]]       |
            └────────────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).floor())

    def ceil(self) -> Self:
        r"""Compute the numerical ceiling.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"values": [1.1, 4.3, -1.3]})
            >>> df = nw.from_native(df_native)
            >>> result = df.with_columns(ceil=nw.col("values").ceil())
            >>> result
            ┌────────────────────────┐
            |   Narwhals DataFrame   |
            |------------------------|
            |pyarrow.Table           |
            |values: double          |
            |ceil: double            |
            |----                    |
            |values: [[1.1,4.3,-1.3]]|
            |ceil: [[2,5,-1]]        |
            └────────────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).ceil())

    def len(self) -> Self:
        r"""Return the number of elements in the column.

        Null values count towards the total.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": ["x", "y", "z"], "b": [1, 2, 1]})
            >>> df = nw.from_native(df_native)
            >>> df.select(
            ...     nw.col("a").filter(nw.col("b") == 1).len().alias("a1"),
            ...     nw.col("a").filter(nw.col("b") == 2).len().alias("a2"),
            ... )
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |       a1  a2     |
            |    0   2   1     |
            └──────────────────┘
        """
        return self._with_aggregation(lambda plx: self._to_compliant_expr(plx).len())

    def clip(
        self,
        lower_bound: IntoExpr | NumericLiteral | TemporalLiteral | None = None,
        upper_bound: IntoExpr | NumericLiteral | TemporalLiteral | None = None,
    ) -> Self:
        r"""Clip values in the Series.

        Arguments:
            lower_bound: Lower bound value. String literals are treated as column names.
            upper_bound: Upper bound value. String literals are treated as column names.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 2, 3]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_clipped=nw.col("a").clip(-1, 3))
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |    a  a_clipped  |
            | 0  1          1  |
            | 1  2          2  |
            | 2  3          3  |
            └──────────────────┘
        """
        return self._with_nary(
            lambda *exprs: exprs[0].clip(
                exprs[1] if lower_bound is not None else None,
                exprs[2] if upper_bound is not None else None,
            ),
            lower_bound,
            upper_bound,
        )

    def first(self) -> Self:
        """Get the first value.

        Notes:
            For lazy backends, this can only be used with `over`. We may introduce
            `min_by` in the future so it can be used as an aggregation.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> data = {"a": [1, 1, 2, 2], "b": ["foo", None, None, "baz"]}
            >>> df_native = pd.DataFrame(data)
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.all().first())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |       a    b     |
            |    0  1  foo     |
            └──────────────────┘

            >>> df.group_by("a").agg(nw.col("b").first())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |       a     b    |
            |    0  1   foo    |
            |    1  2  None    |
            └──────────────────┘
        """
        return self._with_orderable_aggregation(
            lambda plx: self._to_compliant_expr(plx).first()
        )

    def last(self) -> Self:
        """Get the last value.

        Notes:
            For lazy backends, this can only be used with `over`. We may introduce
            `max_by` in the future so it can be used as an aggregation.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> data = {"a": [1, 1, 2, 2], "b": ["foo", None, None, "baz"]}
            >>> df_native = pa.table(data)
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.all().last())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |  pyarrow.Table   |
            |  a: int64        |
            |  b: string       |
            |  ----            |
            |  a: [[2]]        |
            |  b: [["baz"]]    |
            └──────────────────┘

            >>> df.group_by("a").agg(nw.col("b").last())
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |pyarrow.Table     |
            |a: int64          |
            |b: string         |
            |----              |
            |a: [[1,2]]        |
            |b: [[null,"baz"]] |
            └──────────────────┘
        """
        return self._with_orderable_aggregation(
            lambda plx: self._to_compliant_expr(plx).last()
        )

    def mode(self, *, keep: ModeKeepStrategy = "all") -> Self:
        r"""Compute the most occurring value(s).

        Can return multiple values.

        Arguments:
            keep: Whether to keep all modes or any mode found. Remark that `keep='any'`
                is not deterministic for multimodal values.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 1, 2, 3], "b": [1, 1, 2, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.select(nw.col("a").mode()).sort("a")
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |          a       |
            |       0  1       |
            └──────────────────┘
        """
        _supported_keep_values = ("all", "any")
        if keep not in _supported_keep_values:  # pragma: no cover
            msg = f"`keep` must be one of {_supported_keep_values}, found '{keep}'"
            raise ValueError(msg)

        def compliant_expr(plx: Any) -> Any:
            return self._to_compliant_expr(plx).mode(keep=keep)

        if keep == "any":
            return self._with_aggregation(compliant_expr)
        return self._with_filtration(compliant_expr)

    def is_finite(self) -> Self:
        """Returns boolean values indicating which original values are finite.

        Warning:
            pandas handles null values differently from Polars and PyArrow.
            See [null_handling](../concepts/null_handling.md/) for reference.
            `is_finite` will return False for NaN and Null's in the Dask and
            pandas non-nullable backend, while for Polars, PyArrow and pandas
            nullable backends null values are kept as such.

        Examples:
            >>> import polars as pl
            >>> import narwhals as nw
            >>> df_native = pl.DataFrame({"a": [float("nan"), float("inf"), 2.0, None]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(a_is_finite=nw.col("a").is_finite())
            ┌──────────────────────┐
            |  Narwhals DataFrame  |
            |----------------------|
            |shape: (4, 2)         |
            |┌──────┬─────────────┐|
            |│ a    ┆ a_is_finite │|
            |│ ---  ┆ ---         │|
            |│ f64  ┆ bool        │|
            |╞══════╪═════════════╡|
            |│ NaN  ┆ false       │|
            |│ inf  ┆ false       │|
            |│ 2.0  ┆ true        │|
            |│ null ┆ null        │|
            |└──────┴─────────────┘|
            └──────────────────────┘
        """
        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).is_finite()
        )

    def cum_count(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative count of the non-null values in the column.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": ["x", "k", None, "d"]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").cum_count().alias("a_cum_count"),
            ...     nw.col("a").cum_count(reverse=True).alias("a_cum_count_reverse"),
            ... )
            ┌─────────────────────────────────────────┐
            |           Narwhals DataFrame            |
            |-----------------------------------------|
            |      a  a_cum_count  a_cum_count_reverse|
            |0     x            1                    3|
            |1     k            2                    2|
            |2  None            2                    1|
            |3     d            3                    1|
            └─────────────────────────────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).cum_count(reverse=reverse)
        )

    def cum_min(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative min of the non-null values in the column.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [3, 1, None, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").cum_min().alias("a_cum_min"),
            ...     nw.col("a").cum_min(reverse=True).alias("a_cum_min_reverse"),
            ... )
            ┌────────────────────────────────────┐
            |         Narwhals DataFrame         |
            |------------------------------------|
            |     a  a_cum_min  a_cum_min_reverse|
            |0  3.0        3.0                1.0|
            |1  1.0        1.0                1.0|
            |2  NaN        NaN                NaN|
            |3  2.0        1.0                2.0|
            └────────────────────────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).cum_min(reverse=reverse)
        )

    def cum_max(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative max of the non-null values in the column.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 3, None, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").cum_max().alias("a_cum_max"),
            ...     nw.col("a").cum_max(reverse=True).alias("a_cum_max_reverse"),
            ... )
            ┌────────────────────────────────────┐
            |         Narwhals DataFrame         |
            |------------------------------------|
            |     a  a_cum_max  a_cum_max_reverse|
            |0  1.0        1.0                3.0|
            |1  3.0        3.0                3.0|
            |2  NaN        NaN                NaN|
            |3  2.0        3.0                2.0|
            └────────────────────────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).cum_max(reverse=reverse)
        )

    def cum_prod(self, *, reverse: bool = False) -> Self:
        r"""Return the cumulative product of the non-null values in the column.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            reverse: reverse the operation

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1, 3, None, 2]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     nw.col("a").cum_prod().alias("a_cum_prod"),
            ...     nw.col("a").cum_prod(reverse=True).alias("a_cum_prod_reverse"),
            ... )
            ┌──────────────────────────────────────┐
            |          Narwhals DataFrame          |
            |--------------------------------------|
            |     a  a_cum_prod  a_cum_prod_reverse|
            |0  1.0         1.0                 6.0|
            |1  3.0         3.0                 6.0|
            |2  NaN         NaN                 NaN|
            |3  2.0         6.0                 2.0|
            └──────────────────────────────────────┘
        """
        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).cum_prod(reverse=reverse)
        )

    def rolling_sum(
        self, window_size: int, *, min_samples: int | None = None, center: bool = False
    ) -> Self:
        """Apply a rolling sum (moving sum) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their sum.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`
            center: Set the labels at the center of the window.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1.0, 2.0, None, 4.0]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_rolling_sum=nw.col("a").rolling_sum(window_size=3, min_samples=1)
            ... )
            ┌─────────────────────┐
            | Narwhals DataFrame  |
            |---------------------|
            |     a  a_rolling_sum|
            |0  1.0            1.0|
            |1  2.0            3.0|
            |2  NaN            3.0|
            |3  4.0            6.0|
            └─────────────────────┘
        """
        window_size, min_samples_int = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).rolling_sum(
                window_size=window_size, min_samples=min_samples_int, center=center
            )
        )

    def rolling_mean(
        self, window_size: int, *, min_samples: int | None = None, center: bool = False
    ) -> Self:
        """Apply a rolling mean (moving mean) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their mean.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`
            center: Set the labels at the center of the window.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1.0, 2.0, None, 4.0]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_rolling_mean=nw.col("a").rolling_mean(window_size=3, min_samples=1)
            ... )
            ┌──────────────────────┐
            |  Narwhals DataFrame  |
            |----------------------|
            |     a  a_rolling_mean|
            |0  1.0             1.0|
            |1  2.0             1.5|
            |2  NaN             1.5|
            |3  4.0             3.0|
            └──────────────────────┘
        """
        window_size, min_samples = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).rolling_mean(
                window_size=window_size, min_samples=min_samples, center=center
            )
        )

    def rolling_var(
        self,
        window_size: int,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Self:
        """Apply a rolling variance (moving variance) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their variance.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`.
            center: Set the labels at the center of the window.
            ddof: Delta Degrees of Freedom; the divisor for a length N window is N - ddof.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1.0, 2.0, None, 4.0]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_rolling_var=nw.col("a").rolling_var(window_size=3, min_samples=1)
            ... )
            ┌─────────────────────┐
            | Narwhals DataFrame  |
            |---------------------|
            |     a  a_rolling_var|
            |0  1.0            NaN|
            |1  2.0            0.5|
            |2  NaN            0.5|
            |3  4.0            2.0|
            └─────────────────────┘
        """
        window_size, min_samples = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).rolling_var(
                window_size=window_size, min_samples=min_samples, center=center, ddof=ddof
            )
        )

    def rolling_std(
        self,
        window_size: int,
        *,
        min_samples: int | None = None,
        center: bool = False,
        ddof: int = 1,
    ) -> Self:
        """Apply a rolling standard deviation (moving standard deviation) over the values.

        A window of length `window_size` will traverse the values. The resulting values
        will be aggregated to their standard deviation.

        The window at a given row will include the row itself and the `window_size - 1`
        elements before it.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            window_size: The length of the window in number of elements. It must be a
                strictly positive integer.
            min_samples: The number of values in the window that should be non-null before
                computing a result. If set to `None` (default), it will be set equal to
                `window_size`. If provided, it must be a strictly positive integer, and
                less than or equal to `window_size`.
            center: Set the labels at the center of the window.
            ddof: Delta Degrees of Freedom; the divisor for a length N window is N - ddof.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [1.0, 2.0, None, 4.0]})
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     a_rolling_std=nw.col("a").rolling_std(window_size=3, min_samples=1)
            ... )
            ┌─────────────────────┐
            | Narwhals DataFrame  |
            |---------------------|
            |     a  a_rolling_std|
            |0  1.0            NaN|
            |1  2.0       0.707107|
            |2  NaN       0.707107|
            |3  4.0       1.414214|
            └─────────────────────┘
        """
        window_size, min_samples = _validate_rolling_arguments(
            window_size=window_size, min_samples=min_samples
        )

        return self._with_orderable_window(
            lambda plx: self._to_compliant_expr(plx).rolling_std(
                window_size=window_size, min_samples=min_samples, center=center, ddof=ddof
            )
        )

    def rank(self, method: RankMethod = "average", *, descending: bool = False) -> Self:
        """Assign ranks to data, dealing with ties appropriately.

        Notes:
            The resulting dtype may differ between backends.

        Info:
            For lazy backends, this operation must be followed by `Expr.over` with
            `order_by` specified, see [order-dependence](../concepts/order_dependence.md).

        Arguments:
            method: The method used to assign ranks to tied elements.
                The following methods are available (default is 'average')

                - *"average"*: The average of the ranks that would have been assigned to
                    all the tied values is assigned to each value.
                - *"min"*: The minimum of the ranks that would have been assigned to all
                    the tied values is assigned to each value. (This is also referred to
                    as "competition" ranking.)
                - *"max"*: The maximum of the ranks that would have been assigned to all
                    the tied values is assigned to each value.
                - *"dense"*: Like "min", but the rank of the next highest element is
                    assigned the rank immediately after those assigned to the tied elements.
                - *"ordinal"*: All values are given a distinct rank, corresponding to the
                    order that the values occur in the Series.

            descending: Rank in descending order.

        Examples:
            >>> import pandas as pd
            >>> import narwhals as nw
            >>> df_native = pd.DataFrame({"a": [3, 6, 1, 1, 6]})
            >>> df = nw.from_native(df_native)
            >>> result = df.with_columns(rank=nw.col("a").rank(method="dense"))
            >>> result
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |       a  rank    |
            |    0  3   2.0    |
            |    1  6   3.0    |
            |    2  1   1.0    |
            |    3  1   1.0    |
            |    4  6   3.0    |
            └──────────────────┘
        """
        supported_rank_methods = {"average", "min", "max", "dense", "ordinal"}
        if method not in supported_rank_methods:
            msg = (
                "Ranking method must be one of {'average', 'min', 'max', 'dense', 'ordinal'}. "
                f"Found '{method}'"
            )
            raise ValueError(msg)

        return self._with_window(
            lambda plx: self._to_compliant_expr(plx).rank(
                method=method, descending=descending
            )
        )

    def log(self, base: float = math.e) -> Self:
        r"""Compute the logarithm to a given base.

        Arguments:
            base: Given base, defaults to `e`

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"values": [1, 2, 4]})
            >>> df = nw.from_native(df_native)
            >>> result = df.with_columns(
            ...     log=nw.col("values").log(), log_2=nw.col("values").log(base=2)
            ... )
            >>> result
            ┌────────────────────────────────────────────────┐
            |               Narwhals DataFrame               |
            |------------------------------------------------|
            |pyarrow.Table                                   |
            |values: int64                                   |
            |log: double                                     |
            |log_2: double                                   |
            |----                                            |
            |values: [[1,2,4]]                               |
            |log: [[0,0.6931471805599453,1.3862943611198906]]|
            |log_2: [[0,1,2]]                                |
            └────────────────────────────────────────────────┘
        """
        return self._with_elementwise(
            lambda plx: self._to_compliant_expr(plx).log(base=base)
        )

    def exp(self) -> Self:
        r"""Compute the exponent.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"values": [-1, 0, 1]})
            >>> df = nw.from_native(df_native)
            >>> result = df.with_columns(exp=nw.col("values").exp())
            >>> result
            ┌────────────────────────────────────────────────┐
            |               Narwhals DataFrame               |
            |------------------------------------------------|
            |pyarrow.Table                                   |
            |values: int64                                   |
            |exp: double                                     |
            |----                                            |
            |values: [[-1,0,1]]                              |
            |exp: [[0.36787944117144233,1,2.718281828459045]]|
            └────────────────────────────────────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).exp())

    def sqrt(self) -> Self:
        r"""Compute the square root.

        Examples:
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>> df_native = pa.table({"values": [1, 4, 9]})
            >>> df = nw.from_native(df_native)
            >>> result = df.with_columns(sqrt=nw.col("values").sqrt())
            >>> result
            ┌──────────────────┐
            |Narwhals DataFrame|
            |------------------|
            |pyarrow.Table     |
            |values: int64     |
            |sqrt: double      |
            |----              |
            |values: [[1,4,9]] |
            |sqrt: [[1,2,3]]   |
            └──────────────────┘
        """
        return self._with_elementwise(lambda plx: self._to_compliant_expr(plx).sqrt())

    def is_close(
        self,
        other: Self | NumericLiteral,
        *,
        abs_tol: float = 0.0,
        rel_tol: float = 1e-09,
        nans_equal: bool = False,
    ) -> Self:
        r"""Check if this expression is close, i.e. almost equal, to the other expression.

        Two values `a` and `b` are considered close if the following condition holds:

        $$
        |a-b| \le max \{ \text{rel\_tol} \cdot max \{ |a|, |b| \}, \text{abs\_tol} \}
        $$

        Arguments:
            other: Values to compare with.
            abs_tol: Absolute tolerance. This is the maximum allowed absolute difference
                between two values. Must be non-negative.
            rel_tol: Relative tolerance. This is the maximum allowed difference between
                two values, relative to the larger absolute value. Must be in the range
                [0, 1).
            nans_equal: Whether NaN values should be considered equal.

        Notes:
            The implementation of this method is symmetric and mirrors the behavior of
            `math.isclose`. Specifically note that this behavior is different to
            `numpy.isclose`.

        Examples:
            >>> import duckdb
            >>> import pyarrow as pa
            >>> import narwhals as nw
            >>>
            >>> data = {
            ...     "x": [1.0, float("inf"), 1.41, None, float("nan")],
            ...     "y": [1.2, float("inf"), 1.40, None, float("nan")],
            ... }
            >>> _table = pa.table(data)
            >>> df_native = duckdb.table("_table")
            >>> df = nw.from_native(df_native)
            >>> df.with_columns(
            ...     is_close=nw.col("x").is_close(
            ...         nw.col("y"), abs_tol=0.1, nans_equal=True
            ...     )
            ... )
            ┌──────────────────────────────┐
            |      Narwhals LazyFrame      |
            |------------------------------|
            |┌────────┬────────┬──────────┐|
            |│   x    │   y    │ is_close │|
            |│ double │ double │ boolean  │|
            |├────────┼────────┼──────────┤|
            |│    1.0 │    1.2 │ false    │|
            |│    inf │    inf │ true     │|
            |│   1.41 │    1.4 │ true     │|
            |│   NULL │   NULL │ NULL     │|
            |│    nan │    nan │ true     │|
            |└────────┴────────┴──────────┘|
            └──────────────────────────────┘
        """
        if abs_tol < 0:
            msg = f"`abs_tol` must be non-negative but got {abs_tol}"
            raise ComputeError(msg)

        if not (0 <= rel_tol < 1):
            msg = f"`rel_tol` must be in the range [0, 1) but got {rel_tol}"
            raise ComputeError(msg)

        kwargs = {"abs_tol": abs_tol, "rel_tol": rel_tol, "nans_equal": nans_equal}
        return self._with_nary(
            lambda *exprs: exprs[0].is_close(exprs[1], **kwargs), other
        )

    @property
    def str(self) -> ExprStringNamespace[Self]:
        return ExprStringNamespace(self)

    @property
    def dt(self) -> ExprDateTimeNamespace[Self]:
        return ExprDateTimeNamespace(self)

    @property
    def cat(self) -> ExprCatNamespace[Self]:
        return ExprCatNamespace(self)

    @property
    def name(self) -> ExprNameNamespace[Self]:
        return ExprNameNamespace(self)

    @property
    def list(self) -> ExprListNamespace[Self]:
        return ExprListNamespace(self)

    @property
    def struct(self) -> ExprStructNamespace[Self]:
        return ExprStructNamespace(self)


__all__ = ["Expr"]
