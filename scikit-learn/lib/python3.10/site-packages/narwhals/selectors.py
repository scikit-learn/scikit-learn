from __future__ import annotations

from typing import TYPE_CHECKING, Any, NoReturn

from narwhals._expression_parsing import ExprMetadata, combine_metadata
from narwhals._utils import flatten
from narwhals.expr import Expr

if TYPE_CHECKING:
    from collections.abc import Iterable
    from datetime import timezone

    from narwhals.dtypes import DType
    from narwhals.typing import TimeUnit


class Selector(Expr):
    def _to_expr(self) -> Expr:
        return Expr(self._to_compliant_expr, self._metadata)

    def __add__(self, other: Any) -> Expr:  # type: ignore[override]
        if isinstance(other, Selector):
            msg = "unsupported operand type(s) for op: ('Selector' + 'Selector')"
            raise TypeError(msg)
        return self._to_expr() + other  # type: ignore[no-any-return]

    def __or__(self, other: Any) -> Expr:  # type: ignore[override]
        if isinstance(other, Selector):
            return self.__class__(
                lambda plx: self._to_compliant_expr(plx) | other._to_compliant_expr(plx),
                combine_metadata(
                    self,
                    other,
                    str_as_lit=False,
                    allow_multi_output=True,
                    to_single_output=False,
                ),
            )
        return self._to_expr() | other  # type: ignore[no-any-return]

    def __and__(self, other: Any) -> Expr:  # type: ignore[override]
        if isinstance(other, Selector):
            return self.__class__(
                lambda plx: self._to_compliant_expr(plx) & other._to_compliant_expr(plx),
                combine_metadata(
                    self,
                    other,
                    str_as_lit=False,
                    allow_multi_output=True,
                    to_single_output=False,
                ),
            )
        return self._to_expr() & other  # type: ignore[no-any-return]

    def __rsub__(self, other: Any) -> NoReturn:
        raise NotImplementedError

    def __rand__(self, other: Any) -> NoReturn:
        raise NotImplementedError

    def __ror__(self, other: Any) -> NoReturn:
        raise NotImplementedError


def by_dtype(*dtypes: DType | type[DType] | Iterable[DType | type[DType]]) -> Selector:
    """Select columns based on their dtype.

    Arguments:
        dtypes: one or data types to select

    Examples:
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pa.table({"a": [1, 2], "b": ["x", "y"], "c": [4.1, 2.3]})
        >>> df = nw.from_native(df_native)

        Let's select int64 and float64  dtypes and multiply each value by 2:

        >>> df.select(ncs.by_dtype(nw.Int64, nw.Float64) * 2).to_native()
        pyarrow.Table
        a: int64
        c: double
        ----
        a: [[2,4]]
        c: [[8.2,4.6]]
    """
    flattened = flatten(dtypes)
    return Selector(
        lambda plx: plx.selectors.by_dtype(flattened),
        ExprMetadata.selector_multi_unnamed(),
    )


def matches(pattern: str) -> Selector:
    """Select all columns that match the given regex pattern.

    Arguments:
        pattern: A valid regular expression pattern.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pd.DataFrame(
        ...     {"bar": [123, 456], "baz": [2.0, 5.5], "zap": [0, 1]}
        ... )
        >>> df = nw.from_native(df_native)

        Let's select column names containing an 'a', preceded by a character that is not 'z':

        >>> df.select(ncs.matches("[^z]a")).to_native()
           bar  baz
        0  123  2.0
        1  456  5.5
    """
    return Selector(
        lambda plx: plx.selectors.matches(pattern), ExprMetadata.selector_multi_unnamed()
    )


def numeric() -> Selector:
    """Select numeric columns.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": ["x", "y"], "c": [4.1, 2.3]})
        >>> df = nw.from_native(df_native)

        Let's select numeric dtypes and multiply each value by 2:

        >>> df.select(ncs.numeric() * 2).to_native()
        shape: (2, 2)
        ┌─────┬─────┐
        │ a   ┆ c   │
        │ --- ┆ --- │
        │ i64 ┆ f64 │
        ╞═════╪═════╡
        │ 2   ┆ 8.2 │
        │ 4   ┆ 4.6 │
        └─────┴─────┘
    """
    return Selector(
        lambda plx: plx.selectors.numeric(), ExprMetadata.selector_multi_unnamed()
    )


def boolean() -> Selector:
    """Select boolean columns.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": ["x", "y"], "c": [False, True]})
        >>> df = nw.from_native(df_native)

        Let's select boolean dtypes:

        >>> df.select(ncs.boolean())
        ┌──────────────────┐
        |Narwhals DataFrame|
        |------------------|
        |  shape: (2, 1)   |
        |  ┌───────┐       |
        |  │ c     │       |
        |  │ ---   │       |
        |  │ bool  │       |
        |  ╞═══════╡       |
        |  │ false │       |
        |  │ true  │       |
        |  └───────┘       |
        └──────────────────┘
    """
    return Selector(
        lambda plx: plx.selectors.boolean(), ExprMetadata.selector_multi_unnamed()
    )


def string() -> Selector:
    """Select string columns.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": ["x", "y"], "c": [False, True]})
        >>> df = nw.from_native(df_native)

        Let's select string dtypes:

        >>> df.select(ncs.string()).to_native()
        shape: (2, 1)
        ┌─────┐
        │ b   │
        │ --- │
        │ str │
        ╞═════╡
        │ x   │
        │ y   │
        └─────┘
    """
    return Selector(
        lambda plx: plx.selectors.string(), ExprMetadata.selector_multi_unnamed()
    )


def categorical() -> Selector:
    """Select categorical columns.

    Examples:
        >>> import polars as pl
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pl.DataFrame({"a": [1, 2], "b": ["x", "y"], "c": [False, True]})

        Let's convert column "b" to categorical, and then select categorical dtypes:

        >>> df = nw.from_native(df_native).with_columns(
        ...     b=nw.col("b").cast(nw.Categorical())
        ... )
        >>> df.select(ncs.categorical()).to_native()
        shape: (2, 1)
        ┌─────┐
        │ b   │
        │ --- │
        │ cat │
        ╞═════╡
        │ x   │
        │ y   │
        └─────┘
    """
    return Selector(
        lambda plx: plx.selectors.categorical(), ExprMetadata.selector_multi_unnamed()
    )


def all() -> Selector:
    """Select all columns.

    Examples:
        >>> import pandas as pd
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>> df_native = pd.DataFrame({"a": [1, 2], "b": ["x", "y"], "c": [False, True]})
        >>> df = nw.from_native(df_native)

        Let's select all dtypes:

        >>> df.select(ncs.all()).to_native()
           a  b      c
        0  1  x  False
        1  2  y   True
    """
    return Selector(
        lambda plx: plx.selectors.all(), ExprMetadata.selector_multi_unnamed()
    )


def datetime(
    time_unit: TimeUnit | Iterable[TimeUnit] | None = None,
    time_zone: str | timezone | Iterable[str | timezone | None] | None = ("*", None),
) -> Selector:
    """Select all datetime columns, optionally filtering by time unit/zone.

    Arguments:
        time_unit: One (or more) of the allowed timeunit precision strings, "ms", "us",
            "ns" and "s". Omit to select columns with any valid timeunit.
        time_zone: Specify which timezone(s) to select

            * One or more timezone strings, as defined in zoneinfo (to see valid options
                run `import zoneinfo; zoneinfo.available_timezones()` for a full list).
            * Set `None` to select Datetime columns that do not have a timezone.
            * Set `"*"` to select Datetime columns that have *any* timezone.

    Examples:
        >>> from datetime import datetime, timezone
        >>> import pyarrow as pa
        >>> import narwhals as nw
        >>> import narwhals.selectors as ncs
        >>>
        >>> utc_tz = timezone.utc
        >>> data = {
        ...     "tstamp_utc": [
        ...         datetime(2023, 4, 10, 12, 14, 16, 999000, tzinfo=utc_tz),
        ...         datetime(2025, 8, 25, 14, 18, 22, 666000, tzinfo=utc_tz),
        ...     ],
        ...     "tstamp": [
        ...         datetime(2000, 11, 20, 18, 12, 16, 600000),
        ...         datetime(2020, 10, 30, 10, 20, 25, 123000),
        ...     ],
        ...     "numeric": [3.14, 6.28],
        ... }
        >>> df_native = pa.table(data)
        >>> df_nw = nw.from_native(df_native)
        >>> df_nw.select(ncs.datetime()).to_native()
        pyarrow.Table
        tstamp_utc: timestamp[us, tz=UTC]
        tstamp: timestamp[us]
        ----
        tstamp_utc: [[2023-04-10 12:14:16.999000Z,2025-08-25 14:18:22.666000Z]]
        tstamp: [[2000-11-20 18:12:16.600000,2020-10-30 10:20:25.123000]]

        Select only datetime columns that have any time_zone specification:

        >>> df_nw.select(ncs.datetime(time_zone="*")).to_native()
        pyarrow.Table
        tstamp_utc: timestamp[us, tz=UTC]
        ----
        tstamp_utc: [[2023-04-10 12:14:16.999000Z,2025-08-25 14:18:22.666000Z]]
    """
    return Selector(
        lambda plx: plx.selectors.datetime(time_unit=time_unit, time_zone=time_zone),
        ExprMetadata.selector_multi_unnamed(),
    )


__all__ = [
    "all",
    "boolean",
    "by_dtype",
    "categorical",
    "datetime",
    "matches",
    "numeric",
    "string",
]
