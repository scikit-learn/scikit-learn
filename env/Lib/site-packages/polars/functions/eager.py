from __future__ import annotations

import contextlib
from collections.abc import Sequence
from functools import reduce
from itertools import chain
from typing import TYPE_CHECKING, get_args

import polars._reexport as pl
from polars import functions as F
from polars._typing import ConcatMethod
from polars._utils.various import ordered_unique
from polars._utils.wrap import wrap_df, wrap_expr, wrap_ldf, wrap_s
from polars.exceptions import InvalidOperationError

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars.polars as plr

if TYPE_CHECKING:
    from collections.abc import Iterable

    from polars import DataFrame, Expr, LazyFrame, Series
    from polars._typing import FrameType, JoinStrategy, PolarsType


def concat(
    items: Iterable[PolarsType],
    *,
    how: ConcatMethod = "vertical",
    rechunk: bool = False,
    parallel: bool = True,
) -> PolarsType:
    """
    Combine multiple DataFrames, LazyFrames, or Series into a single object.

    Parameters
    ----------
    items
        DataFrames, LazyFrames, or Series to concatenate.
    how : {'vertical', 'vertical_relaxed', 'diagonal', 'diagonal_relaxed', 'horizontal', 'align'}
        Series only support the `vertical` strategy.

        * vertical: Applies multiple `vstack` operations.
        * vertical_relaxed: Same as `vertical`, but additionally coerces columns to
          their common supertype *if* they are mismatched (eg: Int32 → Int64).
        * diagonal: Finds a union between the column schemas and fills missing column
          values with `null`.
        * diagonal_relaxed: Same as `diagonal`, but additionally coerces columns to
          their common supertype *if* they are mismatched (eg: Int32 → Int64).
        * horizontal: Stacks Series from DataFrames horizontally and fills with `null`
          if the lengths don't match.
        * align: Combines frames horizontally, auto-determining the common key columns
          and aligning rows using the same logic as `align_frames`; this behaviour is
          patterned after a full outer join, but does not handle column-name collision.
          (If you need more control, you should use a suitable join method instead).
    rechunk
        Make sure that the result data is in contiguous memory.
    parallel
        Only relevant for LazyFrames. This determines if the concatenated
        lazy computations may be executed in parallel.

    Examples
    --------
    >>> df1 = pl.DataFrame({"a": [1], "b": [3]})
    >>> df2 = pl.DataFrame({"a": [2], "b": [4]})
    >>> pl.concat([df1, df2])  # default is 'vertical' strategy
    shape: (2, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ i64 ┆ i64 │
    ╞═════╪═════╡
    │ 1   ┆ 3   │
    │ 2   ┆ 4   │
    └─────┴─────┘

    >>> df1 = pl.DataFrame({"a": [1], "b": [3]})
    >>> df2 = pl.DataFrame({"a": [2.5], "b": [4]})
    >>> pl.concat([df1, df2], how="vertical_relaxed")  # 'a' coerced into f64
    shape: (2, 2)
    ┌─────┬─────┐
    │ a   ┆ b   │
    │ --- ┆ --- │
    │ f64 ┆ i64 │
    ╞═════╪═════╡
    │ 1.0 ┆ 3   │
    │ 2.5 ┆ 4   │
    └─────┴─────┘

    >>> df_h1 = pl.DataFrame({"l1": [1, 2], "l2": [3, 4]})
    >>> df_h2 = pl.DataFrame({"r1": [5, 6], "r2": [7, 8], "r3": [9, 10]})
    >>> pl.concat([df_h1, df_h2], how="horizontal")
    shape: (2, 5)
    ┌─────┬─────┬─────┬─────┬─────┐
    │ l1  ┆ l2  ┆ r1  ┆ r2  ┆ r3  │
    │ --- ┆ --- ┆ --- ┆ --- ┆ --- │
    │ i64 ┆ i64 ┆ i64 ┆ i64 ┆ i64 │
    ╞═════╪═════╪═════╪═════╪═════╡
    │ 1   ┆ 3   ┆ 5   ┆ 7   ┆ 9   │
    │ 2   ┆ 4   ┆ 6   ┆ 8   ┆ 10  │
    └─────┴─────┴─────┴─────┴─────┘

    >>> df_d1 = pl.DataFrame({"a": [1], "b": [3]})
    >>> df_d2 = pl.DataFrame({"a": [2], "c": [4]})
    >>> pl.concat([df_d1, df_d2], how="diagonal")
    shape: (2, 3)
    ┌─────┬──────┬──────┐
    │ a   ┆ b    ┆ c    │
    │ --- ┆ ---  ┆ ---  │
    │ i64 ┆ i64  ┆ i64  │
    ╞═════╪══════╪══════╡
    │ 1   ┆ 3    ┆ null │
    │ 2   ┆ null ┆ 4    │
    └─────┴──────┴──────┘

    >>> df_a1 = pl.DataFrame({"id": [1, 2], "x": [3, 4]})
    >>> df_a2 = pl.DataFrame({"id": [2, 3], "y": [5, 6]})
    >>> df_a3 = pl.DataFrame({"id": [1, 3], "z": [7, 8]})
    >>> pl.concat([df_a1, df_a2, df_a3], how="align")
    shape: (3, 4)
    ┌─────┬──────┬──────┬──────┐
    │ id  ┆ x    ┆ y    ┆ z    │
    │ --- ┆ ---  ┆ ---  ┆ ---  │
    │ i64 ┆ i64  ┆ i64  ┆ i64  │
    ╞═════╪══════╪══════╪══════╡
    │ 1   ┆ 3    ┆ null ┆ 7    │
    │ 2   ┆ 4    ┆ 5    ┆ null │
    │ 3   ┆ null ┆ 6    ┆ 8    │
    └─────┴──────┴──────┴──────┘
    """  # noqa: W505
    # unpack/standardise (handles generator input)
    elems = list(items)

    if not elems:
        msg = "cannot concat empty list"
        raise ValueError(msg)
    elif len(elems) == 1 and isinstance(
        elems[0], (pl.DataFrame, pl.Series, pl.LazyFrame)
    ):
        return elems[0]

    if how == "align":
        if not isinstance(elems[0], (pl.DataFrame, pl.LazyFrame)):
            msg = f"'align' strategy is not supported for {type(elems[0]).__name__!r}"
            raise TypeError(msg)

        # establish common columns, maintaining the order in which they appear
        all_columns = list(chain.from_iterable(e.collect_schema() for e in elems))
        key = {v: k for k, v in enumerate(ordered_unique(all_columns))}
        common_cols = sorted(
            reduce(
                lambda x, y: set(x) & set(y),  # type: ignore[arg-type, return-value]
                chain(e.collect_schema() for e in elems),
            ),
            key=lambda k: key.get(k, 0),
        )
        # we require at least one key column for 'align'
        if not common_cols:
            msg = "'align' strategy requires at least one common column"
            raise InvalidOperationError(msg)

        # align the frame data using a full outer join with no suffix-resolution
        # (so we raise an error in case of column collision, like "horizontal")
        lf: LazyFrame = reduce(
            lambda x, y: (
                x.join(
                    y,
                    how="full",
                    on=common_cols,
                    suffix="_PL_CONCAT_RIGHT",
                    maintain_order="right_left",
                )
                # Coalesce full outer join columns
                .with_columns(
                    F.coalesce([name, f"{name}_PL_CONCAT_RIGHT"])
                    for name in common_cols
                )
                .drop([f"{name}_PL_CONCAT_RIGHT" for name in common_cols])
            ),
            [df.lazy() for df in elems],
        ).sort(by=common_cols)

        eager = isinstance(elems[0], pl.DataFrame)
        return lf.collect() if eager else lf  # type: ignore[return-value]

    out: Series | DataFrame | LazyFrame | Expr
    first = elems[0]

    if isinstance(first, pl.DataFrame):
        if how == "vertical":
            out = wrap_df(plr.concat_df(elems))
        elif how == "vertical_relaxed":
            out = wrap_ldf(
                plr.concat_lf(
                    [df.lazy() for df in elems],
                    rechunk=rechunk,
                    parallel=parallel,
                    to_supertypes=True,
                )
            ).collect(no_optimization=True)

        elif how == "diagonal":
            out = wrap_df(plr.concat_df_diagonal(elems))
        elif how == "diagonal_relaxed":
            out = wrap_ldf(
                plr.concat_lf_diagonal(
                    [df.lazy() for df in elems],
                    rechunk=rechunk,
                    parallel=parallel,
                    to_supertypes=True,
                )
            ).collect(no_optimization=True)
        elif how == "horizontal":
            out = wrap_df(plr.concat_df_horizontal(elems))
        else:
            allowed = ", ".join(repr(m) for m in get_args(ConcatMethod))
            msg = f"DataFrame `how` must be one of {{{allowed}}}, got {how!r}"
            raise ValueError(msg)

    elif isinstance(first, pl.LazyFrame):
        if how in ("vertical", "vertical_relaxed"):
            return wrap_ldf(
                plr.concat_lf(
                    elems,
                    rechunk=rechunk,
                    parallel=parallel,
                    to_supertypes=how.endswith("relaxed"),
                )
            )
        elif how in ("diagonal", "diagonal_relaxed"):
            return wrap_ldf(
                plr.concat_lf_diagonal(
                    elems,
                    rechunk=rechunk,
                    parallel=parallel,
                    to_supertypes=how.endswith("relaxed"),
                )
            )
        elif how == "horizontal":
            return wrap_ldf(
                plr.concat_lf_horizontal(
                    elems,
                    parallel=parallel,
                )
            )
        else:
            allowed = ", ".join(repr(m) for m in get_args(ConcatMethod))
            msg = f"LazyFrame `how` must be one of {{{allowed}}}, got {how!r}"
            raise ValueError(msg)

    elif isinstance(first, pl.Series):
        if how == "vertical":
            out = wrap_s(plr.concat_series(elems))
        else:
            msg = "Series only supports 'vertical' concat strategy"
            raise ValueError(msg)

    elif isinstance(first, pl.Expr):
        return wrap_expr(plr.concat_expr([e._pyexpr for e in elems], rechunk))
    else:
        msg = f"did not expect type: {type(first).__name__!r} in `concat`"
        raise TypeError(msg)

    if rechunk:
        return out.rechunk()
    return out


def _alignment_join(
    *idx_frames: tuple[int, LazyFrame],
    align_on: list[str],
    how: JoinStrategy = "full",
    descending: bool | Sequence[bool] = False,
) -> LazyFrame:
    """Create a single master frame with all rows aligned on the common key values."""
    # note: can stackoverflow if the join becomes too large, so we
    # collect eagerly when hitting a large enough number of frames
    post_align_collect = len(idx_frames) >= 250

    def join_func(
        idx_x: tuple[int, LazyFrame],
        idx_y: tuple[int, LazyFrame],
    ) -> tuple[int, LazyFrame]:
        (_, x), (y_idx, y) = idx_x, idx_y
        return y_idx, x.join(
            y,
            how=how,
            on=align_on,
            suffix=f":{y_idx}",
            join_nulls=True,
            coalesce=True,
            maintain_order="right_left",
        )

    joined = reduce(join_func, idx_frames)[1].sort(by=align_on, descending=descending)
    if post_align_collect:
        joined = joined.collect(no_optimization=True).lazy()
    return joined


def align_frames(
    *frames: FrameType,
    on: str | Expr | Sequence[str] | Sequence[Expr] | Sequence[str | Expr],
    how: JoinStrategy = "full",
    select: str | Expr | Sequence[str | Expr] | None = None,
    descending: bool | Sequence[bool] = False,
) -> list[FrameType]:
    r"""
    Align a sequence of frames using common values from one or more columns as a key.

    Frames that do not contain the given key values have rows injected (with nulls
    filling the non-key columns), and each resulting frame is sorted by the key.

    The original column order of input frames is not changed unless `select` is
    specified (in which case the final column order is determined from that). In the
    case where duplicate key values exist, the alignment behaviour is determined by
    the given alignment strategy specified in the `how` parameter (by default this
    is a full outer join, but if your data is suitable you can get a large speedup
    by setting `how="left"` instead).

    Note that this function does not result in a joined frame - you receive the same
    number of frames back that you passed in, but each is now aligned by key and has
    the same number of rows.

    Parameters
    ----------
    frames
        Sequence of DataFrames or LazyFrames.
    on
        One or more columns whose unique values will be used to align the frames.
    select
        Optional post-alignment column select to constrain and/or order
        the columns returned from the newly aligned frames.
    descending
        Sort the alignment column values in descending order; can be a single
        boolean or a list of booleans associated with each column in `on`.
    how
        By default the row alignment values are determined using a full outer join
        strategy across all frames; if you know that the first frame contains all
        required keys, you can set `how="left"` for a large performance increase.

    Examples
    --------
    >>> from datetime import date
    >>> df1 = pl.DataFrame(
    ...     {
    ...         "dt": [date(2022, 9, 1), date(2022, 9, 2), date(2022, 9, 3)],
    ...         "x": [3.5, 4.0, 1.0],
    ...         "y": [10.0, 2.5, 1.5],
    ...     }
    ... )
    >>> df2 = pl.DataFrame(
    ...     {
    ...         "dt": [date(2022, 9, 2), date(2022, 9, 3), date(2022, 9, 1)],
    ...         "x": [8.0, 1.0, 3.5],
    ...         "y": [1.5, 12.0, 5.0],
    ...     }
    ... )
    >>> df3 = pl.DataFrame(
    ...     {
    ...         "dt": [date(2022, 9, 3), date(2022, 9, 2)],
    ...         "x": [2.0, 5.0],
    ...         "y": [2.5, 2.0],
    ...     }
    ... )  # doctest: +IGNORE_RESULT
    >>> pl.Config.set_tbl_formatting("UTF8_FULL")  # doctest: +IGNORE_RESULT
    #
    # df1                              df2                              df3
    # shape: (3, 3)                    shape: (3, 3)                    shape: (2, 3)
    # ┌────────────┬─────┬──────┐      ┌────────────┬─────┬──────┐      ┌────────────┬─────┬─────┐
    # │ dt         ┆ x   ┆ y    │      │ dt         ┆ x   ┆ y    │      │ dt         ┆ x   ┆ y   │
    # │ ---        ┆ --- ┆ ---  │      │ ---        ┆ --- ┆ ---  │      │ ---        ┆ --- ┆ --- │
    # │ date       ┆ f64 ┆ f64  │      │ date       ┆ f64 ┆ f64  │      │ date       ┆ f64 ┆ f64 │
    # ╞════════════╪═════╪══════╡      ╞════════════╪═════╪══════╡      ╞════════════╪═════╪═════╡
    # │ 2022-09-01 ┆ 3.5 ┆ 10.0 │\  ,->│ 2022-09-02 ┆ 8.0 ┆ 1.5  │\  ,->│ 2022-09-03 ┆ 2.0 ┆ 2.5 │
    # ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤ \/   ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤ \/   ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌┤
    # │ 2022-09-02 ┆ 4.0 ┆ 2.5  │_/\,->│ 2022-09-03 ┆ 1.0 ┆ 12.0 │_/`-->│ 2022-09-02 ┆ 5.0 ┆ 2.0 │
    # ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤  /\  ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤      └────────────┴─────┴─────┘
    # │ 2022-09-03 ┆ 1.0 ┆ 1.5  │_/  `>│ 2022-09-01 ┆ 3.5 ┆ 5.0  │-//-
    # └────────────┴─────┴──────┘      └────────────┴─────┴──────┘
    ...

    Align frames by the "dt" column:

    >>> af1, af2, af3 = pl.align_frames(
    ...     df1, df2, df3, on="dt"
    ... )  # doctest: +IGNORE_RESULT
    #
    # df1                              df2                              df3
    # shape: (3, 3)                    shape: (3, 3)                    shape: (3, 3)
    # ┌────────────┬─────┬──────┐      ┌────────────┬─────┬──────┐      ┌────────────┬──────┬──────┐
    # │ dt         ┆ x   ┆ y    │      │ dt         ┆ x   ┆ y    │      │ dt         ┆ x    ┆ y    │
    # │ ---        ┆ --- ┆ ---  │      │ ---        ┆ --- ┆ ---  │      │ ---        ┆ ---  ┆ ---  │
    # │ date       ┆ f64 ┆ f64  │      │ date       ┆ f64 ┆ f64  │      │ date       ┆ f64  ┆ f64  │
    # ╞════════════╪═════╪══════╡      ╞════════════╪═════╪══════╡      ╞════════════╪══════╪══════╡
    # │ 2022-09-01 ┆ 3.5 ┆ 10.0 │----->│ 2022-09-01 ┆ 3.5 ┆ 5.0  │----->│ 2022-09-01 ┆ null ┆ null │
    # ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌┼╌╌╌╌╌╌┤
    # │ 2022-09-02 ┆ 4.0 ┆ 2.5  │----->│ 2022-09-02 ┆ 8.0 ┆ 1.5  │----->│ 2022-09-02 ┆ 5.0  ┆ 2.0  │
    # ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌╌╌╌╌╌╌╌┼╌╌╌╌╌╌┼╌╌╌╌╌╌┤
    # │ 2022-09-03 ┆ 1.0 ┆ 1.5  │----->│ 2022-09-03 ┆ 1.0 ┆ 12.0 │----->│ 2022-09-03 ┆ 2.0  ┆ 2.5  │
    # └────────────┴─────┴──────┘      └────────────┴─────┴──────┘      └────────────┴──────┴──────┘
    ...

    Align frames by "dt" using "left" alignment, but keep only cols "x" and "y":

    >>> af1, af2, af3 = pl.align_frames(
    ...     df1, df2, df3, on="dt", select=["x", "y"], how="left"
    ... )  # doctest: +IGNORE_RESULT
    #
    # af1                 af2                 af3
    # shape: (3, 3)       shape: (3, 3)       shape: (3, 3)
    # ┌─────┬──────┐      ┌─────┬──────┐      ┌──────┬──────┐
    # │ x   ┆ y    │      │ x   ┆ y    │      │ x    ┆ y    │
    # │ --- ┆ ---  │      │ --- ┆ ---  │      │ ---  ┆ ---  │
    # │ f64 ┆ f64  │      │ f64 ┆ f64  │      │ f64  ┆ f64  │
    # ╞═════╪══════╡      ╞═════╪══════╡      ╞══════╪══════╡
    # │ 3.5 ┆ 10.0 │      │ 3.5 ┆ 5.0  │      │ null ┆ null │
    # ├╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌╌┼╌╌╌╌╌╌┤
    # │ 4.0 ┆ 2.5  │      │ 8.0 ┆ 1.5  │      │ 5.0  ┆ 2.0  │
    # ├╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌┼╌╌╌╌╌╌┤      ├╌╌╌╌╌╌┼╌╌╌╌╌╌┤
    # │ 1.0 ┆ 1.5  │      │ 1.0 ┆ 12.0 │      │ 2.0  ┆ 2.5  │
    # └─────┴──────┘      └─────┴──────┘      └──────┴──────┘
    ...

    Now data is aligned, and you can easily calculate the row-wise dot product:

    >>> (af1 * af2 * af3).fill_null(0).select(pl.sum_horizontal("*").alias("dot"))
    shape: (3, 1)
    ┌───────┐
    │ dot   │
    │ ---   │
    │ f64   │
    ╞═══════╡
    │ 0.0   │
    ├╌╌╌╌╌╌╌┤
    │ 167.5 │
    ├╌╌╌╌╌╌╌┤
    │ 47.0  │
    └───────┘
    """  # noqa: W505
    if not frames:
        return []

    if len({type(f) for f in frames}) != 1:
        msg = (
            "input frames must be of a consistent type (all LazyFrame or all DataFrame)"
        )
        raise TypeError(msg)

    eager = isinstance(frames[0], pl.DataFrame)
    on = [on] if (isinstance(on, str) or not isinstance(on, Sequence)) else on
    align_on = [(c.meta.output_name() if isinstance(c, pl.Expr) else c) for c in on]

    # create aligned master frame (this is the most expensive part; afterwards
    # we just subselect out the columns representing the component frames)
    idx_frames = [(idx, frame.lazy()) for idx, frame in enumerate(frames)]
    alignment_frame = _alignment_join(
        *idx_frames, align_on=align_on, how=how, descending=descending
    )

    # select-out aligned components from the master frame
    aligned_cols = set(alignment_frame.collect_schema())
    aligned_frames = []
    for idx, lf in idx_frames:
        sfx = f":{idx}"
        df_cols = [
            F.col(f"{c}{sfx}").alias(c) if f"{c}{sfx}" in aligned_cols else F.col(c)
            for c in lf.collect_schema()
        ]
        f = alignment_frame.select(*df_cols)
        if select is not None:
            f = f.select(select)
        aligned_frames.append(f)

    return F.collect_all(aligned_frames) if eager else aligned_frames  # type: ignore[return-value]
