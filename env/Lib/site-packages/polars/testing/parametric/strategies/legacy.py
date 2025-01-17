from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Any

import hypothesis.strategies as st
from hypothesis.errors import InvalidArgument

from polars._utils.deprecation import deprecate_function
from polars.datatypes import is_polars_dtype
from polars.testing.parametric.strategies.core import _COL_LIMIT, column
from polars.testing.parametric.strategies.data import lists
from polars.testing.parametric.strategies.dtype import _instantiate_dtype, dtypes

if TYPE_CHECKING:
    from hypothesis.strategies import SearchStrategy

    from polars._typing import OneOrMoreDataTypes, PolarsDataType


@deprecate_function(
    "Use `column` instead in conjunction with a list comprehension.", version="0.20.26"
)
def columns(
    cols: int | Sequence[str] | None = None,
    *,
    dtype: OneOrMoreDataTypes | None = None,
    min_cols: int = 0,
    max_cols: int = _COL_LIMIT,
    unique: bool = False,
) -> list[column]:
    """
    Define multiple columns for use with the @dataframes strategy.

    .. deprecated:: 0.20.26
        Use :class:`column` instead in conjunction with a list comprehension.

    .. warning::
        This functionality is currently considered **unstable**. It may be
        changed at any point without it being considered a breaking change.

    Generate a fixed sequence of `column` objects suitable for passing to the
    @dataframes strategy, or using standalone (note that this function is not itself
    a strategy).

    Notes
    -----
    Additional control is available by creating a sequence of columns explicitly,
    using the `column` class (an especially useful option is to override the default
    data-generating strategy for a given col/dtype).

    Parameters
    ----------
    cols : {int, [str]}, optional
        integer number of cols to create, or explicit list of column names. if
        omitted a random number of columns (between mincol and max_cols) are
        created.
    dtype : PolarsDataType, optional
        a single dtype for all cols, or list of dtypes (the same length as `cols`).
        if omitted, each generated column is assigned a random dtype.
    min_cols : int, optional
        if not passing an exact size, can set a minimum here (defaults to 0).
    max_cols : int, optional
        if not passing an exact size, can set a maximum value here (defaults to
        MAX_COLS).
    unique : bool, optional
        indicate if the values generated for these columns should be unique
        (per-column).

    Examples
    --------
    >>> from polars.testing.parametric import columns, dataframes
    >>> from hypothesis import given
    >>> @given(dataframes(columns(["x", "y", "z"], unique=True)))  # doctest: +SKIP
    ... def test_unique_xyz(df: pl.DataFrame) -> None:
    ...     assert_something(df)
    """
    # create/assign named columns
    if cols is None:
        cols = st.integers(min_value=min_cols, max_value=max_cols).example()
    if isinstance(cols, int):
        names: Sequence[str] = [f"col{n}" for n in range(cols)]
    else:
        names = cols
    n_cols = len(names)

    if dtype is None:
        dtypes: Sequence[PolarsDataType | None] = [None] * n_cols
    elif is_polars_dtype(dtype):
        dtypes = [dtype] * n_cols
    elif isinstance(dtype, Sequence):
        if (n_dtypes := len(dtype)) != n_cols:
            msg = f"given {n_dtypes} dtypes for {n_cols} names"
            raise InvalidArgument(msg)
        dtypes = dtype
    else:
        msg = f"{dtype!r} is not a valid polars datatype"
        raise InvalidArgument(msg)

    # init list of named/typed columns
    return [column(name=nm, dtype=tp, unique=unique) for nm, tp in zip(names, dtypes)]


@deprecate_function("Use `lists` instead.", version="0.20.26")
def create_list_strategy(
    inner_dtype: PolarsDataType | None = None,
    *,
    select_from: Sequence[Any] | None = None,
    size: int | None = None,
    min_size: int = 0,
    max_size: int | None = None,
    unique: bool = False,
) -> SearchStrategy[list[Any]]:
    """
    Create a strategy for generating Polars :class:`List` data.

    .. deprecated:: 0.20.26
        Use :func:`lists` instead.

    Parameters
    ----------
    inner_dtype : PolarsDataType
        type of the inner list elements (can also be another List).
    select_from : list, optional
        randomly select the innermost values from this list (otherwise
        the default strategy associated with the innermost dtype is used).
    size : int, optional
        if set, generated lists will be of exactly this size (and
        ignore the min_size/max_size params).
    min_size : int, optional
        set the minimum size of the generated lists (default: 0 if unset).
    max_size : int, optional
        set the maximum size of the generated lists (default: 3 if
        min_size is unset or zero, otherwise 2x min_size).
    unique : bool, optional
        ensure that the generated lists contain unique values.

    Examples
    --------
    Create a strategy that generates a list of i32 values:

    >>> from polars.testing.parametric import create_list_strategy
    >>> lst = create_list_strategy(inner_dtype=pl.Int32)  # doctest: +SKIP
    >>> lst.example()  # doctest: +SKIP
    [-11330, 24030, 116]
    """
    if size is not None:
        min_size = max_size = size

    if inner_dtype is None:
        inner_dtype = dtypes().example()
    else:
        inner_dtype = _instantiate_dtype(inner_dtype).example()

    return lists(
        inner_dtype,
        select_from=select_from,
        min_size=min_size,
        max_size=max_size,
        unique=unique,
    )
