from __future__ import annotations

import contextlib
from datetime import date
from typing import TYPE_CHECKING

from polars._utils.deprecation import deprecate_nonkeyword_arguments
from polars._utils.parse import parse_into_expression
from polars._utils.unstable import unstable
from polars._utils.wrap import wrap_expr

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr

if TYPE_CHECKING:
    from collections.abc import Iterable

    from polars import Expr
    from polars._typing import IntoExprColumn


@unstable()
@deprecate_nonkeyword_arguments(allowed_args=["start", "end"], version="1.27.0")
def business_day_count(
    start: date | IntoExprColumn,
    end: date | IntoExprColumn,
    week_mask: Iterable[bool] = (True, True, True, True, True, False, False),
    holidays: Iterable[date] = (),
) -> Expr:
    """
    Count the number of business days between `start` and `end` (not including `end`).

    .. warning::
        This functionality is considered **unstable**. It may be changed
        at any point without it being considered a breaking change.

    .. versionchanged:: 1.27.0
        Parameters after `start` and `end` should now be passed as keyword arguments.

    Parameters
    ----------
    start
        Start dates.
    end
        End dates.
    week_mask
        Which days of the week to count. The default is Monday to Friday.
        If you wanted to count only Monday to Thursday, you would pass
        `(True, True, True, True, False, False, False)`.
    holidays
        Holidays to exclude from the count. The Python package
        `python-holidays <https://github.com/vacanza/python-holidays>`_
        may come in handy here. You can install it with ``pip install holidays``,
        and then, to get all Dutch holidays for years 2020-2024:

        .. code-block:: python

            import holidays

            my_holidays = holidays.country_holidays("NL", years=range(2020, 2025))

        and pass `holidays=my_holidays` when you call `business_day_count`.

    Returns
    -------
    Expr

    Examples
    --------
    >>> from datetime import date
    >>> df = pl.DataFrame(
    ...     {
    ...         "start": [date(2020, 1, 1), date(2020, 1, 2)],
    ...         "end": [date(2020, 1, 2), date(2020, 1, 10)],
    ...     }
    ... )
    >>> df.with_columns(
    ...     business_day_count=pl.business_day_count("start", "end"),
    ... )
    shape: (2, 3)
    ┌────────────┬────────────┬────────────────────┐
    │ start      ┆ end        ┆ business_day_count │
    │ ---        ┆ ---        ┆ ---                │
    │ date       ┆ date       ┆ i32                │
    ╞════════════╪════════════╪════════════════════╡
    │ 2020-01-01 ┆ 2020-01-02 ┆ 1                  │
    │ 2020-01-02 ┆ 2020-01-10 ┆ 6                  │
    └────────────┴────────────┴────────────────────┘

    Note how the business day count is 6 (as opposed a regular day count of 8)
    due to the weekend (2020-01-04 - 2020-01-05) not being counted.

    You can pass a custom weekend - for example, if you only take Sunday off:

    >>> week_mask = (True, True, True, True, True, True, False)
    >>> df.with_columns(
    ...     business_day_count=pl.business_day_count(
    ...         "start", "end", week_mask=week_mask
    ...     ),
    ... )
    shape: (2, 3)
    ┌────────────┬────────────┬────────────────────┐
    │ start      ┆ end        ┆ business_day_count │
    │ ---        ┆ ---        ┆ ---                │
    │ date       ┆ date       ┆ i32                │
    ╞════════════╪════════════╪════════════════════╡
    │ 2020-01-01 ┆ 2020-01-02 ┆ 1                  │
    │ 2020-01-02 ┆ 2020-01-10 ┆ 7                  │
    └────────────┴────────────┴────────────────────┘

    You can also pass a list of holidays to exclude from the count:

    >>> from datetime import date
    >>> holidays = [date(2020, 1, 1), date(2020, 1, 2)]
    >>> df.with_columns(
    ...     business_day_count=pl.business_day_count("start", "end", holidays=holidays)
    ... )
    shape: (2, 3)
    ┌────────────┬────────────┬────────────────────┐
    │ start      ┆ end        ┆ business_day_count │
    │ ---        ┆ ---        ┆ ---                │
    │ date       ┆ date       ┆ i32                │
    ╞════════════╪════════════╪════════════════════╡
    │ 2020-01-01 ┆ 2020-01-02 ┆ 0                  │
    │ 2020-01-02 ┆ 2020-01-10 ┆ 5                  │
    └────────────┴────────────┴────────────────────┘
    """
    start_pyexpr = parse_into_expression(start)
    end_pyexpr = parse_into_expression(end)
    unix_epoch = date(1970, 1, 1)
    return wrap_expr(
        plr.business_day_count(
            start_pyexpr,
            end_pyexpr,
            list(week_mask),
            [(holiday - unix_epoch).days for holiday in holidays],
        )
    )
