from __future__ import annotations

import contextlib
import enum
from datetime import date, datetime, time, timedelta, timezone
from typing import TYPE_CHECKING, Any
from zoneinfo import ZoneInfo

import polars._reexport as pl
from polars._dependencies import (
    _check_for_numpy,
    _check_for_pytz,
    _check_for_torch,
    pytz,
    torch,
)
from polars._dependencies import numpy as np
from polars._utils.wrap import wrap_expr
from polars.datatypes import Date, Datetime, Duration
from polars.datatypes.convert import DataTypeMappings

with contextlib.suppress(ImportError):  # Module not available when building docs
    import polars._plr as plr


if TYPE_CHECKING:
    from polars import Expr
    from polars._typing import PolarsDataType, TimeUnit


def lit(
    value: Any, dtype: PolarsDataType | None = None, *, allow_object: bool = False
) -> Expr:
    """
    Return an expression representing a literal value.

    Parameters
    ----------
    value
        Value that should be used as a `literal`.
    dtype
        The data type of the resulting expression.
        If set to `None` (default), the data type is inferred from the `value` input.
    allow_object
        If type is unknown use an 'object' type.
        By default, we will raise a `ValueException`
        if the type is unknown.

    Notes
    -----
    Expected datatypes:

    - `pl.lit([])` -> empty List<Null>
    - `pl.lit([1, 2, 3])` -> List<i64>
    - `pl.lit(pl.Series([]))`-> empty Series Null
    - `pl.lit(pl.Series([1, 2, 3]))` -> Series Int64
    - `pl.lit(None)` -> Null

    Examples
    --------
    Literal scalar values:

    >>> pl.lit(1)  # doctest: +IGNORE_RESULT
    >>> pl.lit(5.5)  # doctest: +IGNORE_RESULT
    >>> pl.lit(None)  # doctest: +IGNORE_RESULT
    >>> pl.lit("foo_bar")  # doctest: +IGNORE_RESULT
    >>> pl.lit(date(2021, 1, 20))  # doctest: +IGNORE_RESULT
    >>> pl.lit(datetime(2023, 3, 31, 10, 30, 45))  # doctest: +IGNORE_RESULT

    Literal list/Series data (1D):

    >>> pl.lit([1, 2, 3])  # doctest: +SKIP
    >>> pl.lit(pl.Series("x", [1, 2, 3]))  # doctest: +IGNORE_RESULT

    Literal list/Series data (2D):

    >>> pl.lit([[1, 2], [3, 4]])  # doctest: +SKIP
    >>> pl.lit(pl.Series("y", [[1, 2], [3, 4]]))  # doctest: +IGNORE_RESULT
    """
    time_unit: TimeUnit

    if isinstance(value, datetime):
        if dtype == Date:
            return wrap_expr(plr.lit(value.date(), allow_object=False, is_scalar=True))

        # parse time unit
        if dtype is not None and (tu := getattr(dtype, "time_unit", "us")) is not None:
            time_unit = tu  # type: ignore[assignment]
        else:
            time_unit = "us"

        # parse time zone
        dtype_tz = getattr(dtype, "time_zone", None)
        value_tz = value.tzinfo
        if value_tz is None:
            tz = dtype_tz
        else:
            # value has time zone, but dtype does not: keep value time zone
            if dtype_tz is None:
                if isinstance(value_tz, ZoneInfo) or (
                    _check_for_pytz(value_tz)
                    and isinstance(value_tz, pytz.tzinfo.BaseTzInfo)
                    and value_tz.zone is not None
                ):
                    # named timezone
                    tz = str(value_tz)
                else:
                    # fixed offset from UTC (eg: +4:00)
                    value = value.astimezone(timezone.utc)
                    tz = "UTC"

            # dtype and value both have same time zone
            elif str(value_tz) == dtype_tz:
                tz = str(value_tz)

            # given a fixed offset from UTC that matches the dtype tz offset
            elif hasattr(value_tz, "utcoffset") and getattr(
                ZoneInfo(dtype_tz).utcoffset(value), "seconds", 0
            ) == getattr(value_tz.utcoffset(value), "seconds", 1):
                tz = dtype_tz
            else:
                # value has time zone that differs from dtype time zone
                msg = (
                    f"time zone of dtype ({dtype_tz!r}) differs from time zone of "
                    f"value ({value_tz!r})"
                )
                raise TypeError(msg)

        dt_utc = value.replace(tzinfo=timezone.utc)
        expr = wrap_expr(plr.lit(dt_utc, allow_object=False, is_scalar=True)).cast(
            Datetime(time_unit)
        )
        if tz is not None:
            expr = expr.dt.replace_time_zone(
                tz, ambiguous="earliest" if value.fold == 0 else "latest"
            )
        return expr

    elif isinstance(value, timedelta):
        expr = wrap_expr(plr.lit(value, allow_object=False, is_scalar=True))
        if dtype is not None and (tu := getattr(dtype, "time_unit", None)) is not None:
            expr = expr.cast(Duration(tu))
        return expr

    elif isinstance(value, time):
        return wrap_expr(plr.lit(value, allow_object=False, is_scalar=True))

    elif isinstance(value, date):
        if dtype == Datetime:
            time_unit = getattr(dtype, "time_unit", "us") or "us"
            dt_utc = datetime(value.year, value.month, value.day)
            expr = wrap_expr(plr.lit(dt_utc, allow_object=False, is_scalar=True)).cast(
                Datetime(time_unit)
            )
            if (time_zone := getattr(dtype, "time_zone", None)) is not None:
                expr = expr.dt.replace_time_zone(str(time_zone))
            return expr
        else:
            return wrap_expr(plr.lit(value, allow_object=False, is_scalar=True))

    elif isinstance(value, pl.Series):
        value = value._s
        return wrap_expr(plr.lit(value, allow_object, is_scalar=False))

    elif _check_for_numpy(value) and isinstance(value, np.ndarray):
        return lit(pl.Series("literal", value, dtype=dtype))

    elif _check_for_torch(value) and isinstance(value, torch.Tensor):
        return lit(pl.Series("literal", value.numpy(force=False), dtype=dtype))

    elif isinstance(value, (list, tuple)):
        return wrap_expr(
            plr.lit(
                pl.Series("literal", [value], dtype=dtype)._s,
                allow_object,
                is_scalar=True,
            )
        )

    elif isinstance(value, enum.Enum):
        return lit(value.value, dtype=dtype)

    if dtype:
        return wrap_expr(plr.lit(value, allow_object, is_scalar=True)).cast(dtype)

    if _check_for_numpy(value) and isinstance(value, np.generic):
        # note: the item() is a py-native datetime/timedelta when units < 'ns'
        if isinstance(item := value.item(), (date, datetime, timedelta)):
            return lit(item)

        # handle 'ns' units
        if isinstance(item, int) and hasattr(value, "dtype"):
            dtype_name = value.dtype.name
            if dtype_name.startswith("datetime64["):
                time_unit = dtype_name[len("datetime64[") : -1]  # type: ignore[assignment]
                return lit(item).cast(Datetime(time_unit))
            if dtype_name.startswith("timedelta64["):
                time_unit = dtype_name[len("timedelta64[") : -1]  # type: ignore[assignment]
                return lit(item).cast(Duration(time_unit))

        # handle known mappable values
        dtype = DataTypeMappings.NUMPY_KIND_AND_ITEMSIZE_TO_DTYPE.get(
            (value.dtype.kind, value.dtype.itemsize)
        )
        if dtype is not None:
            return lit(value, dtype=dtype)
    else:
        item = value

    return wrap_expr(plr.lit(item, allow_object, is_scalar=True))
