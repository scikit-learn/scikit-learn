from __future__ import annotations

from datetime import datetime, time, timedelta, timezone
from decimal import Context
from functools import lru_cache
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    NoReturn,
    overload,
)
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError

from polars._utils.constants import (
    EPOCH,
    EPOCH_DATE,
    EPOCH_UTC,
    MS_PER_SECOND,
    NS_PER_SECOND,
    SECONDS_PER_DAY,
    SECONDS_PER_HOUR,
    US_PER_SECOND,
)

if TYPE_CHECKING:
    from datetime import date, tzinfo
    from decimal import Decimal

    from polars._typing import TimeUnit


@overload
def parse_as_duration_string(td: None) -> None: ...


@overload
def parse_as_duration_string(td: timedelta | str) -> str: ...


def parse_as_duration_string(td: timedelta | str | None) -> str | None:
    """Parse duration input as a Polars duration string."""
    if td is None or isinstance(td, str):
        return td
    return _timedelta_to_duration_string(td)


def _timedelta_to_duration_string(td: timedelta) -> str:
    """Convert a Python timedelta object to a Polars duration string."""
    # Positive duration
    if td.days >= 0:
        d = f"{td.days}d" if td.days != 0 else ""
        s = f"{td.seconds}s" if td.seconds != 0 else ""
        us = f"{td.microseconds}us" if td.microseconds != 0 else ""
    # Negative, whole days
    elif td.seconds == 0 and td.microseconds == 0:
        return f"{td.days}d"
    # Negative, other
    else:
        corrected_d = td.days + 1
        corrected_seconds = SECONDS_PER_DAY - (td.seconds + (td.microseconds > 0))
        d = f"{corrected_d}d" if corrected_d != 0 else "-"
        s = f"{corrected_seconds}s" if corrected_seconds != 0 else ""
        us = f"{10**6 - td.microseconds}us" if td.microseconds != 0 else ""

    return f"{d}{s}{us}"


def negate_duration_string(duration: str) -> str:
    """Negate a Polars duration string."""
    if duration.startswith("-"):
        return duration[1:]
    else:
        return f"-{duration}"


def date_to_int(d: date) -> int:
    """Convert a Python time object to an integer."""
    return (d - EPOCH_DATE).days


def time_to_int(t: time) -> int:
    """Convert a Python time object to an integer."""
    t = t.replace(tzinfo=timezone.utc)
    seconds = t.hour * SECONDS_PER_HOUR + t.minute * 60 + t.second
    microseconds = t.microsecond
    return seconds * NS_PER_SECOND + microseconds * 1_000


def datetime_to_int(dt: datetime, time_unit: TimeUnit) -> int:
    """Convert a Python datetime object to an integer."""
    # Make sure to use UTC rather than system time zone
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)

    td = dt - EPOCH_UTC
    seconds = td.days * SECONDS_PER_DAY + td.seconds
    microseconds = dt.microsecond

    if time_unit == "us":
        return seconds * US_PER_SECOND + microseconds
    elif time_unit == "ns":
        return seconds * NS_PER_SECOND + microseconds * 1_000
    elif time_unit == "ms":
        return seconds * MS_PER_SECOND + microseconds // 1_000
    else:
        _raise_invalid_time_unit(time_unit)


def timedelta_to_int(td: timedelta, time_unit: TimeUnit) -> int:
    """Convert a Python timedelta object to an integer."""
    seconds = td.days * SECONDS_PER_DAY + td.seconds
    microseconds = td.microseconds

    if time_unit == "us":
        return seconds * US_PER_SECOND + microseconds
    elif time_unit == "ns":
        return seconds * NS_PER_SECOND + microseconds * 1_000
    elif time_unit == "ms":
        return seconds * MS_PER_SECOND + microseconds // 1_000
    else:
        _raise_invalid_time_unit(time_unit)


@lru_cache(256)
def to_py_date(value: int | float) -> date:
    """Convert an integer or float to a Python date object."""
    return EPOCH_DATE + timedelta(days=value)


def to_py_time(value: int) -> time:
    """Convert an integer to a Python time object."""
    # Fast path for 00:00
    if value == 0:
        return time()

    seconds, nanoseconds = divmod(value, NS_PER_SECOND)
    minutes, seconds = divmod(seconds, 60)
    hours, minutes = divmod(minutes, 60)
    return time(
        hour=hours, minute=minutes, second=seconds, microsecond=nanoseconds // 1_000
    )


def to_py_datetime(
    value: int | float,
    time_unit: TimeUnit,
    time_zone: str | None = None,
) -> datetime:
    """Convert an integer or float to a Python datetime object."""
    if time_unit == "us":
        td = timedelta(microseconds=value)
    elif time_unit == "ns":
        td = timedelta(microseconds=value // 1_000)
    elif time_unit == "ms":
        td = timedelta(milliseconds=value)
    else:
        _raise_invalid_time_unit(time_unit)

    if time_zone is None:
        return EPOCH + td
    else:
        dt = EPOCH_UTC + td
        return _localize_datetime(dt, time_zone)


def _localize_datetime(dt: datetime, time_zone: str) -> datetime:
    # zone info installation should already be checked
    tz: ZoneInfo | tzinfo
    try:
        tz = ZoneInfo(time_zone)
    except ZoneInfoNotFoundError:
        # try fixed offset, which is not supported by ZoneInfo
        tz = _parse_fixed_tz_offset(time_zone)

    return dt.astimezone(tz)


# cache here as we have a single tz per column
# and this function will be called on every conversion
@lru_cache(16)
def _parse_fixed_tz_offset(offset: str) -> tzinfo:
    try:
        # use fromisoformat to parse the offset
        dt_offset = datetime.fromisoformat("2000-01-01T00:00:00" + offset)

        # alternatively, we parse the offset ourselves extracting hours and
        # minutes, then we can construct:
        # tzinfo=timezone(timedelta(hours=..., minutes=...))
    except ValueError:
        msg = f"unexpected time zone offset: {offset!r}"
        raise ValueError(msg) from None

    return dt_offset.tzinfo  # type: ignore[return-value]


def to_py_timedelta(value: int | float, time_unit: TimeUnit) -> timedelta:
    """Convert an integer or float to a Python timedelta object."""
    if time_unit == "us":
        return timedelta(microseconds=value)
    elif time_unit == "ns":
        return timedelta(microseconds=value // 1_000)
    elif time_unit == "ms":
        return timedelta(milliseconds=value)
    else:
        _raise_invalid_time_unit(time_unit)


def to_py_decimal(prec: int, value: str) -> Decimal:
    """Convert decimal components to a Python Decimal object."""
    return _create_decimal_with_prec(prec)(value)


@lru_cache(None)
def _create_decimal_with_prec(
    precision: int,
) -> Callable[[str], Decimal]:
    # pre-cache contexts so we don't have to spend time on recreating them every time
    return Context(prec=precision).create_decimal


def _raise_invalid_time_unit(time_unit: Any) -> NoReturn:
    msg = f"`time_unit` must be one of {{'ms', 'us', 'ns'}}, got {time_unit!r}"
    raise ValueError(msg)
