"""Tools for working with the Polars duration string language."""

from __future__ import annotations

import datetime as dt
import re
from typing import TYPE_CHECKING, Literal, cast, get_args

if TYPE_CHECKING:
    from collections.abc import Container, Mapping

    from typing_extensions import TypeAlias

__all__ = ["IntervalUnit"]

IntervalUnit: TypeAlias = Literal["ns", "us", "ms", "s", "m", "h", "d", "mo", "q", "y"]
"""A Polars duration string interval unit.

- 'ns': nanosecond.
- 'us': microsecond.
- 'ms': millisecond.
- 's': second.
- 'm': minute.
- 'h': hour.
- 'd': day.
- 'mo': month.
- 'q': quarter.
- 'y': year.
"""
TimedeltaKwd: TypeAlias = Literal[
    "days", "hours", "minutes", "seconds", "milliseconds", "microseconds"
]

PATTERN_INTERVAL: re.Pattern[str] = re.compile(
    r"^(?P<multiple>-?\d+)(?P<unit>ns|us|ms|mo|m|s|h|d|q|y)\Z"
)
MONTH_MULTIPLES = frozenset([1, 2, 3, 4, 6, 12])
QUARTER_MULTIPLES = frozenset([1, 2, 4])
UNIT_TO_TIMEDELTA: Mapping[IntervalUnit, TimedeltaKwd] = {
    "d": "days",
    "h": "hours",
    "m": "minutes",
    "s": "seconds",
    "ms": "milliseconds",
    "us": "microseconds",
}


class Interval:
    def __init__(self, multiple: int, unit: IntervalUnit, /) -> None:
        self.multiple: int = multiple
        self.unit: IntervalUnit = unit

    def to_timedelta(
        self, *, unsupported: Container[IntervalUnit] = frozenset(("ns", "mo", "q", "y"))
    ) -> dt.timedelta:
        if self.unit in unsupported:  # pragma: no cover
            msg = f"Creating timedelta with {self.unit} unit is not supported."
            raise NotImplementedError(msg)
        kwd = UNIT_TO_TIMEDELTA[self.unit]
        # error: Keywords must be strings (bad mypy)
        return dt.timedelta(**{kwd: self.multiple})  # type: ignore[misc]

    @classmethod
    def parse(cls, every: str) -> Interval:
        multiple, unit = cls._parse(every)
        if unit == "mo" and multiple not in MONTH_MULTIPLES:
            msg = f"Only the following multiples are supported for 'mo' unit: {MONTH_MULTIPLES}.\nGot: {multiple}."
            raise ValueError(msg)
        if unit == "q" and multiple not in QUARTER_MULTIPLES:
            msg = f"Only the following multiples are supported for 'q' unit: {QUARTER_MULTIPLES}.\nGot: {multiple}."
            raise ValueError(msg)
        if unit == "y" and multiple != 1:
            msg = (
                f"Only multiple 1 is currently supported for 'y' unit.\nGot: {multiple}."
            )
            raise ValueError(msg)
        return cls(multiple, unit)

    @classmethod
    def parse_no_constraints(cls, every: str) -> Interval:
        return cls(*cls._parse(every))

    @staticmethod
    def _parse(every: str) -> tuple[int, IntervalUnit]:
        if match := PATTERN_INTERVAL.match(every):
            multiple = int(match["multiple"])
            unit = cast("IntervalUnit", match["unit"])
            return multiple, unit
        msg = (
            f"Invalid `every` string: {every}. Expected string of kind <number><unit>, "
            f"where 'unit' is one of: {get_args(IntervalUnit)}."
        )
        raise ValueError(msg)
