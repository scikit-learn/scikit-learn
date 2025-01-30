from __future__ import annotations

from datetime import timedelta

from polars._utils.convert import parse_as_duration_string


def parse_interval_argument(interval: str | timedelta) -> str:
    """Parse the interval argument as a Polars duration string."""
    if isinstance(interval, timedelta):
        return parse_as_duration_string(interval)

    if " " in interval:
        interval = interval.replace(" ", "")
    return interval.lower()
