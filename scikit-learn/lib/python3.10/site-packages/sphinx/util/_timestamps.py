from __future__ import annotations

import time


def _format_rfc3339_microseconds(timestamp: int, /) -> str:
    """Return an RFC 3339 formatted string representing the given timestamp.

    :param timestamp: The timestamp to format, in microseconds.
    """
    seconds, fraction = divmod(timestamp, 10**6)
    time_tuple = time.gmtime(seconds)
    return time.strftime('%Y-%m-%d %H:%M:%S', time_tuple) + f'.{fraction // 1_000}'
