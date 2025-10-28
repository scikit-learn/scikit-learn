import os
import sys
from typing import Any, Callable


def verbose() -> bool:
    return os.getenv("POLARS_VERBOSE") == "1"


def eprint(*a: Any, **kw: Any) -> None:
    return print(*a, file=sys.stderr, **kw)


def verbose_print_sensitive(create_log_message: Callable[[], str]) -> None:
    if os.getenv("POLARS_VERBOSE_SENSITIVE") == "1":
        # Force the message to be a single line.
        msg = create_log_message().replace("\n", "")
        print(f"[SENSITIVE]: {msg}", file=sys.stderr)
