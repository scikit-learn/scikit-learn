from contextlib import contextmanager
from typing import Optional, Tuple, Iterator

# These are global mutable state. Don't add anything here unless there's a very
# good reason.

# Value varies by file being processed
strict_optional = False
find_occurrences: Optional[Tuple[str, str]] = None


@contextmanager
def strict_optional_set(value: bool) -> Iterator[None]:
    global strict_optional
    saved = strict_optional
    strict_optional = value
    try:
        yield
    finally:
        strict_optional = saved
