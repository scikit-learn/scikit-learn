from __future__ import annotations

import inspect
import os
from functools import wraps
from typing import TYPE_CHECKING, Callable, TypeVar

from polars._utils.various import issue_warning
from polars.exceptions import UnstableWarning

if TYPE_CHECKING:
    import sys

    if sys.version_info >= (3, 10):
        from typing import ParamSpec
    else:
        from typing_extensions import ParamSpec

    P = ParamSpec("P")
    T = TypeVar("T")


def issue_unstable_warning(message: str | None = None) -> None:
    """
    Issue a warning for use of unstable functionality.

    The `warn_unstable` setting must be enabled, otherwise no warning is issued.

    Parameters
    ----------
    message
        The message associated with the warning.

    See Also
    --------
    Config.warn_unstable
    """
    warnings_enabled = bool(int(os.environ.get("POLARS_WARN_UNSTABLE", 0)))
    if not warnings_enabled:
        return

    if message is None:
        message = "This functionality is considered unstable."
    message += (
        " It may be changed at any point without it being considered a breaking change."
    )

    issue_warning(message, UnstableWarning)


def unstable() -> Callable[[Callable[P, T]], Callable[P, T]]:
    """Decorator to mark a function as unstable."""

    def decorate(function: Callable[P, T]) -> Callable[P, T]:
        @wraps(function)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> T:
            issue_unstable_warning(f"`{function.__name__}` is considered unstable.")
            return function(*args, **kwargs)

        wrapper.__signature__ = inspect.signature(function)  # type: ignore[attr-defined]
        return wrapper

    return decorate
