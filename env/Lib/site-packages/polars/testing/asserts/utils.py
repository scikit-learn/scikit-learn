from __future__ import annotations

from typing import Any, NoReturn


def raise_assertion_error(
    objects: str, detail: str, left: Any, right: Any, *, cause: Exception | None = None
) -> NoReturn:
    """Raise a detailed assertion error."""
    __tracebackhide__ = True
    msg = f"{objects} are different ({detail})\n[left]:  {left}\n[right]: {right}"
    raise AssertionError(msg) from cause
