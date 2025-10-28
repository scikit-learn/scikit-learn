from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

if TYPE_CHECKING:
    from collections.abc import Iterator


@contextmanager
def PortableTemporaryFile(
    mode: str = "w+b",
    *,
    buffering: int = -1,
    encoding: str | None = None,
    newline: str | None = None,
    suffix: str | None = None,
    prefix: str | None = None,
    dir: str | Path | None = None,
    delete: bool = True,
    errors: str | None = None,
) -> Iterator[Any]:
    """
    Slightly more resilient version of the standard `NamedTemporaryFile`.

    Plays better with Windows when using the 'delete' option.
    """
    from tempfile import NamedTemporaryFile

    params = cast(
        Any,
        {
            "mode": mode,
            "buffering": buffering,
            "encoding": encoding,
            "newline": newline,
            "suffix": suffix,
            "prefix": prefix,
            "dir": dir,
            "delete": False,
            "errors": errors,
        },
    )

    with NamedTemporaryFile(**params) as tmp:
        try:
            yield tmp
        finally:
            tmp.close()
            if delete:
                Path(tmp.name).unlink(missing_ok=True)
