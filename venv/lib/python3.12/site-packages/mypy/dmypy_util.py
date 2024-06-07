"""Shared code between dmypy.py and dmypy_server.py.

This should be pretty lightweight and not depend on other mypy code (other than ipc).
"""

from __future__ import annotations

import io
import json
from types import TracebackType
from typing import Any, Final, Iterable, Iterator, TextIO

from mypy.ipc import IPCBase

DEFAULT_STATUS_FILE: Final = ".dmypy.json"


def receive(connection: IPCBase) -> Any:
    """Receive single JSON data frame from a connection.

    Raise OSError if the data received is not valid JSON or if it is
    not a dict.
    """
    bdata = connection.read()
    if not bdata:
        raise OSError("No data received")
    try:
        data = json.loads(bdata)
    except Exception as e:
        raise OSError("Data received is not valid JSON") from e
    if not isinstance(data, dict):
        raise OSError(f"Data received is not a dict ({type(data)})")
    return data


def send(connection: IPCBase, data: Any) -> None:
    """Send data to a connection encoded and framed.

    The data must be JSON-serializable. We assume that a single send call is a
    single frame to be sent on the connect.
    """
    connection.write(json.dumps(data))


class WriteToConn(TextIO):
    """Helper class to write to a connection instead of standard output."""

    def __init__(self, server: IPCBase, output_key: str, isatty: bool) -> None:
        self.server = server
        self.output_key = output_key
        self._isatty = isatty

    def __enter__(self) -> TextIO:
        return self

    def __exit__(
        self,
        t: type[BaseException] | None,
        value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        pass

    def __iter__(self) -> Iterator[str]:
        raise io.UnsupportedOperation

    def __next__(self) -> str:
        raise io.UnsupportedOperation

    def close(self) -> None:
        pass

    def fileno(self) -> int:
        raise OSError

    def flush(self) -> None:
        pass

    def isatty(self) -> bool:
        return self._isatty

    def read(self, n: int = 0) -> str:
        raise io.UnsupportedOperation

    def readable(self) -> bool:
        return False

    def readline(self, limit: int = 0) -> str:
        raise io.UnsupportedOperation

    def readlines(self, hint: int = 0) -> list[str]:
        raise io.UnsupportedOperation

    def seek(self, offset: int, whence: int = 0) -> int:
        raise io.UnsupportedOperation

    def seekable(self) -> bool:
        return False

    def tell(self) -> int:
        raise io.UnsupportedOperation

    def truncate(self, size: int | None = 0) -> int:
        raise io.UnsupportedOperation

    def write(self, output: str) -> int:
        resp: dict[str, Any] = {}
        resp[self.output_key] = output
        send(self.server, resp)
        return len(output)

    def writable(self) -> bool:
        return True

    def writelines(self, lines: Iterable[str]) -> None:
        for s in lines:
            self.write(s)
