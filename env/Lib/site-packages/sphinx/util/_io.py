from __future__ import annotations

from typing import TYPE_CHECKING

from sphinx.util.console import strip_escape_sequences

if TYPE_CHECKING:
    from typing import Protocol

    class SupportsWrite(Protocol):
        def write(self, text: str, /) -> int | None: ...  # NoQA: E704


class TeeStripANSI:
    """File-like object writing to two streams."""

    def __init__(
        self,
        stream_term: SupportsWrite,
        stream_file: SupportsWrite,
    ) -> None:
        self.stream_term = stream_term
        self.stream_file = stream_file

    def write(self, text: str, /) -> None:
        self.stream_term.write(text)
        self.stream_file.write(strip_escape_sequences(text))

    def flush(self) -> None:
        if hasattr(self.stream_term, 'flush'):
            self.stream_term.flush()
        if hasattr(self.stream_file, 'flush'):
            self.stream_file.flush()
