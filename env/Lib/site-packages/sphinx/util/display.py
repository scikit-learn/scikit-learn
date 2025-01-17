from __future__ import annotations

import functools

from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.console import bold, color_terminal

if False:
    from collections.abc import Callable, Iterable, Iterator
    from types import TracebackType
    from typing import Any, TypeVar

    from typing_extensions import ParamSpec

    T = TypeVar('T')
    P = ParamSpec('P')
    R = TypeVar('R')

logger = logging.getLogger(__name__)


def display_chunk(chunk: Any) -> str:
    if isinstance(chunk, list | tuple):
        if len(chunk) == 1:
            return str(chunk[0])
        return f'{chunk[0]} .. {chunk[-1]}'
    return str(chunk)


def status_iterator(
    iterable: Iterable[T],
    summary: str,
    color: str = 'darkgreen',
    length: int = 0,
    verbosity: int = 0,
    stringify_func: Callable[[Any], str] = display_chunk,
) -> Iterator[T]:
    # printing on a single line requires ANSI control sequences
    single_line = verbosity < 1 and color_terminal()
    bold_summary = bold(summary)
    if length == 0:
        logger.info(bold_summary, nonl=True)
        for item in iterable:
            logger.info(stringify_func(item) + ' ', nonl=True, color=color)
            yield item
    else:
        for i, item in enumerate(iterable, start=1):
            if single_line:
                # clear the entire line ('Erase in Line')
                logger.info('\x1b[2K', nonl=True)
            logger.info(f'{bold_summary}[{i / length: >4.0%}] ', nonl=True)  # NoQA: G004
            # Emit the string representation of ``item``
            logger.info(stringify_func(item), nonl=True, color=color)
            # If in single-line mode, emit a carriage return to move the cursor
            # to the start of the line.
            # If not, emit a newline to move the cursor to the next line.
            logger.info('\r' * single_line, nonl=single_line)
            yield item
    logger.info('')


class SkipProgressMessage(Exception):
    pass


class progress_message:
    def __init__(self, message: str, *, nonl: bool = True) -> None:
        self.message = message
        self.nonl = nonl

    def __enter__(self) -> None:
        logger.info(bold(self.message + '... '), nonl=self.nonl)

    def __exit__(
        self,
        typ: type[BaseException] | None,
        val: BaseException | None,
        tb: TracebackType | None,
    ) -> bool:
        prefix = '' if self.nonl else bold(self.message + ': ')
        if isinstance(val, SkipProgressMessage):
            logger.info(prefix + __('skipped'))
            if val.args:
                logger.info(*val.args)
            return True
        elif val:
            logger.info(prefix + __('failed'))
        else:
            logger.info(prefix + __('done'))

        return False

    def __call__(self, f: Callable[P, R]) -> Callable[P, R]:
        @functools.wraps(f)
        def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:  # type: ignore[return]
            with self:
                return f(*args, **kwargs)

        return wrapper
