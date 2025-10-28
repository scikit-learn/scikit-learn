# Copyright (c) Amber Brown, 2015
# See LICENSE for details.

"""
Responsible for writing the built news fragments to a news file without
affecting existing content.
"""

from __future__ import annotations

import sys

from pathlib import Path
from typing import Any


if sys.version_info < (3, 10):
    # Compatibility shim for newline parameter to write_text, added in 3.10
    def _newline_write_text(path: Path, content: str, **kwargs: Any) -> None:
        with path.open("w", **kwargs) as strm:  # pragma: no branch
            strm.write(content)

else:

    def _newline_write_text(path: Path, content: str, **kwargs: Any) -> None:
        path.write_text(content, **kwargs)


def append_to_newsfile(
    directory: str,
    filename: str,
    start_string: str,
    top_line: str,
    content: str,
    single_file: bool,
) -> None:
    """
    Write *content* to *directory*/*filename* behind *start_string*.

    Double-check *top_line* (i.e. the release header) is not already in the
    file.

    if *single_file* is True, add it to an existing file, otherwise create a
    fresh one.
    """
    news_file = Path(directory) / filename

    header, prev_body = _figure_out_existing_content(
        news_file, start_string, single_file
    )

    if top_line and top_line in prev_body:
        raise ValueError("It seems you've already produced newsfiles for this version?")

    _newline_write_text(
        news_file,
        # If there is no previous body that means we're writing a brand new news file.
        # We don't want extra whitespace at the end of this new file.
        header + (content + prev_body if prev_body else content.rstrip() + "\n"),
        encoding="utf-8",
        # Leave newlines alone. This probably leads to inconsistent newlines,
        # because we've loaded existing content with universal newlines, but that's
        # the original behavior.
        newline="",
    )


def _figure_out_existing_content(
    news_file: Path, start_string: str, single_file: bool
) -> tuple[str, str]:
    """
    Try to read *news_file* and split it into header (everything before
    *start_string*) and the old body (everything after *start_string*).

    If there's no *start_string*, return empty header.

    Empty file and per-release files have neither.
    """
    if not single_file or not news_file.exists():
        # Per-release news files always start empty.
        # Non-existent files have no existing content.
        return "", ""

    content = Path(news_file).read_text(encoding="utf-8")

    t = content.split(start_string, 1)
    if len(t) == 2:
        return f"{t[0].rstrip()}\n\n{start_string}\n", t[1].lstrip()

    return "", content.lstrip()
