"""Docutils utility functions for parsing text."""

from __future__ import annotations

import contextlib
from typing import TYPE_CHECKING

from docutils.nodes import Element, Node
from docutils.statemachine import StringList, string2lines

if TYPE_CHECKING:
    from collections.abc import Iterator

    from docutils.parsers.rst.states import RSTState


def nested_parse_to_nodes(
    state: RSTState,
    text: str | StringList,
    *,
    source: str = '<generated text>',
    offset: int = 0,
    allow_section_headings: bool = True,
    keep_title_context: bool = False,
) -> list[Node]:  # Element | nodes.Text
    """Parse *text* into nodes.

    :param state:
        The state machine state. Must be a subclass of ``RSTState``.
    :param text:
        Text, in string form. ``StringList`` is also accepted.
    :param source:
        The text's source, used when creating a new ``StringList``.
    :param offset:
        The offset of the content.
    :param allow_section_headings:
        Are titles (sections) allowed in *text*?
        Note that this option bypasses Docutils' usual checks on
        doctree structure, and misuse of this option can lead to
        an incoherent doctree. In Docutils, section nodes should
        only be children of ``Structural`` nodes, which includes
        ``document``, ``section``, and ``sidebar`` nodes.
    :param keep_title_context:
        If this is False (the default), then *content* is parsed as if it were
        an independent document, meaning that title decorations (e.g. underlines)
        do not need to match the surrounding document.
        This is useful when the parsed content comes from
        a completely different context, such as docstrings.
        If this is True, then title underlines must match those in
        the surrounding document, otherwise the behaviour is undefined.

    .. versionadded:: 7.4
    """
    document = state.document
    content = _text_to_string_list(
        text, source=source, tab_width=document.settings.tab_width
    )
    node = Element()  # Anonymous container for parsing
    node.document = document

    if keep_title_context:
        state.nested_parse(content, offset, node, match_titles=allow_section_headings)
    else:
        with _fresh_title_style_context(state):
            state.nested_parse(
                content, offset, node, match_titles=allow_section_headings
            )
    return node.children


@contextlib.contextmanager
def _fresh_title_style_context(state: RSTState) -> Iterator[None]:
    # hack around title style bookkeeping
    memo = state.memo
    surrounding_title_styles: list[str | tuple[str, str]] = memo.title_styles
    surrounding_section_level: int = memo.section_level
    # clear current title styles
    memo.title_styles = []
    memo.section_level = 0
    try:
        yield
    finally:
        # reset title styles
        memo.title_styles = surrounding_title_styles
        memo.section_level = surrounding_section_level


def _text_to_string_list(
    text: str | StringList, /, *, source: str, tab_width: int
) -> StringList:
    # Doesn't really belong in this module, but avoids circular imports.
    if isinstance(text, StringList):
        return text
    content = string2lines(text, tab_width, convert_whitespace=True)
    return StringList(content, source=source)
