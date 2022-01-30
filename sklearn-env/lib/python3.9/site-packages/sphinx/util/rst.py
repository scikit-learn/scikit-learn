"""
    sphinx.util.rst
    ~~~~~~~~~~~~~~~

    reST helper functions.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re
from collections import defaultdict
from contextlib import contextmanager
from typing import Dict, Generator
from unicodedata import east_asian_width

from docutils.parsers.rst import roles
from docutils.parsers.rst.languages import en as english
from docutils.statemachine import StringList
from docutils.utils import Reporter
from jinja2 import Environment

from sphinx.locale import __
from sphinx.util import docutils, logging

try:
    from jinja2.utils import pass_environment
except ImportError:
    from jinja2 import environmentfilter as pass_environment


logger = logging.getLogger(__name__)

docinfo_re = re.compile(':\\w+:.*?')
symbols_re = re.compile(r'([!-\-/:-@\[-`{-~])')  # symbols without dot(0x2e)
SECTIONING_CHARS = ['=', '-', '~']

# width of characters
WIDECHARS: Dict[str, str] = defaultdict(lambda: "WF")  # WF: Wide + Full-width
WIDECHARS["ja"] = "WFA"  # In Japanese, Ambiguous characters also have double width


def escape(text: str) -> str:
    text = symbols_re.sub(r'\\\1', text)
    text = re.sub(r'^\.', r'\.', text)  # escape a dot at top
    return text


def textwidth(text: str, widechars: str = 'WF') -> int:
    """Get width of text."""
    def charwidth(char: str, widechars: str) -> int:
        if east_asian_width(char) in widechars:
            return 2
        else:
            return 1

    return sum(charwidth(c, widechars) for c in text)


@pass_environment
def heading(env: Environment, text: str, level: int = 1) -> str:
    """Create a heading for *level*."""
    assert level <= 3
    width = textwidth(text, WIDECHARS[env.language])
    sectioning_char = SECTIONING_CHARS[level - 1]
    return '%s\n%s' % (text, sectioning_char * width)


@contextmanager
def default_role(docname: str, name: str) -> Generator[None, None, None]:
    if name:
        dummy_reporter = Reporter('', 4, 4)
        role_fn, _ = roles.role(name, english, 0, dummy_reporter)
        if role_fn:
            docutils.register_role('', role_fn)
        else:
            logger.warning(__('default role %s not found'), name, location=docname)

    yield

    docutils.unregister_role('')


def prepend_prolog(content: StringList, prolog: str) -> None:
    """Prepend a string to content body as prolog."""
    if prolog:
        pos = 0
        for line in content:
            if docinfo_re.match(line):
                pos += 1
            else:
                break

        if pos > 0:
            # insert a blank line after docinfo
            content.insert(pos, '', '<generated>', 0)
            pos += 1

        # insert prolog (after docinfo if exists)
        for lineno, line in enumerate(prolog.splitlines()):
            content.insert(pos + lineno, line, '<rst_prolog>', lineno)

        content.insert(pos + lineno + 1, '', '<generated>', 0)


def append_epilog(content: StringList, epilog: str) -> None:
    """Append a string to content body as epilog."""
    if epilog:
        if 0 < len(content):
            source, lineno = content.info(-1)
        else:
            source = '<generated>'
            lineno = 0
        content.append('', source, lineno + 1)
        for lineno, line in enumerate(epilog.splitlines()):
            content.append(line, '<rst_epilog>', lineno)
