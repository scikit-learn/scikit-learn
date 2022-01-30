"""
    sphinx.util.matching
    ~~~~~~~~~~~~~~~~~~~~

    Pattern-matching utility functions for Sphinx.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import re
from typing import Callable, Dict, Iterable, List, Match, Optional, Pattern

from sphinx.util.osutil import canon_path


def _translate_pattern(pat: str) -> str:
    """Translate a shell-style glob pattern to a regular expression.

    Adapted from the fnmatch module, but enhanced so that single stars don't
    match slashes.
    """
    i, n = 0, len(pat)
    res = ''
    while i < n:
        c = pat[i]
        i += 1
        if c == '*':
            if i < n and pat[i] == '*':
                # double star matches slashes too
                i += 1
                res = res + '.*'
            else:
                # single star doesn't match slashes
                res = res + '[^/]*'
        elif c == '?':
            # question mark doesn't match slashes too
            res = res + '[^/]'
        elif c == '[':
            j = i
            if j < n and pat[j] == '!':
                j += 1
            if j < n and pat[j] == ']':
                j += 1
            while j < n and pat[j] != ']':
                j += 1
            if j >= n:
                res = res + '\\['
            else:
                stuff = pat[i:j].replace('\\', '\\\\')
                i = j + 1
                if stuff[0] == '!':
                    # negative pattern mustn't match slashes too
                    stuff = '^/' + stuff[1:]
                elif stuff[0] == '^':
                    stuff = '\\' + stuff
                res = '%s[%s]' % (res, stuff)
        else:
            res += re.escape(c)
    return res + '$'


def compile_matchers(patterns: List[str]) -> List[Callable[[str], Optional[Match[str]]]]:
    return [re.compile(_translate_pattern(pat)).match for pat in patterns]


class Matcher:
    """A pattern matcher for Multiple shell-style glob patterns.

    Note: this modifies the patterns to work with copy_asset().
          For example, "**/index.rst" matches with "index.rst"
    """

    def __init__(self, patterns: List[str]) -> None:
        expanded = [pat[3:] for pat in patterns if pat.startswith('**/')]
        self.patterns = compile_matchers(patterns + expanded)

    def __call__(self, string: str) -> bool:
        return self.match(string)

    def match(self, string: str) -> bool:
        string = canon_path(string)
        return any(pat(string) for pat in self.patterns)


DOTFILES = Matcher(['**/.*'])


_pat_cache: Dict[str, Pattern] = {}


def patmatch(name: str, pat: str) -> Optional[Match[str]]:
    """Return if name matches the regular expression (pattern)
    ``pat```. Adapted from fnmatch module."""
    if pat not in _pat_cache:
        _pat_cache[pat] = re.compile(_translate_pattern(pat))
    return _pat_cache[pat].match(name)


def patfilter(names: Iterable[str], pat: str) -> List[str]:
    """Return the subset of the list ``names`` that match
    the regular expression (pattern) ``pat``.

    Adapted from fnmatch module.
    """
    if pat not in _pat_cache:
        _pat_cache[pat] = re.compile(_translate_pattern(pat))
    match = _pat_cache[pat].match
    return list(filter(match, names))
