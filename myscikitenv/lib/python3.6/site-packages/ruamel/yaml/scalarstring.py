# coding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from ruamel.yaml.compat import text_type
from ruamel.yaml.anchor import Anchor

if False:  # MYPY
    from typing import Text, Any, Dict, List  # NOQA

__all__ = [
    'ScalarString',
    'LiteralScalarString',
    'FoldedScalarString',
    'SingleQuotedScalarString',
    'DoubleQuotedScalarString',
    'PlainScalarString',
    # PreservedScalarString is the old name, as it was the first to be preserved on rt,
    # use LiteralScalarString instead
    'PreservedScalarString',
]


class ScalarString(text_type):
    __slots__ = Anchor.attrib

    def __new__(cls, *args, **kw):
        # type: (Any, Any) -> Any
        anchor = kw.pop('anchor', None)  # type: ignore
        ret_val = text_type.__new__(cls, *args, **kw)  # type: ignore
        if anchor is not None:
            ret_val.yaml_set_anchor(anchor, always_dump=True)
        return ret_val

    def replace(self, old, new, maxreplace=-1):
        # type: (Any, Any, int) -> Any
        return type(self)((text_type.replace(self, old, new, maxreplace)))

    @property
    def anchor(self):
        # type: () -> Any
        if not hasattr(self, Anchor.attrib):
            setattr(self, Anchor.attrib, Anchor())
        return getattr(self, Anchor.attrib)

    def yaml_anchor(self, any=False):
        # type: (bool) -> Any
        if not hasattr(self, Anchor.attrib):
            return None
        if any or self.anchor.always_dump:
            return self.anchor
        return None

    def yaml_set_anchor(self, value, always_dump=False):
        # type: (Any, bool) -> None
        self.anchor.value = value
        self.anchor.always_dump = always_dump


class LiteralScalarString(ScalarString):
    __slots__ = 'comment'  # the comment after the | on the first line

    style = '|'

    def __new__(cls, value, anchor=None):
        # type: (Text, Any) -> Any
        return ScalarString.__new__(cls, value, anchor=anchor)


PreservedScalarString = LiteralScalarString


class FoldedScalarString(ScalarString):
    __slots__ = ('fold_pos', 'comment')  # the comment after the > on the first line

    style = '>'

    def __new__(cls, value, anchor=None):
        # type: (Text, Any) -> Any
        return ScalarString.__new__(cls, value, anchor=anchor)


class SingleQuotedScalarString(ScalarString):
    __slots__ = ()

    style = "'"

    def __new__(cls, value, anchor=None):
        # type: (Text, Any) -> Any
        return ScalarString.__new__(cls, value, anchor=anchor)


class DoubleQuotedScalarString(ScalarString):
    __slots__ = ()

    style = '"'

    def __new__(cls, value, anchor=None):
        # type: (Text, Any) -> Any
        return ScalarString.__new__(cls, value, anchor=anchor)


class PlainScalarString(ScalarString):
    __slots__ = ()

    style = ''

    def __new__(cls, value, anchor=None):
        # type: (Text, Any) -> Any
        return ScalarString.__new__(cls, value, anchor=anchor)


def preserve_literal(s):
    # type: (Text) -> Text
    return LiteralScalarString(s.replace('\r\n', '\n').replace('\r', '\n'))


def walk_tree(base, map=None):
    # type: (Any, Any) -> None
    """
    the routine here walks over a simple yaml tree (recursing in
    dict values and list items) and converts strings that
    have multiple lines to literal scalars

    You can also provide an explicit (ordered) mapping for multiple transforms
    (first of which is executed):
        map = ruamel.yaml.compat.ordereddict
        map['\n'] = preserve_literal
        map[':'] = SingleQuotedScalarString
        walk_tree(data, map=map)
    """
    from ruamel.yaml.compat import string_types, MutableMapping, MutableSequence

    if map is None:
        map = {'\n': preserve_literal}

    if isinstance(base, MutableMapping):
        for k in base:
            v = base[k]  # type: Text
            if isinstance(v, string_types):
                for ch in map:
                    if ch in v:
                        base[k] = map[ch](v)
                        break
            else:
                walk_tree(v)
    elif isinstance(base, MutableSequence):
        for idx, elem in enumerate(base):
            if isinstance(elem, string_types):
                for ch in map:
                    if ch in elem:  # type: ignore
                        base[idx] = map[ch](elem)
                        break
            else:
                walk_tree(elem)
