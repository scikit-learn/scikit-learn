# coding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

from .compat import no_limit_int  # NOQA
from ruamel.yaml.anchor import Anchor

if False:  # MYPY
    from typing import Text, Any, Dict, List  # NOQA

__all__ = ['ScalarInt', 'BinaryInt', 'OctalInt', 'HexInt', 'HexCapsInt', 'DecimalInt']


class ScalarInt(no_limit_int):
    def __new__(cls, *args, **kw):
        # type: (Any, Any, Any) -> Any
        width = kw.pop('width', None)  # type: ignore
        underscore = kw.pop('underscore', None)  # type: ignore
        anchor = kw.pop('anchor', None)  # type: ignore
        v = no_limit_int.__new__(cls, *args, **kw)  # type: ignore
        v._width = width
        v._underscore = underscore
        if anchor is not None:
            v.yaml_set_anchor(anchor, always_dump=True)
        return v

    def __iadd__(self, a):  # type: ignore
        # type: (Any) -> Any
        x = type(self)(self + a)
        x._width = self._width  # type: ignore
        x._underscore = (  # type: ignore
            self._underscore[:] if self._underscore is not None else None  # type: ignore
        )  # NOQA
        return x

    def __ifloordiv__(self, a):  # type: ignore
        # type: (Any) -> Any
        x = type(self)(self // a)
        x._width = self._width  # type: ignore
        x._underscore = (  # type: ignore
            self._underscore[:] if self._underscore is not None else None  # type: ignore
        )  # NOQA
        return x

    def __imul__(self, a):  # type: ignore
        # type: (Any) -> Any
        x = type(self)(self * a)
        x._width = self._width  # type: ignore
        x._underscore = (  # type: ignore
            self._underscore[:] if self._underscore is not None else None  # type: ignore
        )  # NOQA
        return x

    def __ipow__(self, a):  # type: ignore
        # type: (Any) -> Any
        x = type(self)(self ** a)
        x._width = self._width  # type: ignore
        x._underscore = (  # type: ignore
            self._underscore[:] if self._underscore is not None else None  # type: ignore
        )  # NOQA
        return x

    def __isub__(self, a):  # type: ignore
        # type: (Any) -> Any
        x = type(self)(self - a)
        x._width = self._width  # type: ignore
        x._underscore = (  # type: ignore
            self._underscore[:] if self._underscore is not None else None  # type: ignore
        )  # NOQA
        return x

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


class BinaryInt(ScalarInt):
    def __new__(cls, value, width=None, underscore=None, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        return ScalarInt.__new__(cls, value, width=width, underscore=underscore, anchor=anchor)


class OctalInt(ScalarInt):
    def __new__(cls, value, width=None, underscore=None, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        return ScalarInt.__new__(cls, value, width=width, underscore=underscore, anchor=anchor)


# mixed casing of A-F is not supported, when loading the first non digit
# determines the case


class HexInt(ScalarInt):
    """uses lower case (a-f)"""

    def __new__(cls, value, width=None, underscore=None, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        return ScalarInt.__new__(cls, value, width=width, underscore=underscore, anchor=anchor)


class HexCapsInt(ScalarInt):
    """uses upper case (A-F)"""

    def __new__(cls, value, width=None, underscore=None, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        return ScalarInt.__new__(cls, value, width=width, underscore=underscore, anchor=anchor)


class DecimalInt(ScalarInt):
    """needed if anchor"""

    def __new__(cls, value, width=None, underscore=None, anchor=None):
        # type: (Any, Any, Any, Any) -> Any
        return ScalarInt.__new__(cls, value, width=width, underscore=underscore, anchor=anchor)
