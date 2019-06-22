# coding: utf-8

from __future__ import print_function, absolute_import, division, unicode_literals

import sys
from .compat import no_limit_int  # NOQA
from ruamel.yaml.anchor import Anchor

if False:  # MYPY
    from typing import Text, Any, Dict, List  # NOQA

__all__ = ['ScalarFloat', 'ExponentialFloat', 'ExponentialCapsFloat']


class ScalarFloat(float):
    def __new__(cls, *args, **kw):
        # type: (Any, Any, Any) -> Any
        width = kw.pop('width', None)  # type: ignore
        prec = kw.pop('prec', None)  # type: ignore
        m_sign = kw.pop('m_sign', None)  # type: ignore
        m_lead0 = kw.pop('m_lead0', 0)  # type: ignore
        exp = kw.pop('exp', None)  # type: ignore
        e_width = kw.pop('e_width', None)  # type: ignore
        e_sign = kw.pop('e_sign', None)  # type: ignore
        underscore = kw.pop('underscore', None)  # type: ignore
        anchor = kw.pop('anchor', None)  # type: ignore
        v = float.__new__(cls, *args, **kw)  # type: ignore
        v._width = width
        v._prec = prec
        v._m_sign = m_sign
        v._m_lead0 = m_lead0
        v._exp = exp
        v._e_width = e_width
        v._e_sign = e_sign
        v._underscore = underscore
        if anchor is not None:
            v.yaml_set_anchor(anchor, always_dump=True)
        return v

    def __iadd__(self, a):  # type: ignore
        # type: (Any) -> Any
        return float(self) + a
        x = type(self)(self + a)
        x._width = self._width
        x._underscore = self._underscore[:] if self._underscore is not None else None  # NOQA
        return x

    def __ifloordiv__(self, a):  # type: ignore
        # type: (Any) -> Any
        return float(self) // a
        x = type(self)(self // a)
        x._width = self._width
        x._underscore = self._underscore[:] if self._underscore is not None else None  # NOQA
        return x

    def __imul__(self, a):  # type: ignore
        # type: (Any) -> Any
        return float(self) * a
        x = type(self)(self * a)
        x._width = self._width
        x._underscore = self._underscore[:] if self._underscore is not None else None  # NOQA
        x._prec = self._prec  # check for others
        return x

    def __ipow__(self, a):  # type: ignore
        # type: (Any) -> Any
        return float(self) ** a
        x = type(self)(self ** a)
        x._width = self._width
        x._underscore = self._underscore[:] if self._underscore is not None else None  # NOQA
        return x

    def __isub__(self, a):  # type: ignore
        # type: (Any) -> Any
        return float(self) - a
        x = type(self)(self - a)
        x._width = self._width
        x._underscore = self._underscore[:] if self._underscore is not None else None  # NOQA
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

    def dump(self, out=sys.stdout):
        # type: (Any) -> Any
        out.write(
            'ScalarFloat({}| w:{}, p:{}, s:{}, lz:{}, _:{}|{}, w:{}, s:{})\n'.format(
                self,
                self._width,  # type: ignore
                self._prec,  # type: ignore
                self._m_sign,  # type: ignore
                self._m_lead0,  # type: ignore
                self._underscore,  # type: ignore
                self._exp,  # type: ignore
                self._e_width,  # type: ignore
                self._e_sign,  # type: ignore
            )
        )


class ExponentialFloat(ScalarFloat):
    def __new__(cls, value, width=None, underscore=None):
        # type: (Any, Any, Any) -> Any
        return ScalarFloat.__new__(cls, value, width=width, underscore=underscore)


class ExponentialCapsFloat(ScalarFloat):
    def __new__(cls, value, width=None, underscore=None):
        # type: (Any, Any, Any) -> Any
        return ScalarFloat.__new__(cls, value, width=width, underscore=underscore)
