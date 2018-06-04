"""
The base classes for the styling.
"""
from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from collections import namedtuple
from six import with_metaclass

__all__ = (
    'Attrs',
    'DEFAULT_ATTRS',
    'ANSI_COLOR_NAMES',
    'Style',
    'DynamicStyle',
)


#: Style attributes.
Attrs = namedtuple('Attrs', 'color bgcolor bold underline italic blink reverse')
"""
:param color: Hexadecimal string. E.g. '000000' or Ansi color name: e.g. 'ansiblue'
:param bgcolor: Hexadecimal string. E.g. 'ffffff' or Ansi color name: e.g. 'ansired'
:param bold: Boolean
:param underline: Boolean
:param italic: Boolean
:param blink: Boolean
:param reverse: Boolean
"""

#: The default `Attrs`.
DEFAULT_ATTRS = Attrs(color=None, bgcolor=None, bold=False, underline=False,
                      italic=False, blink=False, reverse=False)


#: ``Attrs.bgcolor/fgcolor`` can be in either 'ffffff' format, or can be any of
#: the following in case we want to take colors from the 8/16 color palette.
#: Usually, in that case, the terminal application allows to configure the RGB
#: values for these names.
ANSI_COLOR_NAMES = [
    'ansiblack', 'ansiwhite', 'ansidefault',

    # Low intensity.
    'ansired', 'ansigreen', 'ansiyellow', 'ansiblue', 'ansifuchsia', 'ansiturquoise', 'ansilightgray',

    # High intensity. (Not supported everywhere.)
    'ansidarkgray', 'ansidarkred', 'ansidarkgreen', 'ansibrown', 'ansidarkblue',
    'ansipurple', 'ansiteal',
]


class Style(with_metaclass(ABCMeta, object)):
    """
    Abstract base class for prompt_toolkit styles.
    """
    @abstractmethod
    def get_attrs_for_token(self, token):
        """
        Return :class:`.Attrs` for the given token.
        """

    @abstractmethod
    def invalidation_hash(self):
        """
        Invalidation hash for the style. When this changes over time, the
        renderer knows that something in the style changed, and that everything
        has to be redrawn.
        """


class DynamicStyle(Style):
    """
    Style class that can dynamically returns an other Style.

    :param get_style: Callable that returns a :class:`.Style` instance.
    """
    def __init__(self, get_style):
        self.get_style = get_style

    def get_attrs_for_token(self, token):
        style = self.get_style()
        assert isinstance(style, Style)

        return style.get_attrs_for_token(token)

    def invalidation_hash(self):
        return self.get_style().invalidation_hash()
