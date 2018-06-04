"""
Clipboard for command line interface.
"""
from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from six import with_metaclass
import six

from prompt_toolkit.selection import SelectionType

__all__ = (
    'Clipboard',
    'ClipboardData',
)


class ClipboardData(object):
    """
    Text on the clipboard.

    :param text: string
    :param type: :class:`~prompt_toolkit.selection.SelectionType`
    """
    def __init__(self, text='', type=SelectionType.CHARACTERS):
        assert isinstance(text, six.string_types)
        assert type in (SelectionType.CHARACTERS, SelectionType.LINES, SelectionType.BLOCK)

        self.text = text
        self.type = type


class Clipboard(with_metaclass(ABCMeta, object)):
    """
    Abstract baseclass for clipboards.
    (An implementation can be in memory, it can share the X11 or Windows
    keyboard, or can be persistent.)
    """
    @abstractmethod
    def set_data(self, data):
        """
        Set data to the clipboard.

        :param data: :class:`~.ClipboardData` instance.
        """

    def set_text(self, text):  # Not abstract.
        """
        Shortcut for setting plain text on clipboard.
        """
        assert isinstance(text, six.string_types)
        self.set_data(ClipboardData(text))

    def rotate(self):
        """
        For Emacs mode, rotate the kill ring.
        """

    @abstractmethod
    def get_data(self):
        """
        Return clipboard data.
        """
