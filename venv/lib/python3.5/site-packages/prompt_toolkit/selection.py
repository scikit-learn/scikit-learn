"""
Data structures for the selection.
"""
from __future__ import unicode_literals

__all__ = (
    'SelectionType',
    'PasteMode',
    'SelectionState',
)


class SelectionType(object):
    """
    Type of selection.
    """
    #: Characters. (Visual in Vi.)
    CHARACTERS = 'CHARACTERS'

    #: Whole lines. (Visual-Line in Vi.)
    LINES = 'LINES'

    #: A block selection. (Visual-Block in Vi.)
    BLOCK = 'BLOCK'


class PasteMode(object):
    EMACS = 'EMACS'  # Yank like emacs.
    VI_AFTER = 'VI_AFTER'  # When pressing 'p' in Vi.
    VI_BEFORE = 'VI_BEFORE'  # When pressing 'P' in Vi.


class SelectionState(object):
    """
    State of the current selection.

    :param original_cursor_position: int
    :param type: :class:`~.SelectionType`
    """
    def __init__(self, original_cursor_position=0, type=SelectionType.CHARACTERS):
        self.original_cursor_position = original_cursor_position
        self.type = type

    def __repr__(self):
        return '%s(original_cursor_position=%r, type=%r)' % (
            self.__class__.__name__,
            self.original_cursor_position, self.type)
