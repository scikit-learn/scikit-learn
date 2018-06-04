from .enums import IncrementalSearchDirection
from .filters import to_simple_filter

__all__ = (
    'SearchState',
)


class SearchState(object):
    """
    A search 'query'.
    """
    __slots__ = ('text', 'direction', 'ignore_case')

    def __init__(self, text='', direction=IncrementalSearchDirection.FORWARD, ignore_case=False):
        ignore_case = to_simple_filter(ignore_case)

        self.text = text
        self.direction = direction
        self.ignore_case = ignore_case

    def __repr__(self):
        return '%s(%r, direction=%r, ignore_case=%r)' % (
            self.__class__.__name__, self.text, self.direction, self.ignore_case)

    def __invert__(self):
        """
        Create a new SearchState where backwards becomes forwards and the other
        way around.
        """
        if self.direction == IncrementalSearchDirection.BACKWARD:
            direction = IncrementalSearchDirection.FORWARD
        else:
            direction = IncrementalSearchDirection.BACKWARD

        return SearchState(text=self.text, direction=direction, ignore_case=self.ignore_case)
