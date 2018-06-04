from __future__ import unicode_literals
from .base import DEFAULT_ATTRS, Attrs

__all__ = (
    'split_token_in_parts',
    'merge_attrs',
)


def split_token_in_parts(token):
    """
    Take a Token, and turn it in a list of tokens, by splitting
    it on ':' (taking that as a separator.)
    """
    result = []
    current = []
    for part in token + (':', ):
        if part == ':':
            if current:
                result.append(tuple(current))
                current = []
        else:
            current.append(part)

    return result


def merge_attrs(list_of_attrs):
    """
    Take a list of :class:`.Attrs` instances and merge them into one.
    Every `Attr` in the list can override the styling of the previous one.
    """
    result = DEFAULT_ATTRS

    for attr in list_of_attrs:
        result = Attrs(
            color=attr.color or result.color,
            bgcolor=attr.bgcolor or result.bgcolor,
            bold=attr.bold or result.bold,
            underline=attr.underline or result.underline,
            italic=attr.italic or result.italic,
            blink=attr.blink or result.blink,
            reverse=attr.reverse or result.reverse)

    return result
