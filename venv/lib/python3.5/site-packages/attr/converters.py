"""
Commonly useful converters.
"""

from __future__ import absolute_import, division, print_function


def optional(converter):
    """
    A converter that allows an attribute to be optional. An optional attribute
    is one which can be set to ``None``.

    :param callable converter: the converter that is used for non-``None``
        values.

    ..  versionadded:: 17.1.0
    """

    def optional_converter(val):
        if val is None:
            return None
        return converter(val)

    return optional_converter
