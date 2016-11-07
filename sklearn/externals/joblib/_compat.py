"""
Compatibility layer for Python 3/Python 2 single codebase
"""
import sys

PY3_OR_LATER = sys.version_info[0] >= 3
PY26 = sys.version_info[:2] == (2, 6)
PY27 = sys.version_info[:2] == (2, 7)

try:
    _basestring = basestring
    _bytes_or_unicode = (str, unicode)
except NameError:
    _basestring = str
    _bytes_or_unicode = (bytes, str)


def with_metaclass(meta, *bases):
    """Create a base class with a metaclass."""
    return meta("NewBase", bases, {})
