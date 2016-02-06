"""
Compatibility layer for Python 3/Python 2 single codebase
"""
import sys

PY3_OR_LATER = sys.version_info[0] >= 3

try:
    _basestring = basestring
    _bytes_or_unicode = (str, unicode)
except NameError:
    _basestring = str
    _bytes_or_unicode = (bytes, str)
