"""
Compatibility layer for Python 3/Python 2 single codebase
"""

try:
    _basestring = basestring
except NameError:
    _basestring = str

