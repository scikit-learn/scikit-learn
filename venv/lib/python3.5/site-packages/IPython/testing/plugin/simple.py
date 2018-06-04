"""Simple example using doctests.

This file just contains doctests both using plain python and IPython prompts.
All tests should be loaded by nose.
"""

def pyfunc():
    """Some pure python tests...

    >>> pyfunc()
    'pyfunc'

    >>> import os

    >>> 2+3
    5

    >>> for i in range(3):
    ...     print(i, end=' ')
    ...     print(i+1, end=' ')
    ...
    0 1 1 2 2 3 
    """
    return 'pyfunc'


def ipyfunc2():
    """Some pure python tests...

    >>> 1+1
    2
    """
    return 'pyfunc2'
