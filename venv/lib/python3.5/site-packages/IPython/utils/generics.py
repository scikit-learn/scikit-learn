# encoding: utf-8
"""Generic functions for extending IPython.

See http://pypi.python.org/pypi/simplegeneric.
"""

from IPython.core.error import TryNext
from simplegeneric import generic


@generic
def inspect_object(obj):
    """Called when you do obj?"""
    raise TryNext


@generic
def complete_object(obj, prev_completions):
    """Custom completer dispatching for python objects.

    Parameters
    ----------
    obj : object
        The object to complete.
    prev_completions : list
        List of attributes discovered so far.

    This should return the list of attributes in obj. If you only wish to
    add to the attributes already discovered normally, return
    own_attrs + prev_completions.
    """
    raise TryNext


