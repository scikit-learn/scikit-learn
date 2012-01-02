"""
Test the func_inspect module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import nose
import tempfile
import functools

from ..func_inspect import filter_args, get_func_name, get_func_code
from ..memory import Memory


###############################################################################
# Module-level functions, for tests
def f(x, y=0):
    pass


def f2(x):
    pass


# Create a Memory object to test decorated functions.
# We should be careful not to call the decorated functions, so that
# cache directories are not created in the temp dir.
mem = Memory(cachedir=tempfile.gettempdir())


@mem.cache
def g(x):
    return x


def h(x, y=0, *args, **kwargs):
    pass


def i(x=1):
    pass


def j(x, y, **kwargs):
    pass


def k(*args, **kwargs):
    pass


class Klass(object):

    def f(self, x):
        return x


###############################################################################
# Tests

def test_filter_args():
    yield nose.tools.assert_equal, filter_args(f, [], 1), {'x': 1, 'y': 0}
    yield nose.tools.assert_equal, filter_args(f, ['x'], 1), {'y': 0}
    yield nose.tools.assert_equal, filter_args(f, ['y'], 0), {'x': 0}
    yield nose.tools.assert_equal, filter_args(f, ['y'], 0, y=1), {'x': 0}
    yield nose.tools.assert_equal, filter_args(f, ['x', 'y'], 0), {}
    yield nose.tools.assert_equal, filter_args(f, [], 0, y=1), {'x': 0, 'y': 1}
    yield nose.tools.assert_equal, filter_args(f, ['y'], x=2, y=1), {'x': 2}

    yield nose.tools.assert_equal, filter_args(i, [], 2), {'x': 2}
    yield nose.tools.assert_equal, filter_args(f2, [], x=1), {'x': 1}


def test_filter_args_method():
    obj = Klass()
    nose.tools.assert_equal(filter_args(obj.f, [], 1),
        {'x': 1, 'self': obj})


def test_filter_varargs():
    yield nose.tools.assert_equal, filter_args(h, [], 1), \
                            {'x': 1, 'y': 0, '*': [], '**': {}}
    yield nose.tools.assert_equal, filter_args(h, [], 1, 2, 3, 4), \
                            {'x': 1, 'y': 2, '*': [3, 4], '**': {}}
    yield nose.tools.assert_equal, filter_args(h, [], 1, 25, ee=2), \
                            {'x': 1, 'y': 25, '*': [], '**': {'ee': 2}}
    yield nose.tools.assert_equal, filter_args(h, ['*'], 1, 2, 25, ee=2), \
                            {'x': 1, 'y': 2, '**': {'ee': 2}}


def test_filter_kwargs():
    nose.tools.assert_equal(filter_args(k, [], 1, 2, ee=2),
                            {'*': [1, 2], '**': {'ee': 2}})
    nose.tools.assert_equal(filter_args(k, [], 3, 4),
                            {'*': [3, 4], '**': {}})


def test_filter_args_2():
    nose.tools.assert_equal(filter_args(j, [], 1, 2, ee=2),
                            {'x': 1, 'y': 2, '**': {'ee': 2}})

    nose.tools.assert_raises(ValueError, filter_args, f, 'a', None)
    # Check that we capture an undefined argument
    nose.tools.assert_raises(ValueError, filter_args, f, ['a'], None)
    ff = functools.partial(f, 1)
    # filter_args has to special-case partial
    nose.tools.assert_equal(filter_args(ff, [], 1),
                            {'*': [1], '**': {}})
    nose.tools.assert_equal(filter_args(ff, ['y'], 1),
                            {'*': [1], '**': {}})


def test_func_name():
    yield nose.tools.assert_equal, 'f', get_func_name(f)[1]
    # Check that we are not confused by the decoration
    yield nose.tools.assert_equal, 'g', get_func_name(g)[1]


def test_func_inspect_errors():
    """ Check that func_inspect is robust and will work on weird objects
    """
    nose.tools.assert_equal(get_func_name('a'.lower)[-1], 'lower')
    nose.tools.assert_equal(get_func_code('a'.lower)[1:], (None, -1))
    ff = lambda x: x
    nose.tools.assert_equal(get_func_name(ff, win_characters=False)[-1],
                                                            '<lambda>')
    nose.tools.assert_equal(get_func_code(ff)[1],
                                    __file__.replace('.pyc', '.py'))
    # Simulate a function defined in __main__
    ff.__module__ = '__main__'
    nose.tools.assert_equal(get_func_name(ff, win_characters=False)[-1],
                                                            '<lambda>')
    nose.tools.assert_equal(get_func_code(ff)[1],
                                    __file__.replace('.pyc', '.py'))


def test_bound_methods():
    """ Make sure that calling the same method on two different instances
        of the same class does resolv to different signatures.
    """
    a = Klass()
    b = Klass()
    nose.tools.assert_not_equal(filter_args(a.f, [], 1),
                                filter_args(b.f, [], 1))


def test_filter_args_error_msg():
    """ Make sure that filter_args returns decent error messages, for the
        sake of the user.
    """
    nose.tools.assert_raises(ValueError, filter_args, f, [])
