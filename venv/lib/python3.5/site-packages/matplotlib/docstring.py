from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from matplotlib import cbook
import sys
import types


class Substitution(object):
    """
    A decorator to take a function's docstring and perform string
    substitution on it.

    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter)

    Usage: construct a docstring.Substitution with a sequence or
    dictionary suitable for performing substitution; then
    decorate a suitable function with the constructed object. e.g.

    sub_author_name = Substitution(author='Jason')

    @sub_author_name
    def some_function(x):
        "%(author)s wrote this function"

    # note that some_function.__doc__ is now "Jason wrote this function"

    One can also use positional arguments.

    sub_first_last_names = Substitution('Edgar Allen', 'Poe')

    @sub_first_last_names
    def some_function(x):
        "%s %s wrote the Raven"
    """
    def __init__(self, *args, **kwargs):
        assert not (len(args) and len(kwargs)), \
                "Only positional or keyword args are allowed"
        self.params = args or kwargs

    def __call__(self, func):
        func.__doc__ = func.__doc__ and func.__doc__ % self.params
        return func

    def update(self, *args, **kwargs):
        "Assume self.params is a dict and update it with supplied args"
        self.params.update(*args, **kwargs)

    @classmethod
    def from_params(cls, params):
        """
        In the case where the params is a mutable sequence (list or
        dictionary) and it may change before this class is called, one may
        explicitly use a reference to the params rather than using *args or
        **kwargs which will copy the values and not reference them.
        """
        result = cls()
        result.params = params
        return result


class Appender(object):
    """
    A function decorator that will append an addendum to the docstring
    of the target function.

    This decorator should be robust even if func.__doc__ is None
    (for example, if -OO was passed to the interpreter).

    Usage: construct a docstring.Appender with a string to be joined to
    the original docstring. An optional 'join' parameter may be supplied
    which will be used to join the docstring and addendum. e.g.

    add_copyright = Appender("Copyright (c) 2009", join='\n')

    @add_copyright
    def my_dog(has='fleas'):
        "This docstring will have a copyright below"
        pass
    """
    def __init__(self, addendum, join=''):
        self.addendum = addendum
        self.join = join

    def __call__(self, func):
        docitems = [func.__doc__, self.addendum]
        func.__doc__ = func.__doc__ and self.join.join(docitems)
        return func


def dedent(func):
    "Dedent a docstring (if present)"
    func.__doc__ = func.__doc__ and cbook.dedent(func.__doc__)
    return func


def copy(source):
    "Copy a docstring from another source function (if present)"
    def do_copy(target):
        if source.__doc__:
            target.__doc__ = source.__doc__
        return target
    return do_copy

# create a decorator that will house the various documentation that
#  is reused throughout matplotlib
interpd = Substitution()


def dedent_interpd(func):
    """A special case of the interpd that first performs a dedent on
    the incoming docstring"""
    if isinstance(func, types.MethodType) and not six.PY3:
        func = func.im_func
    return interpd(dedent(func))


def copy_dedent(source):
    """A decorator that will copy the docstring from the source and
    then dedent it"""
    # note the following is ugly because "Python is not a functional
    # language" - GVR. Perhaps one day, functools.compose will exist.
    #  or perhaps not.
    #  http://mail.python.org/pipermail/patches/2007-February/021687.html
    return lambda target: dedent(copy(source)(target))
