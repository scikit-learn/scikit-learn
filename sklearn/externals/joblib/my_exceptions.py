"""
Exceptions
"""
# Author: Gael Varoquaux < gael dot varoquaux at normalesup dot org >
# Copyright: 2010, Gael Varoquaux
# License: BSD 3 clause

import sys


class JoblibException(Exception):
    """ A simple exception with an error message that you can get to.
    """

    def __init__(self, message):
        self.message = message

    def __reduce__(self):
        # For pickling
        return self.__class__, (self.message,), {}

    def __repr__(self):
        return '%s\n%s\n%s\n%s' % (
                    self.__class__.__name__,
                    75 * '_',
                    self.message,
                    75 * '_')

    __str__ = __repr__


class TransportableException(JoblibException):
    """ An exception containing all the info to wrap an original
        exception and recreate it.
    """

    def __init__(self, message, etype):
        self.message = message
        self.etype = etype

    def __reduce__(self):
        # For pickling
        return self.__class__, (self.message, self.etype), {}


_exception_mapping = dict()


def _mk_exception(exception, name=None):
    # Create an exception inheriting from both JoblibException
    # and that exception
    if name is None:
        name = exception.__name__
    this_name = 'Joblib%s' % name
    if this_name in _exception_mapping:
        # Avoid creating twice the same exception
        this_exception = _exception_mapping[this_name]
    else:
        this_exception = type(this_name, (exception, JoblibException),
                    dict(__repr__=JoblibException.__repr__,
                         __str__=JoblibException.__str__),
                    )
        _exception_mapping[this_name] = this_exception
    return this_exception, this_name


def _mk_common_exceptions():
    namespace = dict()
    if sys.version_info[0] == 3:
        import builtins as _builtin_exceptions
        common_exceptions = filter(
            lambda x: x.endswith('Error'),
            dir(_builtin_exceptions))
    else:
        import exceptions as _builtin_exceptions
        common_exceptions = dir(_builtin_exceptions)

    for name in common_exceptions:
        obj = getattr(_builtin_exceptions, name)
        if isinstance(obj, type) and issubclass(obj, BaseException):
            try:
                this_obj, this_name = _mk_exception(obj, name=name)
                namespace[this_name] = this_obj
            except TypeError:
                # Cannot create a consistent method resolution order:
                # a class that we can't subclass properly, probably
                # BaseException
                pass
    return namespace


# Updating module locals so that the exceptions pickle right. AFAIK this
# works only at module-creation time
locals().update(_mk_common_exceptions())
