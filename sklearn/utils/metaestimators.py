"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD

from operator import attrgetter
from functools import update_wrapper
from .validation import check_is_fitted

__all__ = ['if_delegate_has_method']


class _IffHasAttrDescriptor(object):
    """Implements a conditional property using the descriptor protocol.

    Using this class to create a decorator will raise an ``AttributeError``
    if the delegate doesn't have an attribute ``attribute_name``.

    This allows ducktyping of the decorated method based on
    ``delegate.attribute_name``

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, delegate_name, attribute_name):
        self.fn = fn
        self.delegate_name = delegate_name
        self.attribute_name = attribute_name

        # update the docstring of the descriptor
        update_wrapper(self, fn)

    def __get__(self, obj, type=None):
        # raise an AttributeError if the attribute is not present on the object
        if obj is not None:
            # delegate only on instances, not the classes.
            # this is to allow access to the docstrings.
            check_is_fitted(obj, self.delegate_name)
            delegate = attrgetter(self.delegate_name)(obj)
            getattr(delegate, self.attribute_name)

        # lambda, but not partial, allows help() to work with update_wrapper
        out = lambda *args, **kwargs: self.fn(obj, *args, **kwargs)
        # update the docstring of the returned function
        update_wrapper(out, self.fn)
        return out


def if_delegate_has_method(delegate):
    """Create a decorator for methods that are delegated to a sub-estimator

    This enables ducktyping by hasattr returning True according to the
    sub-estimator.

    Parameters
    ----------
    delegate : string
        Name of the sub-estimator that can be accessed as an attribute of the
        base object.

    """

    return lambda fn: _IffHasAttrDescriptor(fn, delegate,
                                            attribute_name=fn.__name__)
