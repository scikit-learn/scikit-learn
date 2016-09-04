"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD

from operator import attrgetter
from functools import update_wrapper
from .validation import check_is_fitted

__all__ = ['if_delegate_has_method',
           'if_fitted_delegate_has_method']


class _IffHasAttrDescriptor(object):
    """Implements a conditional property using the descriptor protocol.

    Using this class to create a decorator will raise an ``AttributeError``
    if the delegate doesn't have an attribute ``attribute_name``.

    This allows ducktyping of the decorated method based on
    ``delegate.attribute_name``

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, delegate_name, attribute_name, fn_check_obj):
        self.fn = fn
        self.delegate_name = delegate_name
        self.attribute_name = attribute_name
        self.fn_check_obj = fn_check_obj

        # update the docstring of the descriptor
        update_wrapper(self, fn)

    def __get__(self, obj, type=None):
        # raise an AttributeError if the attribute is not present on the object
        if obj is not None:
            # delegate only on instances, not the classes.
            # this is to allow access to the docstrings.
            if self.fn_check_obj is not None:
                self.fn_check_obj(obj)
            delegate = attrgetter(self.delegate_name)(obj)
            getattr(delegate, self.attribute_name)

        # lambda, but not partial, allows help() to work with update_wrapper
        out = lambda *args, **kwargs: self.fn(obj, *args, **kwargs)
        # update the docstring of the returned function
        update_wrapper(out, self.fn)
        return out


def if_delegate_has_method(delegate, fn_check_obj=None):
    """Create a decorator for methods that are delegated to a sub-estimator

    This enables ducktyping by hasattr returning True according to the
    sub-estimator.

    Parameters
    ----------
    delegate : string
        Name of the sub-estimator that can be accessed as an attribute of the
        base object.

    fn_check_obj : function
        fn_check_obj(obj) with obj as the metaestimator instance.
        The function will be executed before ducktyping. Typically used for
        checking if the estimator has been fitted.

    """

    return lambda fn: _IffHasAttrDescriptor(fn, delegate,
                                            attribute_name=fn.__name__,
                                            fn_check_obj=fn_check_obj)


def if_fitted_delegate_has_method(delegate, *args, **kwargs):
    """Check if metaestimator has been fitted before ducktyping

    Parameters
    ----------
    delegate : string
        Name of the sub-estimator that can be accessed as an attribute of the
        base object.

    """
    fn = lambda obj: check_is_fitted(obj, *args, **kwargs)
    return if_delegate_has_method(delegate, fn_check_obj=fn)
