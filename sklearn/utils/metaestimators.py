"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# Licence: BSD

from operator import attrgetter
from functools import update_wrapper


__all__ = ['if_delegate_has_method']


class _IffHasAttrDescriptor(object):
    """Implements a conditional property using the descriptor protocol.

    Using this class to create a decorator will raise an ``AttributeError``
    if the ``attribute_name`` is not present on the base object.

    This allows ducktyping of the decorated method based on ``attribute_name``.

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, attribute_name):
        self.fn = fn
        self.get_attribute = attrgetter(attribute_name)
        # update the docstring of the descriptor
        update_wrapper(self, fn)

    def __get__(self, obj, type=None):
        # raise an AttributeError if the attribute is not present on the object
        if obj is not None:
            # delegate only on instances, not the classes.
            # this is to allow access to the docstrings.
            self.get_attribute(obj)
        # lambda, but not partial, allows help() to work with update_wrapper
        out = lambda *args, **kwargs: self.fn(obj, *args, **kwargs)
        # update the docstring of the returned function
        update_wrapper(out, self.fn)
        return out


def if_delegate_has_method(delegate):
    """Create a decorator for methods that are delegated to a sub-estimator

    This enables ducktyping by hasattr returning True according to the
    sub-estimator.

    >>> from sklearn.utils.metaestimators import if_delegate_has_method
    >>>
    >>>
    >>> class MetaEst(object):
    ...     def __init__(self, sub_est):
    ...         self.sub_est = sub_est
    ...
    ...     @if_delegate_has_method(delegate='sub_est')
    ...     def predict(self, X):
    ...         return self.sub_est.predict(X)
    ...
    >>> class HasPredict(object):
    ...     def predict(self, X):
    ...         return X.sum(axis=1)
    ...
    >>> class HasNoPredict(object):
    ...     pass
    ...
    >>> hasattr(MetaEst(HasPredict()), 'predict')
    True
    >>> hasattr(MetaEst(HasNoPredict()), 'predict')
    False
    """
    return lambda fn: _IffHasAttrDescriptor(fn, '%s.%s' % (delegate, fn.__name__))
