"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# Licence: BSD

from operator import attrgetter
from functools import update_wrapper


__all__ = ['if_delegate_has_method']


class _IffHasAttrDescriptor(object):
    """Implements a conditional propery using the descriptor protocol.

    Using this class to create a decorator will raise an attribute error
    if the ``attr`` is not present on the base object.

    This allows ducktyping of the decorated method based on ``attr``.

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, attr):
        self.fn = fn
        self.get_attr = attrgetter(attr)
        # update the docstring of the descriptor
        update_wrapper(self, fn)

    def __get__(self, obj, type=None):
        # raise an AttributeError if the attribute is not present on the object
        if obj is not None:
            # delegate only on instances, not the classes.
            # this is to allow access to the docstrings.
            self.get_attr(obj)
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
    >>> iff_sub_est_has_attr = if_delegate_has_method(delegate='sub_est')
    >>>
    >>> class MetaEst(object):
    ...     def __init__(self, sub_est):
    ...         self.sub_est = sub_est
    ...
    ...     @iff_sub_est_has_attr
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
