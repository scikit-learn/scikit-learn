"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD

from operator import attrgetter
from functools import update_wrapper


__all__ = ['if_delegate_has_method']


def _hasattr_nested(obj, attr):
    """Check if obj has an attribute

    Unlike `getattr`, this function allows attribute name to contain
    dots, e.g. `_hasattr_nested(obj, 'name.first')`
    """
    parts = attr.split(".")
    for part in parts:
        if hasattr(obj, part):
            obj = getattr(obj, part)
        else:
            return False
    else:
        return True


class _IffHasAttrDescriptor(object):
    """Implements a conditional property using the descriptor protocol.

    Using this class to create a decorator will raise an ``AttributeError``
    if the all items in ``attribute_name`` is not present on the base object.
    attribute_name can be a single string or a tuple of strings. The
    ``AttributeError`` raised will indicate the last item is not present on
    the base object

    This allows ducktyping of the decorated method based on ``attribute_name``.

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, attribute_name):
        self.fn = fn
        self.attribute_name = attribute_name

        # update the docstring of the descriptor
        update_wrapper(self, fn)

    def __get__(self, obj, type=None):
        # raise an AttributeError if the attribute is not present on the object
        if obj is not None:
            # delegate only on instances, not the classes.
            # this is to allow access to the docstrings.
            obj_hasattr = [_hasattr_nested(obj, item) for item in self.attribute_name]
            if not any(obj_hasattr):
                attrgetter(self.attribute_name[-1])(obj)
        # lambda, but not partial, allows help() to work with update_wrapper
        out = lambda *args, **kwargs: self.fn(obj, *args, **kwargs)
        # update the docstring of the returned function
        update_wrapper(out, self.fn)
        return out


def if_delegate_has_method(delegate):
    """Create a decorator for methods that are delegated to a sub-estimator

    Delegate can be string or a tuple of strings

    This enables ducktyping by hasattr returning True according to the
    sub-estimator.

    >>> from sklearn.utils.metaestimators import if_delegate_has_method
    >>>
    >>>
    >>> class MetaEst(object):
    ...     def __init__(self, sub_est, better_sub_est=None):
    ...         self.sub_est = sub_est
    ...         self.better_sub_est = better_sub_est
    ...
    ...     @if_delegate_has_method(delegate='sub_est')
    ...     def predict(self, X):
    ...         return self.sub_est.predict(X)
    ...
    ...     @if_delegate_has_method(delegate=('sub_est', 'better_sub_est'))
    ...     def predict_cond(self, X):
    ...         if self.better_sub_est is not None:
    ...             return self.better_sub_est.predict_cond(X)
    ...         else:
    ...             return self.sub_est.predict_cond(X)
    ...
    >>> class HasPredict(object):
    ...     def predict(self, X):
    ...         return X.sum(axis=1)
    ...
    ...     def predict_cond(self, X):
    ...         return X.sum(axis=0)
    ...
    >>> class HasNoPredict(object):
    ...     pass
    ...
    >>> hasattr(MetaEst(HasPredict()), 'predict')
    True
    >>> hasattr(MetaEst(HasNoPredict()), 'predict')
    False
    >>> hasattr(MetaEst(HasNoPredict(), HasNoPredict()), 'predict_cond')
    False
    >>> hasattr(MetaEst(HasPredict(), HasNoPredict()), 'predict_cond')
    True
    >>> hasattr(MetaEst(HasNoPredict(), HasPredict()), 'predict_cond')
    True
    >>> hasattr(MetaEst(HasPredict(), HasPredict()), 'predict_cond')
    True
    """
    if not isinstance(delegate, tuple):
        delegate = tuple([delegate])

    def func(fn):
        attrs = []
        for item in delegate:
            attrs.append('%s.%s' % (item, fn.__name__))
        return _IffHasAttrDescriptor(fn, tuple(attrs))

    return func
