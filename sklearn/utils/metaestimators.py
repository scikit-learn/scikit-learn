"""Utilities for meta-estimators"""
# Author: Joel Nothman
# Licence: BSD

from operator import attrgetter
from functools import partial, update_wrapper


__all__ = ['make_delegation_decorator']


class _IffHasAttrDescriptor(object):
    def __init__(self, fn, attr):
        self.fn = fn
        self.get_attr = attrgetter(attr)

    def __get__(self, obj, type=None):
        self.get_attr(obj)
        # lambda, but not partial, allows help() to work with update_wrapper
        out = lambda *args, **kwargs: self.fn(obj, *args, **kwargs)
        update_wrapper(out, self.fn)
        return out


def _iff_attr_has_method(prefix, fn, name=None):
    if not callable(fn):
        return partial(_iff_attr_has_method, prefix, name=fn)
    if name is None:
        name = fn.__name__
    attr = '%s.%s' % (prefix, name)
    out = _IffHasAttrDescriptor(fn, attr)
    update_wrapper(out, fn)
    return out


def make_delegation_decorator(prefix):
    """Create a decorator for methods that are delegated to a sub-estimator

    This enables ducktyping by hasattr returning True according to the
    sub-estimator.

    >>> from sklearn.utils.metaestimators import make_delegation_decorator
    >>>
    >>> iff_sub_est_has_attr = make_delegation_decorator('sub_est')
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
    return partial(_iff_attr_has_method, prefix)
