"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD

from operator import attrgetter
from functools import update_wrapper


__all__ = ['if_delegate_has_method']


class _IffHasAttrDescriptor(object):
    """Implements a conditional property using the descriptor protocol.

    Using this class to create a decorator will raise an ``AttributeError``
    if none of the delegates (specified in ``delegate_names``) is an attribute
    of the base object or none of the delegates has an attribute
    ``method_name``.

    This allows ducktyping of the decorated method based on
    ``delegate.method_name`` where ``delegate`` is the first item in
    ``delegate_names`` that is an attribute of the base object.

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, delegate_names, method_name):
        self.fn = fn
        self.delegate_names = delegate_names
        self.method_name = method_name

        # update the docstring of the descriptor
        update_wrapper(self, fn)

    def __get__(self, obj, type=None):
        # raise an AttributeError if the attribute is not present on the object
        if obj is not None:
            # delegate only on instances, not the classes.
            # this is to allow access to the docstrings.
            for delegate_name in self.delegate_names:
                try:
                    delegate = attrgetter(delegate_name)(obj)
                except AttributeError:
                    continue
                else:
                    getattr(delegate, self.method_name)
                    break
            else:
                attrgetter(self.delegate_names[-1])(obj)

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
    delegate : string, list of strings or tuple of strings
        Name of the sub-estimator that can be accessed as an attribute of the
        base object. If a list or a tuple of names are provided, the first
        sub-estimator that is an attribute of the base object  will be used.

    Examples
    --------
    >>> from sklearn.utils.metaestimators import if_delegate_has_method
    >>>
    >>>
    >>> class MetaEst(object):
    ...     def __init__(self, sub_est, better_sub_est=None):
    ...         self.sub_est = sub_est
    ...         self.better_sub_est = better_sub_est
    ...
    ...     @if_delegate_has_method(delegate='sub_est')
    ...     def predict(self):
    ...         pass
    ...
    >>> class MetaEstTestTuple(MetaEst):
    ...     @if_delegate_has_method(delegate=('sub_est', 'better_sub_est'))
    ...     def predict(self):
    ...         pass
    ...
    >>> class MetaEstTestList(MetaEst):
    ...     @if_delegate_has_method(delegate=['sub_est', 'better_sub_est'])
    ...     def predict(self):
    ...         pass
    ...
    >>> class HasPredict(object):
    ...     def predict(self):
    ...         pass
    ...
    >>> class HasNoPredict(object):
    ...     pass
    ...
    >>> hasattr(MetaEst(HasPredict()), 'predict')
    True
    >>> hasattr(MetaEst(HasNoPredict()), 'predict')
    False
    >>> hasattr(MetaEstTestTuple(HasNoPredict(), HasNoPredict()), 'predict')
    False
    >>> hasattr(MetaEstTestTuple(HasPredict(), HasNoPredict()), 'predict')
    True
    >>> hasattr(MetaEstTestTuple(HasNoPredict(), HasPredict()), 'predict')
    False
    >>> hasattr(MetaEstTestList(HasNoPredict(), HasPredict()), 'predict')
    False
    >>> hasattr(MetaEstTestList(HasPredict(), HasPredict()), 'predict')
    True

    """
    if isinstance(delegate, list):
        delegate = tuple(delegate)
    if not isinstance(delegate, tuple):
        delegate = (delegate,)

    return lambda fn: _IffHasAttrDescriptor(fn, delegate,
                                            method_name=fn.__name__)
