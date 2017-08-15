"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD

from abc import ABCMeta, abstractmethod
from operator import attrgetter
from functools import update_wrapper
import re
import fnmatch

import numpy as np

from ..utils import safe_indexing
from ..externals import six
from ..base import BaseEstimator

__all__ = ['if_delegate_has_method', 'check_routing']


class _BaseComposition(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Handles parameter management for classifiers composed of named estimators.
    """
    @abstractmethod
    def __init__(self):
        pass

    def _get_params(self, attr, deep=True):
        out = super(_BaseComposition, self).get_params(deep=False)
        if not deep:
            return out
        estimators = getattr(self, attr)
        out.update(estimators)
        for name, estimator in estimators:
            if estimator is None:
                continue
            for key, value in six.iteritems(estimator.get_params(deep=True)):
                out['%s__%s' % (name, key)] = value
        return out

    def _set_params(self, attr, **params):
        # Ensure strict ordering of parameter setting:
        # 1. All steps
        if attr in params:
            setattr(self, attr, params.pop(attr))
        # 2. Step replacement
        names, _ = zip(*getattr(self, attr))
        for name in list(six.iterkeys(params)):
            if '__' not in name and name in names:
                self._replace_estimator(attr, name, params.pop(name))
        # 3. Step parameters and other initilisation arguments
        super(_BaseComposition, self).set_params(**params)
        return self

    def _replace_estimator(self, attr, name, new_val):
        # assumes `name` is a valid estimator name
        new_estimators = getattr(self, attr)[:]
        for i, (estimator_name, _) in enumerate(new_estimators):
            if estimator_name == name:
                new_estimators[i] = (name, new_val)
                break
        setattr(self, attr, new_estimators)

    def _validate_names(self, names):
        if len(set(names)) != len(names):
            raise ValueError('Names provided are not unique: '
                             '{0!r}'.format(list(names)))
        invalid_names = set(names).intersection(self.get_params(deep=False))
        if invalid_names:
            raise ValueError('Estimator names conflict with constructor '
                             'arguments: {0!r}'.format(sorted(invalid_names)))
        invalid_names = [name for name in names if '__' in name]
        if invalid_names:
            raise ValueError('Estimator names must not contain __: got '
                             '{0!r}'.format(invalid_names))


class _IffHasAttrDescriptor(object):
    """Implements a conditional property using the descriptor protocol.

    Using this class to create a decorator will raise an ``AttributeError``
    if none of the delegates (specified in ``delegate_names``) is an attribute
    of the base object or the first found delegate does not have an attribute
    ``attribute_name``.

    This allows ducktyping of the decorated method based on
    ``delegate.attribute_name``. Here ``delegate`` is the first item in
    ``delegate_names`` for which ``hasattr(object, delegate) is True``.

    See https://docs.python.org/3/howto/descriptor.html for an explanation of
    descriptors.
    """
    def __init__(self, fn, delegate_names, attribute_name):
        self.fn = fn
        self.delegate_names = delegate_names
        self.attribute_name = attribute_name

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
                    getattr(delegate, self.attribute_name)
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
        sub-estimator that is an attribute of the base object will be used.

    """
    if isinstance(delegate, list):
        delegate = tuple(delegate)
    if not isinstance(delegate, tuple):
        delegate = (delegate,)

    return lambda fn: _IffHasAttrDescriptor(fn, delegate,
                                            attribute_name=fn.__name__)


def _safe_split(estimator, X, y, indices, train_indices=None):
    """Create subset of dataset and properly handle kernels.

    Slice X, y according to indices for cross-validation, but take care of
    precomputed kernel-matrices or pairwise affinities / distances.

    If ``estimator._pairwise is True``, X needs to be square and
    we slice rows and columns. If ``train_indices`` is not None,
    we slice rows using ``indices`` (assumed the test set) and columns
    using ``train_indices``, indicating the training set.

    Labels y will always be sliced only along the last axis.

    Parameters
    ----------
    estimator : object
        Estimator to determine whether we should slice only rows or rows and
        columns.

    X : array-like, sparse matrix or iterable
        Data to be sliced. If ``estimator._pairwise is True``,
        this needs to be a square array-like or sparse matrix.

    y : array-like, sparse matrix or iterable
        Targets to be sliced.

    indices : array of int
        Rows to select from X and y.
        If ``estimator._pairwise is True`` and ``train_indices is None``
        then ``indices`` will also be used to slice columns.

    train_indices : array of int or None, default=None
        If ``estimator._pairwise is True`` and ``train_indices is not None``,
        then ``train_indices`` will be use to slice the columns of X.

    Returns
    -------
    X_sliced : array-like, sparse matrix or list
        Sliced data.

    y_sliced : array-like, sparse matrix or list
        Sliced targets.

    """
    if getattr(estimator, "_pairwise", False):
        if not hasattr(X, "shape"):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        # X is a precomputed square kernel matrix
        if X.shape[0] != X.shape[1]:
            raise ValueError("X should be a square kernel matrix")
        if train_indices is None:
            X_subset = X[np.ix_(indices, indices)]
        else:
            X_subset = X[np.ix_(indices, train_indices)]
    else:
        X_subset = safe_indexing(X, indices)

    if y is not None:
        y_subset = safe_indexing(y, indices)
    else:
        y_subset = None

    return X_subset, y_subset


class _Router:
    """Matches destinations and props according to a specification

    Parameters
    ----------
    routing : dict
        Each key should be a string naming the input property, or '*'.
        Each value should be a tuple ``(destination, dest_prop)``
        or a list thereof. Here ``destination`` is a pattern to match
        destinations as specified by the meta-estimator using this router.
        Patterns are matched with :mod:`fnmatch`.  ``dest_prop`` is
        an identifier for the prop to be passed to the destination with,
        or '*' to match names with the input prop.
    dest_regex : str
        Used to help validate the destinations

    Examples
    --------
    >>> props = {'foo': [1, 2, 3], 'bar': [4, 5, 6]}
    >>> dests = ['dest1', 'dest2']
    >>> r = _Router({'*': [('*', '*')]})
    >>> (dest1_props, dest2_props), remainder = r(props, dests)
    >>> sorted(dest1_props.items())
    [('bar', [4, 5, 6]), ('foo', [1, 2, 3])]
    >>> sorted(dest2_props.items())
    [('bar', [4, 5, 6]), ('foo', [1, 2, 3])]
    >>> sorted(remainder)
    []

    >>> r = _Router({'foo': [('*', '*')]})
    >>> (dest1_props, dest2_props), remainder = r(props, dests)
    >>> sorted(dest1_props.items())
    [('foo', [1, 2, 3])]
    >>> sorted(dest2_props.items())
    [('foo', [1, 2, 3])]
    >>> sorted(remainder)
    ['bar']

    >>> r = _Router({'foo': [('d*1', '*')]})
    >>> (dest1_props, dest2_props), remainder = r(props, dests)
    >>> sorted(dest1_props.items())
    [('foo', [1, 2, 3])]
    >>> sorted(dest2_props.items())
    []
    >>> sorted(remainder)
    ['bar']

    >>> r = _Router({'foo': [('d*1', 'goo')]})
    >>> (dest1_props, dest2_props), remainder = r(props, dests)
    >>> sorted(dest1_props.items())
    [('goo', [1, 2, 3])]
    >>> sorted(dest2_props.items())
    []
    >>> sorted(remainder)
    ['bar']
    """

    def __init__(self, routing, dest_regex='.*'):
        routing = {prop: dests if isinstance(dests, list) else [dests]
                   for prop, dests in routing.items()}
        self._validate(routing, dest_regex)
        self.routing = routing

    @staticmethod
    def _validate(routing, dest_regex):
        for prop, dests in routing.items():
            for dest_tuple in dests:
                if not isinstance(dest_tuple, tuple):
                    raise ValueError('Each routing destination should be '
                                     '(destination, dest_prop)')
                dest, dest_prop = dest_tuple
                if not isinstance(dest, str):
                    raise ValueError('Expected string for routing '
                                     'destination, got %r' % dest)
                if not isinstance(dest_prop, str) \
                   or not re.match(r'\*$|^[A-Za-z_]\w*$', dest_prop):
                    raise ValueError('Expected identifier or "*" for '
                                     'routing destination prop, got %r'
                                     % dest_prop)
                if not re.match(dest_regex, dest):
                    raise ValueError('Routing destination was not of expected '
                                     'form. %r does not match %r'
                                     % (dest, dest_regex))

            if not isinstance(prop, str):
                raise ValueError('Expected string for prop name, '
                                 'got %r' % prop)

    def __call__(self, props, dests):
        out = {dest: {} for dest in dests}
        remainder = set()
        for prop, val in props.items():
            if val is None:
                continue

            def _set(routes):
                n_matched = 0
                for dest_pattern, dest_prop in routes:
                    if dest_prop == '*':
                        dest_prop = prop
                    for dest in fnmatch.filter(dests, dest_pattern):
                        n_matched += 1
                        if dest_prop in out[dest]:
                            raise ValueError('Already have a prop %r '
                                             'for destination %r'
                                             % (dest_prop, dest))
                        out[dest][dest_prop] = val
                return n_matched

            n_matched = _set(self.routing.get(prop, []))
            n_matched += _set(self.routing.get('*', []))
            if not n_matched:
                remainder.add(prop)

        return [out[dest] for dest in dests], remainder


def check_routing(routing, default=None, dest_regex='.*'):
    if routing is None:
        if default is None:
            raise ValueError('Routing must be specified')
        routing = default
    if not callable(routing):
        return _Router(routing, dest_regex)
    return routing
