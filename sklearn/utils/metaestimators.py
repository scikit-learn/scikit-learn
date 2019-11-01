"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD

from abc import ABCMeta, abstractmethod
from operator import attrgetter
from functools import update_wrapper
from collections import defaultdict
import copy
import six
import numpy as np

from ..utils import _safe_indexing
from ..base import BaseEstimator

__all__ = ['if_delegate_has_method']


class _BaseComposition(BaseEstimator, metaclass=ABCMeta):
    """Handles parameter management for classifiers composed of named estimators.
    """
    @abstractmethod
    def __init__(self):
        pass

    def _get_params(self, attr, deep=True):
        out = super().get_params(deep=deep)
        if not deep:
            return out
        estimators = getattr(self, attr)
        out.update(estimators)
        for name, estimator in estimators:
            if hasattr(estimator, 'get_params'):
                for key, value in estimator.get_params(deep=True).items():
                    out['%s__%s' % (name, key)] = value
        return out

    def _set_params(self, attr, **params):
        # Ensure strict ordering of parameter setting:
        # 1. All steps
        if attr in params:
            setattr(self, attr, params.pop(attr))
        # 2. Step replacement
        items = getattr(self, attr)
        names = []
        if items:
            names, _ = zip(*items)
        for name in list(params.keys()):
            if '__' not in name and name in names:
                self._replace_estimator(attr, name, params.pop(name))
        # 3. Step parameters and other initialisation arguments
        super().set_params(**params)
        return self

    def _replace_estimator(self, attr, name, new_val):
        # assumes `name` is a valid estimator name
        new_estimators = list(getattr(self, attr))
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


class _IffHasAttrDescriptor:
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

    Labels y will always be indexed only along the first axis.

    Parameters
    ----------
    estimator : object
        Estimator to determine whether we should slice only rows or rows and
        columns.

    X : array-like, sparse matrix or iterable
        Data to be indexed. If ``estimator._pairwise is True``,
        this needs to be a square array-like or sparse matrix.

    y : array-like, sparse matrix or iterable
        Targets to be indexed.

    indices : array of int
        Rows to select from X and y.
        If ``estimator._pairwise is True`` and ``train_indices is None``
        then ``indices`` will also be used to slice columns.

    train_indices : array of int or None, default=None
        If ``estimator._pairwise is True`` and ``train_indices is not None``,
        then ``train_indices`` will be use to slice the columns of X.

    Returns
    -------
    X_subset : array-like, sparse matrix or list
        Indexed data.

    y_subset : array-like, sparse matrix or list
        Indexed targets.

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
        X_subset = _safe_indexing(X, indices)

    if y is not None:
        y_subset = _safe_indexing(y, indices)
    else:
        y_subset = None

    return X_subset, y_subset


class _NonePolicy:
    def apply(self, props, unused):
        return {}

    def __repr__(self):
        return 'None'


class _AllPolicy:
    def apply(self, props, unused):
        unused.clear()
        return props

    def update(self, other):
        pass

    def __repr__(self):
        return repr('*')


class _ExcludePolicy:
    def __init__(self, exclusions):
        # exclusions is a set of strings
        self.exclusions = exclusions

    def apply(self, props, unused):
        out = {k: v for k, v in props.items()
               if k not in self.exclusions}
        unused.intersection_update(self.exclusions)
        return out

    def update(self, other):
        self.exclusions.update(other)

    def __repr__(self):
        return repr(['-' + k for k in sorted(self.exclusions)])


class _IncludePolicy:
    def __init__(self, inclusions):
        # inclusions maps {tgt: src}
        self.inclusions = inclusions

    def apply(self, props, unused):
        unused.difference_update(self.inclusions.values())
        return {tgt: props[src]
                for tgt, src in self.inclusions.items()
                if src in props}

    def update(self, other):
        intersection = set(self.inclusions) & set(other.inclusions)
        if intersection:
            raise ValueError('Target property names {!r} are used multiple '
                             'times in the routing policy for the same '
                             'destination'.format(sorted(intersection)))
        self.inclusions.update(other.inclusions)

    def __repr__(self):
        return '{%s}' % ', '.join('{!r}: {!r}'.format(tgt, src)
                                  for tgt, src
                                  in self.inclusions.items())


class _Router:
    """Matches sample props to destinations according to a routing policy
    Parameters
    ----------
    routing : dict {dest: props} for dest str, props {str, list, dict}
        User-defined routing policy.
        Maps each destination string to the properties that should be provided
        to that destination. Props may be:
        - '*': provide all properties
        - a list of property names to include
        - a list of property names, each prefixed by '-', to exclude; all
          others will be provided
        - a single name to include or exclude
        - a dict mapping each target property name to its source property name
    dests : list of {str, iterable of str}
        The ordered destinations for this router. If a set of strings
        is provided for each entry, any of these strings will route
        parameters to the destination.
        Usually this is fixed for a metaestimator.
    Notes
    -----
    Abstracting away destination names/aliases in this way allows for providing
    syntactic sugar to users (e.g. Pipeline can declare '*' or 'steps' as an
    alias for "provide these props to ``fit`` of all steps).  While this may be
    instead facilitated by some string-based pattern matching, the present
    approach is more explicit and ensures backwards compatibility can be
    maintained.
    """

    def __init__(self, routing, dests):
        # can immediately:
        #     * check that all routes have valid dests
        #     * consolidate routing across aliases, such that each must be
        #       either a blacklist of length 0 or more or a mapping
        alias_to_idx = defaultdict(list)
        for i, aliases in enumerate(dests):
            if isinstance(aliases, six.string_types):
                aliases = [aliases]
            for alias in aliases:
                alias_to_idx[alias].append(i)
        alias_to_idx.default_factory = None

        policies = [None] * len(dests)

        # Sorted so error messages are deterministic
        for dest, props in sorted(routing.items()):
            if props == '*':
                policy = _AllPolicy()
            else:
                if isinstance(props, six.string_types):
                    props = [props]

                if isinstance(props, dict):
                    policy = _IncludePolicy(props)
                else:
                    minuses = [prop[:1] == '-' for prop in props]
                    if all(minuses):
                        policy = _ExcludePolicy({prop[1:] for prop in props})
                    elif any(minuses):
                        raise ValueError('Routing props should either all '
                                         'start with "-" or none should start '
                                         'with "-". Got a mix for %r' % dest)
                    else:
                        policy = _IncludePolicy({prop: prop for prop in props})

            # raises KeyError if unknown dest
            for idx in alias_to_idx[dest]:
                if policies[idx] is None:
                    policies[idx] = copy.deepcopy(policy)
                else:
                    if type(policies[idx]) is not type(policy):
                        raise ValueError('When handling routing for '
                                         'destination {!r}, found a mix of '
                                         'inclusion, exclusion and pass all '
                                         'policies.'.format(dest))
                    policies[idx].update(policy)

        self.policies = [_NonePolicy() if policy is None else policy
                         for policy in policies]

    def __call__(self, props):
        """Apply the routing policy to the given sample props
        Parameters
        ----------
        props : dict
        Returns
        -------
        dest_props : list of dicts
            Props to be passed to each destination declared in the constructor.
        unused : set of str
            Names of props that were not routed anywhere.
        """
        unused = set(props)
        out = [policy.apply(props, unused)
               for policy in self.policies]
        return out, unused


def check_routing(routing, dests, default=None):
    """Validates a prop_routing parameter and returns a router
    A router is a function which takes a dict of sample properties and returns
    a tuple ``(dest_props, unused)``, defined by:
        dest_props : list of dicts
            Props to be passed to each destination declared by ``dests``.
        unused : set of str
            Names of props that were not routed anywhere.
    Parameters
    ----------
    routing : dict of {dest: props}, callable or None if default is given
        User-defined routing policy.  A callable is returned unchanged.
        Maps each destination string to the properties that should be provided
        to that destination. Props may be:
        - ``'*'``: provide all properties
        - a list of property names to include
        - a list of property names, each prefixed by ``-``, to exclude; all
          others will be provided
        - a single name to include or exclude (prefixed by ``-``)
        - a dict mapping each target property name to its source property name
    dests : list of {str, iterable of str}
        The ordered destination names for the router.  If a set of strings is
        provided for each entry, any of these aliases will route parameters to
        the destination.  The same string may appear in multiple dests entries.
        This should be fixed for each metaestimator.
    default : dict or callable, optional
        This replaces ``routing`` where routing is None.
        This should be fixed for each metaestimator.
    Returns
    -------
    router : callable with signature ``props -> (dest_props, unused)``
    Examples
    --------
    >>> from sklearn.utils.metaestimators import check_routing
    >>> props = {'foo': [1, 2], 'bar': [3, 4]}
    >>> dests = [['d1', 'all'], ['d2', 'all']]
    >>> # route 'foo' to d1 and 'bar' to d2
    >>> router = check_routing({'d1': 'foo', 'd2': 'bar'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'bar': [3, 4]}
    >>> list(unused)
    []
    >>> # rename 'bar' to 'baz'
    >>> router = check_routing({'d1': 'foo', 'd2': {'baz': 'bar'}}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'baz': [3, 4]}
    >>> list(unused)
    []
    >>> # d2 takes all but foo
    >>> router = check_routing({'d1': 'foo', 'd2': '-foo'}, dests)
    >>> router(props)[0]
    [{'foo': [1, 2]}, {'bar': [3, 4]}]
    >>> # d1 takes all; d2 takes none
    >>> router = check_routing({'d1': '*'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> sorted(d1_props.items())
    [('bar', [3, 4]), ('foo', [1, 2])]
    >>> d2_props
    {}
    >>> # the 'all' alias distributes to both dests
    >>> router = check_routing({'all': 'foo', 'd1': 'bar'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> sorted(d1_props.items())
    [('bar', [3, 4]), ('foo', [1, 2])]
    >>> d2_props
    {'foo': [1, 2]}
    >>> # an unused prop
    >>> router = check_routing({'all': 'foo'}, dests)
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'foo': [1, 2]}
    >>> list(unused)
    ['bar']
    >>> # a default: both get foo
    >>> router = check_routing(None, dests, default={'all': 'foo'})
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {'foo': [1, 2]}
    >>> # an overridden default: only d1 gets foo
    >>> router = check_routing({'d1': 'foo'}, dests, default={'all': 'foo'})
    >>> (d1_props, d2_props), unused = router(props)
    >>> d1_props
    {'foo': [1, 2]}
    >>> d2_props
    {}
    """
    if routing is None:
        if default is None:
            raise ValueError('Routing must be specified')
        routing = default
    if not callable(routing):
        return _Router(routing, dests)
    return routing
