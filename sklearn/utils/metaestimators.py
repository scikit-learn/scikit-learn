"""Utilities for meta-estimators"""
# Author: Joel Nothman
#         Andreas Mueller
# License: BSD
from typing import List, Any
import warnings

from abc import ABCMeta, abstractmethod
from operator import attrgetter
import numpy as np
from contextlib import suppress

from ..utils import _safe_indexing
from ..utils._tags import _safe_tags
from ..base import BaseEstimator
from ._available_if import available_if, _AvailableIfDescriptor

__all__ = ["available_if", "if_delegate_has_method"]


class _BaseComposition(BaseEstimator, metaclass=ABCMeta):
    """Handles parameter management for classifiers composed of named estimators."""

    steps: List[Any]

    @abstractmethod
    def __init__(self):
        pass

    def _get_params(self, attr, deep=True):
        out = super().get_params(deep=deep)
        if not deep:
            return out

        estimators = getattr(self, attr)
        try:
            out.update(estimators)
        except (TypeError, ValueError):
            # Ignore TypeError for cases where estimators is not a list of
            # (name, estimator) and ignore ValueError when the list is not
            # formatted correctly. This is to prevent errors when calling
            # `set_params`. `BaseEstimator.set_params` calls `get_params` which
            # can error for invalid values for `estimators`.
            return out

        for name, estimator in estimators:
            if hasattr(estimator, "get_params"):
                for key, value in estimator.get_params(deep=True).items():
                    out["%s__%s" % (name, key)] = value
        return out

    def _set_params(self, attr, **params):
        # Ensure strict ordering of parameter setting:
        # 1. All steps
        if attr in params:
            setattr(self, attr, params.pop(attr))
        # 2. Replace items with estimators in params
        items = getattr(self, attr)
        if isinstance(items, list) and items:
            # Get item names used to identify valid names in params
            # `zip` raises a TypeError when `items` does not contains
            # elements of length 2
            with suppress(TypeError):
                item_names, _ = zip(*items)
                for name in list(params.keys()):
                    if "__" not in name and name in item_names:
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
            raise ValueError("Names provided are not unique: {0!r}".format(list(names)))
        invalid_names = set(names).intersection(self.get_params(deep=False))
        if invalid_names:
            raise ValueError(
                "Estimator names conflict with constructor arguments: {0!r}".format(
                    sorted(invalid_names)
                )
            )
        invalid_names = [name for name in names if "__" in name]
        if invalid_names:
            raise ValueError(
                "Estimator names must not contain __: got {0!r}".format(invalid_names)
            )


# TODO(1.3) remove
class _IffHasAttrDescriptor(_AvailableIfDescriptor):
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
        super().__init__(fn, self._check, attribute_name)
        self.delegate_names = delegate_names

    def _check(self, obj):
        warnings.warn(
            "if_delegate_has_method was deprecated in version 1.1 and will be "
            "removed in version 1.3. Use available_if instead.",
            FutureWarning,
        )

        delegate = None
        for delegate_name in self.delegate_names:
            try:
                delegate = attrgetter(delegate_name)(obj)
                break
            except AttributeError:
                continue

        if delegate is None:
            return False
        # raise original AttributeError
        getattr(delegate, self.attribute_name)

        return True


# TODO(1.3) remove
def if_delegate_has_method(delegate):
    """Create a decorator for methods that are delegated to a sub-estimator.

    .. deprecated:: 1.3
        `if_delegate_has_method` is deprecated in version 1.1 and will be removed in
        version 1.3. Use `available_if` instead.

    This enables ducktyping by hasattr returning True according to the
    sub-estimator.

    Parameters
    ----------
    delegate : str, list of str or tuple of str
        Name of the sub-estimator that can be accessed as an attribute of the
        base object. If a list or a tuple of names are provided, the first
        sub-estimator that is an attribute of the base object will be used.

    Returns
    -------
    callable
        Callable makes the decorated method available if the delegate
        has a method with the same name as the decorated method.
    """
    if isinstance(delegate, list):
        delegate = tuple(delegate)
    if not isinstance(delegate, tuple):
        delegate = (delegate,)

    return lambda fn: _IffHasAttrDescriptor(fn, delegate, attribute_name=fn.__name__)


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
    if _safe_tags(estimator, key="pairwise"):
        if not hasattr(X, "shape"):
            raise ValueError(
                "Precomputed kernels or affinity matrices have "
                "to be passed as arrays or sparse matrices."
            )
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
