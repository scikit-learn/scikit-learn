"""
The :mod:`sklearn.model_selection._split` module includes classes and
functions to split the data based on a preset strategy.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>,
#         Olivier Grisel <olivier.grisel@ensta.org>
#         Raghav RV <rvraghav93@gmail.com>
# License: BSD 3 clause

from collections.abc import Iterable
import warnings
from itertools import chain, combinations
from math import ceil, floor
import numbers
from abc import ABCMeta, abstractmethod
from inspect import signature

import numpy as np

from ..utils import indexable, check_random_state, safe_indexing
from ..utils import _approximate_mode
from ..utils.validation import _num_samples, column_or_1d
from ..utils.validation import check_array
from ..utils.multiclass import type_of_target
from ..utils.fixes import comb
from ..base import _pprint

__all__ = ['BaseCrossValidator',
           'KFold',
           'GroupKFold',
           'LeaveOneGroupOut',
           'LeaveOneOut',
           'LeavePGroupsOut',
           'LeavePOut',
           'RepeatedStratifiedKFold',
           'RepeatedKFold',
           'ShuffleSplit',
           'GroupShuffleSplit',
           'StratifiedKFold',
           'StratifiedShuffleSplit',
           'PredefinedSplit',
           'train_test_split',
           'check_cv']


class BaseCrossValidator(metaclass=ABCMeta):
    """Base class for all cross-validators

    Implementations must define `_iter_test_masks` or `_iter_test_indices`.
    """

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, of length n_samples
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        X, y, groups = indexable(X, y, groups)
        indices = np.arange(_num_samples(X))
        for test_index in self._iter_test_masks(X, y, groups):
            train_index = indices[np.logical_not(test_index)]
            test_index = indices[test_index]
            yield train_index, test_index

    # Since subclasses must implement either _iter_test_masks or
    # _iter_test_indices, neither can be abstract.
    def _iter_test_masks(self, X=None, y=None, groups=None):
        """Generates boolean masks corresponding to test sets.

        By default, delegates to _iter_test_indices(X, y, groups)
        """
        for test_index in self._iter_test_indices(X, y, groups):
            test_mask = np.zeros(_num_samples(X), dtype=np.bool)
            test_mask[test_index] = True
            yield test_mask

    def _iter_test_indices(self, X=None, y=None, groups=None):
        """Generates integer indices corresponding to test sets."""
        raise NotImplementedError

    @abstractmethod
    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator"""

    def __repr__(self):
        return _build_repr(self)


class LeaveOneOut(BaseCrossValidator):
    """Leave-One-Out cross-validator

    Provides train/test indices to split data in train/test sets. Each
    sample is used once as a test set (singleton) while the remaining
    samples form the training set.

    Note: ``LeaveOneOut()`` is equivalent to ``KFold(n_splits=n)`` and
    ``LeavePOut(p=1)`` where ``n`` is the number of samples.

    Due to the high number of test sets (which is the same as the
    number of samples) this cross-validation method can be very costly.
    For large datasets one should favor :class:`KFold`, :class:`ShuffleSplit`
    or :class:`StratifiedKFold`.

    Read more in the :ref:`User Guide <cross_validation>`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import LeaveOneOut
    >>> X = np.array([[1, 2], [3, 4]])
    >>> y = np.array([1, 2])
    >>> loo = LeaveOneOut()
    >>> loo.get_n_splits(X)
    2
    >>> print(loo)
    LeaveOneOut()
    >>> for train_index, test_index in loo.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    ...    print(X_train, X_test, y_train, y_test)
    TRAIN: [1] TEST: [0]
    [[3 4]] [[1 2]] [2] [1]
    TRAIN: [0] TEST: [1]
    [[1 2]] [[3 4]] [1] [2]

    See also
    --------
    LeaveOneGroupOut
        For splitting the data according to explicit, domain-specific
        stratification of the dataset.

    GroupKFold: K-fold iterator variant with non-overlapping groups.
    """

    def _iter_test_indices(self, X, y=None, groups=None):
        n_samples = _num_samples(X)
        if n_samples <= 1:
            raise ValueError(
                'Cannot perform LeaveOneOut with n_samples={}.'.format(
                    n_samples)
            )
        return range(n_samples)

    def get_n_splits(self, X, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if X is None:
            raise ValueError("The 'X' parameter should not be None.")
        return _num_samples(X)


class LeavePOut(BaseCrossValidator):
    """Leave-P-Out cross-validator

    Provides train/test indices to split data in train/test sets. This results
    in testing on all distinct samples of size p, while the remaining n - p
    samples form the training set in each iteration.

    Note: ``LeavePOut(p)`` is NOT equivalent to
    ``KFold(n_splits=n_samples // p)`` which creates non-overlapping test sets.

    Due to the high number of iterations which grows combinatorically with the
    number of samples this cross-validation method can be very costly. For
    large datasets one should favor :class:`KFold`, :class:`StratifiedKFold`
    or :class:`ShuffleSplit`.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    p : int
        Size of the test sets. Must be strictly less than the number of
        samples.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import LeavePOut
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 3, 4])
    >>> lpo = LeavePOut(2)
    >>> lpo.get_n_splits(X)
    6
    >>> print(lpo)
    LeavePOut(p=2)
    >>> for train_index, test_index in lpo.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [2 3] TEST: [0 1]
    TRAIN: [1 3] TEST: [0 2]
    TRAIN: [1 2] TEST: [0 3]
    TRAIN: [0 3] TEST: [1 2]
    TRAIN: [0 2] TEST: [1 3]
    TRAIN: [0 1] TEST: [2 3]
    """

    def __init__(self, p):
        self.p = p

    def _iter_test_indices(self, X, y=None, groups=None):
        n_samples = _num_samples(X)
        if n_samples <= self.p:
            raise ValueError(
                'p={} must be strictly less than the number of '
                'samples={}'.format(self.p, n_samples)
            )
        for combination in combinations(range(n_samples), self.p):
            yield np.array(combination)

    def get_n_splits(self, X, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.
        """
        if X is None:
            raise ValueError("The 'X' parameter should not be None.")
        return int(comb(_num_samples(X), self.p, exact=True))


class _BaseKFold(BaseCrossValidator, metaclass=ABCMeta):
    """Base class for KFold, GroupKFold, and StratifiedKFold"""

    @abstractmethod
    def __init__(self, n_splits, shuffle, random_state):
        if not isinstance(n_splits, numbers.Integral):
            raise ValueError('The number of folds must be of Integral type. '
                             '%s of type %s was passed.'
                             % (n_splits, type(n_splits)))
        n_splits = int(n_splits)

        if n_splits <= 1:
            raise ValueError(
                "k-fold cross-validation requires at least one"
                " train/test split by setting n_splits=2 or more,"
                " got n_splits={0}.".format(n_splits))

        if not isinstance(shuffle, bool):
            raise TypeError("shuffle must be True or False;"
                            " got {0}".format(shuffle))

        self.n_splits = n_splits
        self.shuffle = shuffle
        self.random_state = random_state

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        X, y, groups = indexable(X, y, groups)
        n_samples = _num_samples(X)
        if self.n_splits > n_samples:
            raise ValueError(
                ("Cannot have number of splits n_splits={0} greater"
                 " than the number of samples: n_samples={1}.")
                .format(self.n_splits, n_samples))

        for train, test in super().split(X, y, groups):
            yield train, test

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return self.n_splits


class KFold(_BaseKFold):
    """K-Folds cross-validator

    Provides train/test indices to split data in train/test sets. Split
    dataset into k consecutive folds (without shuffling by default).

    Each fold is then used once as a validation while the k - 1 remaining
    folds form the training set.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

        .. versionchanged:: 0.22
            ``n_splits`` default value changed from 3 to 5.

    shuffle : boolean, optional
        Whether to shuffle the data before splitting into batches.

    random_state : int, RandomState instance or None, optional, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used when ``shuffle`` == True.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import KFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4])
    >>> kf = KFold(n_splits=2)
    >>> kf.get_n_splits(X)
    2
    >>> print(kf)
    KFold(n_splits=2, random_state=None, shuffle=False)
    >>> for train_index, test_index in kf.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [2 3] TEST: [0 1]
    TRAIN: [0 1] TEST: [2 3]

    Notes
    -----
    The first ``n_samples % n_splits`` folds have size
    ``n_samples // n_splits + 1``, other folds have size
    ``n_samples // n_splits``, where ``n_samples`` is the number of samples.

    Randomized CV splitters may return different results for each call of
    split. You can make the results identical by setting ``random_state``
    to an integer.

    See also
    --------
    StratifiedKFold
        Takes group information into account to avoid building folds with
        imbalanced class distributions (for binary or multiclass
        classification tasks).

    GroupKFold: K-fold iterator variant with non-overlapping groups.

    RepeatedKFold: Repeats K-Fold n times.
    """

    def __init__(self, n_splits=5, shuffle=False,
                 random_state=None):
        super().__init__(n_splits, shuffle, random_state)

    def _iter_test_indices(self, X, y=None, groups=None):
        n_samples = _num_samples(X)
        indices = np.arange(n_samples)
        if self.shuffle:
            check_random_state(self.random_state).shuffle(indices)

        n_splits = self.n_splits
        fold_sizes = np.full(n_splits, n_samples // n_splits, dtype=np.int)
        fold_sizes[:n_samples % n_splits] += 1
        current = 0
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            yield indices[start:stop]
            current = stop


class GroupKFold(_BaseKFold):
    """K-fold iterator variant with non-overlapping groups.

    The same group will not appear in two different folds (the number of
    distinct groups has to be at least equal to the number of folds).

    The folds are approximately balanced in the sense that the number of
    distinct groups is approximately the same in each fold.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

        .. versionchanged:: 0.22
            ``n_splits`` default value changed from 3 to 5.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import GroupKFold
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 3, 4])
    >>> groups = np.array([0, 0, 2, 2])
    >>> group_kfold = GroupKFold(n_splits=2)
    >>> group_kfold.get_n_splits(X, y, groups)
    2
    >>> print(group_kfold)
    GroupKFold(n_splits=2)
    >>> for train_index, test_index in group_kfold.split(X, y, groups):
    ...     print("TRAIN:", train_index, "TEST:", test_index)
    ...     X_train, X_test = X[train_index], X[test_index]
    ...     y_train, y_test = y[train_index], y[test_index]
    ...     print(X_train, X_test, y_train, y_test)
    ...
    TRAIN: [0 1] TEST: [2 3]
    [[1 2]
     [3 4]] [[5 6]
     [7 8]] [1 2] [3 4]
    TRAIN: [2 3] TEST: [0 1]
    [[5 6]
     [7 8]] [[1 2]
     [3 4]] [3 4] [1 2]

    See also
    --------
    LeaveOneGroupOut
        For splitting the data according to explicit domain-specific
        stratification of the dataset.
    """
    def __init__(self, n_splits=5):
        super().__init__(n_splits, shuffle=False, random_state=None)

    def _iter_test_indices(self, X, y, groups):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, ensure_2d=False, dtype=None)

        unique_groups, groups = np.unique(groups, return_inverse=True)
        n_groups = len(unique_groups)

        if self.n_splits > n_groups:
            raise ValueError("Cannot have number of splits n_splits=%d greater"
                             " than the number of groups: %d."
                             % (self.n_splits, n_groups))

        # Weight groups by their number of occurrences
        n_samples_per_group = np.bincount(groups)

        # Distribute the most frequent groups first
        indices = np.argsort(n_samples_per_group)[::-1]
        n_samples_per_group = n_samples_per_group[indices]

        # Total weight of each fold
        n_samples_per_fold = np.zeros(self.n_splits)

        # Mapping from group index to fold index
        group_to_fold = np.zeros(len(unique_groups))

        # Distribute samples by adding the largest weight to the lightest fold
        for group_index, weight in enumerate(n_samples_per_group):
            lightest_fold = np.argmin(n_samples_per_fold)
            n_samples_per_fold[lightest_fold] += weight
            group_to_fold[indices[group_index]] = lightest_fold

        indices = group_to_fold[groups]

        for f in range(self.n_splits):
            yield np.where(indices == f)[0]

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,), optional
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        return super().split(X, y, groups)


class StratifiedKFold(_BaseKFold):
    """Stratified K-Folds cross-validator

    Provides train/test indices to split data in train/test sets.

    This cross-validation object is a variation of KFold that returns
    stratified folds. The folds are made by preserving the percentage of
    samples for each class.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

        .. versionchanged:: 0.22
            ``n_splits`` default value changed from 3 to 5.

    shuffle : boolean, optional
        Whether to shuffle each class's samples before splitting into batches.

    random_state : int, RandomState instance or None, optional, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`. Used when ``shuffle`` == True.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import StratifiedKFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> skf = StratifiedKFold(n_splits=2)
    >>> skf.get_n_splits(X, y)
    2
    >>> print(skf)
    StratifiedKFold(n_splits=2, random_state=None, shuffle=False)
    >>> for train_index, test_index in skf.split(X, y):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 3] TEST: [0 2]
    TRAIN: [0 2] TEST: [1 3]

    Notes
    -----
    The implementation is designed to:

    * Generate test sets such that all contain the same distribution of
      classes, or as close as possible.
    * Be invariant to class label: relabelling ``y = ["Happy", "Sad"]`` to
      ``y = [1, 0]`` should not change the indices generated.
    * Preserve order dependencies in the dataset ordering, when
      ``shuffle=False``: all samples from class k in some test set were
      contiguous in y, or separated in y by samples from classes other than k.
    * Generate test sets where the smallest and largest differ by at most one
      sample.

    .. versionchanged:: 0.22
        The previous implementation did not follow the last constraint.

    See also
    --------
    RepeatedStratifiedKFold: Repeats Stratified K-Fold n times.
    """

    def __init__(self, n_splits=5, shuffle=False, random_state=None):
        super().__init__(n_splits, shuffle, random_state)

    def _make_test_folds(self, X, y=None):
        rng = check_random_state(self.random_state)
        y = np.asarray(y)
        type_of_target_y = type_of_target(y)
        allowed_target_types = ('binary', 'multiclass')
        if type_of_target_y not in allowed_target_types:
            raise ValueError(
                'Supported target types are: {}. Got {!r} instead.'.format(
                    allowed_target_types, type_of_target_y))

        y = column_or_1d(y)

        _, y_idx, y_inv = np.unique(y, return_index=True, return_inverse=True)
        # y_inv encodes y according to lexicographic order. We invert y_idx to
        # map the classes so that they are encoded by order of appearance:
        # 0 represents the first label appearing in y, 1 the second, etc.
        _, class_perm = np.unique(y_idx, return_inverse=True)
        y_encoded = class_perm[y_inv]

        n_classes = len(y_idx)
        y_counts = np.bincount(y_encoded)
        min_groups = np.min(y_counts)
        if np.all(self.n_splits > y_counts):
            raise ValueError("n_splits=%d cannot be greater than the"
                             " number of members in each class."
                             % (self.n_splits))
        if self.n_splits > min_groups:
            warnings.warn(("The least populated class in y has only %d"
                           " members, which is less than n_splits=%d."
                           % (min_groups, self.n_splits)), UserWarning)

        # Determine the optimal number of samples from each class in each fold,
        # using round robin over the sorted y. (This can be done direct from
        # counts, but that code is unreadable.)
        y_order = np.sort(y_encoded)
        allocation = np.asarray(
            [np.bincount(y_order[i::self.n_splits], minlength=n_classes)
             for i in range(self.n_splits)])

        # To maintain the data order dependencies as best as possible within
        # the stratification constraint, we assign samples from each class in
        # blocks (and then mess that up when shuffle=True).
        test_folds = np.empty(len(y), dtype='i')
        for k in range(n_classes):
            # since the kth column of allocation stores the number of samples
            # of class k in each test set, this generates blocks of fold
            # indices corresponding to the allocation for class k.
            folds_for_class = np.arange(self.n_splits).repeat(allocation[:, k])
            if self.shuffle:
                rng.shuffle(folds_for_class)
            test_folds[y_encoded == k] = folds_for_class
        return test_folds

    def _iter_test_masks(self, X, y=None, groups=None):
        test_folds = self._make_test_folds(X, y)
        for i in range(self.n_splits):
            yield test_folds == i

    def split(self, X, y, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

            Note that providing ``y`` is sufficient to generate the splits and
            hence ``np.zeros(n_samples)`` may be used as a placeholder for
            ``X`` instead of actual training data.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.
            Stratification is done based on the y labels.

        groups : object
            Always ignored, exists for compatibility.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.

        Notes
        -----
        Randomized CV splitters may return different results for each call of
        split. You can make the results identical by setting ``random_state``
        to an integer.
        """
        y = check_array(y, ensure_2d=False, dtype=None)
        return super().split(X, y, groups)


class TimeSeriesSplit(_BaseKFold):
    """Time Series cross-validator

    Provides train/test indices to split time series data samples
    that are observed at fixed time intervals, in train/test sets.
    In each split, test indices must be higher than before, and thus shuffling
    in cross validator is inappropriate.

    This cross-validation object is a variation of :class:`KFold`.
    In the kth split, it returns first k folds as train set and the
    (k+1)th fold as test set.

    Note that unlike standard cross-validation methods, successive
    training sets are supersets of those that come before them.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of splits. Must be at least 2.

        .. versionchanged:: 0.22
            ``n_splits`` default value changed from 3 to 5.

    max_train_size : int, optional
        Maximum size for a single training set.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import TimeSeriesSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4, 5, 6])
    >>> tscv = TimeSeriesSplit()
    >>> print(tscv)
    TimeSeriesSplit(max_train_size=None, n_splits=5)
    >>> for train_index, test_index in tscv.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [0] TEST: [1]
    TRAIN: [0 1] TEST: [2]
    TRAIN: [0 1 2] TEST: [3]
    TRAIN: [0 1 2 3] TEST: [4]
    TRAIN: [0 1 2 3 4] TEST: [5]

    Notes
    -----
    The training set has size ``i * n_samples // (n_splits + 1)
    + n_samples % (n_splits + 1)`` in the ``i``th split,
    with a test set of size ``n_samples//(n_splits + 1)``,
    where ``n_samples`` is the number of samples.
    """
    def __init__(self, n_splits=5, max_train_size=None):
        super().__init__(n_splits, shuffle=False, random_state=None)
        self.max_train_size = max_train_size

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Always ignored, exists for compatibility.

        groups : array-like, with shape (n_samples,)
            Always ignored, exists for compatibility.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        X, y, groups = indexable(X, y, groups)
        n_samples = _num_samples(X)
        n_splits = self.n_splits
        n_folds = n_splits + 1
        if n_folds > n_samples:
            raise ValueError(
                ("Cannot have number of folds ={0} greater"
                 " than the number of samples: {1}.").format(n_folds,
                                                             n_samples))
        indices = np.arange(n_samples)
        test_size = (n_samples // n_folds)
        test_starts = range(test_size + n_samples % n_folds,
                            n_samples, test_size)
        for test_start in test_starts:
            if self.max_train_size and self.max_train_size < test_start:
                yield (indices[test_start - self.max_train_size:test_start],
                       indices[test_start:test_start + test_size])
            else:
                yield (indices[:test_start],
                       indices[test_start:test_start + test_size])


class LeaveOneGroupOut(BaseCrossValidator):
    """Leave One Group Out cross-validator

    Provides train/test indices to split data according to a third-party
    provided group. This group information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the groups could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    Read more in the :ref:`User Guide <cross_validation>`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import LeaveOneGroupOut
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 1, 2])
    >>> groups = np.array([1, 1, 2, 2])
    >>> logo = LeaveOneGroupOut()
    >>> logo.get_n_splits(X, y, groups)
    2
    >>> logo.get_n_splits(groups=groups)  # 'groups' is always required
    2
    >>> print(logo)
    LeaveOneGroupOut()
    >>> for train_index, test_index in logo.split(X, y, groups):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    ...    print(X_train, X_test, y_train, y_test)
    TRAIN: [2 3] TEST: [0 1]
    [[5 6]
     [7 8]] [[1 2]
     [3 4]] [1 2] [1 2]
    TRAIN: [0 1] TEST: [2 3]
    [[1 2]
     [3 4]] [[5 6]
     [7 8]] [1 2] [1 2]

    """

    def _iter_test_masks(self, X, y, groups):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        # We make a copy of groups to avoid side-effects during iteration
        groups = check_array(groups, copy=True, ensure_2d=False, dtype=None)
        unique_groups = np.unique(groups)
        if len(unique_groups) <= 1:
            raise ValueError(
                "The groups parameter contains fewer than 2 unique groups "
                "(%s). LeaveOneGroupOut expects at least 2." % unique_groups)
        for i in unique_groups:
            yield groups == i

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set. This 'groups' parameter must always be specified to
            calculate the number of splits, though the other parameters can be
            omitted.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, ensure_2d=False, dtype=None)
        return len(np.unique(groups))

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, of length n_samples, optional
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        return super().split(X, y, groups)


class LeavePGroupsOut(BaseCrossValidator):
    """Leave P Group(s) Out cross-validator

    Provides train/test indices to split data according to a third-party
    provided group. This group information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the groups could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    The difference between LeavePGroupsOut and LeaveOneGroupOut is that
    the former builds the test sets with all the samples assigned to
    ``p`` different values of the groups while the latter uses samples
    all assigned the same groups.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_groups : int
        Number of groups (``p``) to leave out in the test split.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import LeavePGroupsOut
    >>> X = np.array([[1, 2], [3, 4], [5, 6]])
    >>> y = np.array([1, 2, 1])
    >>> groups = np.array([1, 2, 3])
    >>> lpgo = LeavePGroupsOut(n_groups=2)
    >>> lpgo.get_n_splits(X, y, groups)
    3
    >>> lpgo.get_n_splits(groups=groups)  # 'groups' is always required
    3
    >>> print(lpgo)
    LeavePGroupsOut(n_groups=2)
    >>> for train_index, test_index in lpgo.split(X, y, groups):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    ...    print(X_train, X_test, y_train, y_test)
    TRAIN: [2] TEST: [0 1]
    [[5 6]] [[1 2]
     [3 4]] [1] [1 2]
    TRAIN: [1] TEST: [0 2]
    [[3 4]] [[1 2]
     [5 6]] [2] [1 1]
    TRAIN: [0] TEST: [1 2]
    [[1 2]] [[3 4]
     [5 6]] [1] [2 1]

    See also
    --------
    GroupKFold: K-fold iterator variant with non-overlapping groups.
    """

    def __init__(self, n_groups):
        self.n_groups = n_groups

    def _iter_test_masks(self, X, y, groups):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, copy=True, ensure_2d=False, dtype=None)
        unique_groups = np.unique(groups)
        if self.n_groups >= len(unique_groups):
            raise ValueError(
                "The groups parameter contains fewer than (or equal to) "
                "n_groups (%d) numbers of unique groups (%s). LeavePGroupsOut "
                "expects that at least n_groups + 1 (%d) unique groups be "
                "present" % (self.n_groups, unique_groups, self.n_groups + 1))
        combi = combinations(range(len(unique_groups)), self.n_groups)
        for indices in combi:
            test_index = np.zeros(_num_samples(X), dtype=np.bool)
            for l in unique_groups[np.array(indices)]:
                test_index[groups == l] = True
            yield test_index

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set. This 'groups' parameter must always be specified to
            calculate the number of splits, though the other parameters can be
            omitted.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, ensure_2d=False, dtype=None)
        return int(comb(len(np.unique(groups)), self.n_groups, exact=True))

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, of length n_samples, optional
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        return super().split(X, y, groups)


class _RepeatedSplits(metaclass=ABCMeta):
    """Repeated splits for an arbitrary randomized CV splitter.

    Repeats splits for cross-validators n times with different randomization
    in each repetition.

    Parameters
    ----------
    cv : callable
        Cross-validator class.

    n_repeats : int, default=10
        Number of times cross-validator needs to be repeated.

    random_state : int, RandomState instance or None, optional, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    **cvargs : additional params
        Constructor parameters for cv. Must not contain random_state
        and shuffle.
    """
    def __init__(self, cv, n_repeats=10, random_state=None, **cvargs):
        if not isinstance(n_repeats, numbers.Integral):
            raise ValueError("Number of repetitions must be of Integral type.")

        if n_repeats <= 0:
            raise ValueError("Number of repetitions must be greater than 0.")

        if any(key in cvargs for key in ('random_state', 'shuffle')):
            raise ValueError(
                "cvargs must not contain random_state or shuffle.")

        self.cv = cv
        self.n_repeats = n_repeats
        self.random_state = random_state
        self.cvargs = cvargs

    def split(self, X, y=None, groups=None):
        """Generates indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, of length n_samples
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        n_repeats = self.n_repeats
        rng = check_random_state(self.random_state)

        for idx in range(n_repeats):
            cv = self.cv(random_state=rng, shuffle=True,
                         **self.cvargs)
            for train_index, test_index in cv.split(X, y, groups):
                yield train_index, test_index

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.
            ``np.zeros(n_samples)`` may be used as a placeholder.

        y : object
            Always ignored, exists for compatibility.
            ``np.zeros(n_samples)`` may be used as a placeholder.

        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        rng = check_random_state(self.random_state)
        cv = self.cv(random_state=rng, shuffle=True,
                     **self.cvargs)
        return cv.get_n_splits(X, y, groups) * self.n_repeats

    def __repr__(self):
        return _build_repr(self)


class RepeatedKFold(_RepeatedSplits):
    """Repeated K-Fold cross validator.

    Repeats K-Fold n times with different randomization in each repetition.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

    n_repeats : int, default=10
        Number of times cross-validator needs to be repeated.

    random_state : int, RandomState instance or None, optional, default=None
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import RepeatedKFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> rkf = RepeatedKFold(n_splits=2, n_repeats=2, random_state=2652124)
    >>> for train_index, test_index in rkf.split(X):
    ...     print("TRAIN:", train_index, "TEST:", test_index)
    ...     X_train, X_test = X[train_index], X[test_index]
    ...     y_train, y_test = y[train_index], y[test_index]
    ...
    TRAIN: [0 1] TEST: [2 3]
    TRAIN: [2 3] TEST: [0 1]
    TRAIN: [1 2] TEST: [0 3]
    TRAIN: [0 3] TEST: [1 2]

    Notes
    -----
    Randomized CV splitters may return different results for each call of
    split. You can make the results identical by setting ``random_state``
    to an integer.

    See also
    --------
    RepeatedStratifiedKFold: Repeats Stratified K-Fold n times.
    """
    def __init__(self, n_splits=5, n_repeats=10, random_state=None):
        super().__init__(
            KFold, n_repeats, random_state, n_splits=n_splits)


class RepeatedStratifiedKFold(_RepeatedSplits):
    """Repeated Stratified K-Fold cross validator.

    Repeats Stratified K-Fold n times with different randomization in each
    repetition.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default=5
        Number of folds. Must be at least 2.

    n_repeats : int, default=10
        Number of times cross-validator needs to be repeated.

    random_state : None, int or RandomState, default=None
        Random state to be used to generate random state for each
        repetition.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import RepeatedStratifiedKFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> rskf = RepeatedStratifiedKFold(n_splits=2, n_repeats=2,
    ...     random_state=36851234)
    >>> for train_index, test_index in rskf.split(X, y):
    ...     print("TRAIN:", train_index, "TEST:", test_index)
    ...     X_train, X_test = X[train_index], X[test_index]
    ...     y_train, y_test = y[train_index], y[test_index]
    ...
    TRAIN: [1 2] TEST: [0 3]
    TRAIN: [0 3] TEST: [1 2]
    TRAIN: [1 3] TEST: [0 2]
    TRAIN: [0 2] TEST: [1 3]

    Notes
    -----
    Randomized CV splitters may return different results for each call of
    split. You can make the results identical by setting ``random_state``
    to an integer.

    See also
    --------
    RepeatedKFold: Repeats K-Fold n times.
    """
    def __init__(self, n_splits=5, n_repeats=10, random_state=None):
        super().__init__(
            StratifiedKFold, n_repeats, random_state, n_splits=n_splits)


class BaseShuffleSplit(metaclass=ABCMeta):
    """Base class for ShuffleSplit and StratifiedShuffleSplit"""

    def __init__(self, n_splits=10, test_size=None, train_size=None,
                 random_state=None):
        self.n_splits = n_splits
        self.test_size = test_size
        self.train_size = train_size
        self.random_state = random_state
        self._default_test_size = 0.1

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.

        Notes
        -----
        Randomized CV splitters may return different results for each call of
        split. You can make the results identical by setting ``random_state``
        to an integer.
        """
        X, y, groups = indexable(X, y, groups)
        for train, test in self._iter_indices(X, y, groups):
            yield train, test

    @abstractmethod
    def _iter_indices(self, X, y=None, groups=None):
        """Generate (train, test) indices"""

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return self.n_splits

    def __repr__(self):
        return _build_repr(self)


class ShuffleSplit(BaseShuffleSplit):
    """Random permutation cross-validator

    Yields indices to split data into training and test sets.

    Note: contrary to other cross-validation strategies, random splits
    do not guarantee that all folds will be different, although this is
    still very likely for sizeable datasets.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default 10
        Number of re-shuffling & splitting iterations.

    test_size : float, int, None, default=None
        If float, should be between 0.0 and 1.0 and represent the proportion
        of the dataset to include in the test split. If int, represents the
        absolute number of test samples. If None, the value is set to the
        complement of the train size. If ``train_size`` is also None, it will
        be set to 0.1.

    train_size : float, int, or None, default=None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import ShuffleSplit
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8], [3, 4], [5, 6]])
    >>> y = np.array([1, 2, 1, 2, 1, 2])
    >>> rs = ShuffleSplit(n_splits=5, test_size=.25, random_state=0)
    >>> rs.get_n_splits(X)
    5
    >>> print(rs)
    ShuffleSplit(n_splits=5, random_state=0, test_size=0.25, train_size=None)
    >>> for train_index, test_index in rs.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    TRAIN: [1 3 0 4] TEST: [5 2]
    TRAIN: [4 0 2 5] TEST: [1 3]
    TRAIN: [1 2 4 0] TEST: [3 5]
    TRAIN: [3 4 1 0] TEST: [5 2]
    TRAIN: [3 5 1 0] TEST: [2 4]
    >>> rs = ShuffleSplit(n_splits=5, train_size=0.5, test_size=.25,
    ...                   random_state=0)
    >>> for train_index, test_index in rs.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    TRAIN: [1 3 0] TEST: [5 2]
    TRAIN: [4 0 2] TEST: [1 3]
    TRAIN: [1 2 4] TEST: [3 5]
    TRAIN: [3 4 1] TEST: [5 2]
    TRAIN: [3 5 1] TEST: [2 4]
    """
    def __init__(self, n_splits=10, test_size=None, train_size=None,
                 random_state=None):
        super().__init__(
            n_splits=n_splits,
            test_size=test_size,
            train_size=train_size,
            random_state=random_state)
        self._default_test_size = 0.1

    def _iter_indices(self, X, y=None, groups=None):
        n_samples = _num_samples(X)
        n_train, n_test = _validate_shuffle_split(
            n_samples, self.test_size, self.train_size,
            default_test_size=self._default_test_size)

        rng = check_random_state(self.random_state)
        for i in range(self.n_splits):
            # random partition
            permutation = rng.permutation(n_samples)
            ind_test = permutation[:n_test]
            ind_train = permutation[n_test:(n_test + n_train)]
            yield ind_train, ind_test


class GroupShuffleSplit(ShuffleSplit):
    '''Shuffle-Group(s)-Out cross-validation iterator

    Provides randomized train/test indices to split data according to a
    third-party provided group. This group information can be used to encode
    arbitrary domain specific stratifications of the samples as integers.

    For instance the groups could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    The difference between LeavePGroupsOut and GroupShuffleSplit is that
    the former generates splits using all subsets of size ``p`` unique groups,
    whereas GroupShuffleSplit generates a user-determined number of random
    test splits, each with a user-determined fraction of unique groups.

    For example, a less computationally intensive alternative to
    ``LeavePGroupsOut(p=10)`` would be
    ``GroupShuffleSplit(test_size=10, n_splits=100)``.

    Note: The parameters ``test_size`` and ``train_size`` refer to groups, and
    not to samples, as in ShuffleSplit.


    Parameters
    ----------
    n_splits : int (default 5)
        Number of re-shuffling & splitting iterations.

    test_size : float, int, None, optional (default=None)
        If float, should be between 0.0 and 1.0 and represent the proportion
        of groups to include in the test split (rounded up). If int,
        represents the absolute number of test groups. If None, the value is
        set to the complement of the train size. By default, the value is set
        to 0.2.
        The default will change in version 0.21. It will remain 0.2 only
        if ``train_size`` is unspecified, otherwise it will complement
        the specified ``train_size``.

    train_size : float, int, or None, default is None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the groups to include in the train split. If
        int, represents the absolute number of train groups. If None,
        the value is automatically set to the complement of the test size.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import GroupShuffleSplit
    >>> X = np.ones(shape=(8, 2))
    >>> y = np.ones(shape=(8, 1))
    >>> groups = np.array([1, 1, 2, 2, 2, 3, 3, 3])
    >>> print(groups.shape)
    (8,)
    >>> gss = GroupShuffleSplit(n_splits=2, train_size=.7, random_state=42)
    >>> gss.get_n_splits()
    2
    >>> for train_idx, test_idx in gss.split(X, y, groups):
    ...    print("TRAIN:", train_idx, "TEST:", test_idx)
    TRAIN: [2 3 4 5 6 7] TEST: [0 1]
    TRAIN: [0 1 5 6 7] TEST: [2 3 4]
    '''

    def __init__(self, n_splits=5, test_size=None, train_size=None,
                 random_state=None):
        super().__init__(
            n_splits=n_splits,
            test_size=test_size,
            train_size=train_size,
            random_state=random_state)
        self._default_test_size = 0.2

    def _iter_indices(self, X, y, groups):
        if groups is None:
            raise ValueError("The 'groups' parameter should not be None.")
        groups = check_array(groups, ensure_2d=False, dtype=None)
        classes, group_indices = np.unique(groups, return_inverse=True)
        for group_train, group_test in super()._iter_indices(X=classes):
            # these are the indices of classes in the partition
            # invert them into data indices

            train = np.flatnonzero(np.in1d(group_indices, group_train))
            test = np.flatnonzero(np.in1d(group_indices, group_test))

            yield train, test

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,), optional
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.

        Notes
        -----
        Randomized CV splitters may return different results for each call of
        split. You can make the results identical by setting ``random_state``
        to an integer.
        """
        return super().split(X, y, groups)


class StratifiedShuffleSplit(BaseShuffleSplit):
    """Stratified ShuffleSplit cross-validator

    Provides train/test indices to split data in train/test sets.

    This cross-validation object is a merge of StratifiedKFold and
    ShuffleSplit, which returns stratified randomized folds. The folds
    are made by preserving the percentage of samples for each class.

    Note: like the ShuffleSplit strategy, stratified random splits
    do not guarantee that all folds will be different, although this is
    still very likely for sizeable datasets.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_splits : int, default 10
        Number of re-shuffling & splitting iterations.

    test_size : float, int, None, optional (default=None)
        If float, should be between 0.0 and 1.0 and represent the proportion
        of the dataset to include in the test split. If int, represents the
        absolute number of test samples. If None, the value is set to the
        complement of the train size. If ``train_size`` is also None, it will
        be set to 0.1.

    train_size : float, int, or None, default is None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import StratifiedShuffleSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 0, 1, 1, 1])
    >>> sss = StratifiedShuffleSplit(n_splits=5, test_size=0.5, random_state=0)
    >>> sss.get_n_splits(X, y)
    5
    >>> print(sss)
    StratifiedShuffleSplit(n_splits=5, random_state=0, ...)
    >>> for train_index, test_index in sss.split(X, y):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [5 2 3] TEST: [4 1 0]
    TRAIN: [5 1 4] TEST: [0 2 3]
    TRAIN: [5 0 2] TEST: [4 3 1]
    TRAIN: [4 1 0] TEST: [2 3 5]
    TRAIN: [0 5 1] TEST: [3 4 2]
    """

    def __init__(self, n_splits=10, test_size=None, train_size=None,
                 random_state=None):
        super().__init__(
            n_splits=n_splits,
            test_size=test_size,
            train_size=train_size,
            random_state=random_state)
        self._default_test_size = 0.1

    def _iter_indices(self, X, y, groups=None):
        n_samples = _num_samples(X)
        y = check_array(y, ensure_2d=False, dtype=None)
        n_train, n_test = _validate_shuffle_split(
            n_samples, self.test_size, self.train_size,
            default_test_size=self._default_test_size)

        if y.ndim == 2:
            # for multi-label y, map each distinct row to a string repr
            # using join because str(row) uses an ellipsis if len(row) > 1000
            y = np.array([' '.join(row.astype('str')) for row in y])

        classes, y_indices = np.unique(y, return_inverse=True)
        n_classes = classes.shape[0]

        class_counts = np.bincount(y_indices)
        if np.min(class_counts) < 2:
            raise ValueError("The least populated class in y has only 1"
                             " member, which is too few. The minimum"
                             " number of groups for any class cannot"
                             " be less than 2.")

        if n_train < n_classes:
            raise ValueError('The train_size = %d should be greater or '
                             'equal to the number of classes = %d' %
                             (n_train, n_classes))
        if n_test < n_classes:
            raise ValueError('The test_size = %d should be greater or '
                             'equal to the number of classes = %d' %
                             (n_test, n_classes))

        # Find the sorted list of instances for each class:
        # (np.unique above performs a sort, so code is O(n logn) already)
        class_indices = np.split(np.argsort(y_indices, kind='mergesort'),
                                 np.cumsum(class_counts)[:-1])

        rng = check_random_state(self.random_state)

        for _ in range(self.n_splits):
            # if there are ties in the class-counts, we want
            # to make sure to break them anew in each iteration
            n_i = _approximate_mode(class_counts, n_train, rng)
            class_counts_remaining = class_counts - n_i
            t_i = _approximate_mode(class_counts_remaining, n_test, rng)

            train = []
            test = []

            for i in range(n_classes):
                permutation = rng.permutation(class_counts[i])
                perm_indices_class_i = class_indices[i].take(permutation,
                                                             mode='clip')

                train.extend(perm_indices_class_i[:n_i[i]])
                test.extend(perm_indices_class_i[n_i[i]:n_i[i] + t_i[i]])

            train = rng.permutation(train)
            test = rng.permutation(test)

            yield train, test

    def split(self, X, y, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

            Note that providing ``y`` is sufficient to generate the splits and
            hence ``np.zeros(n_samples)`` may be used as a placeholder for
            ``X`` instead of actual training data.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.
            Stratification is done based on the y labels.

        groups : object
            Always ignored, exists for compatibility.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.

        Notes
        -----
        Randomized CV splitters may return different results for each call of
        split. You can make the results identical by setting ``random_state``
        to an integer.
        """
        y = check_array(y, ensure_2d=False, dtype=None)
        return super().split(X, y, groups)


def _validate_shuffle_split(n_samples, test_size, train_size,
                            default_test_size=None):
    """
    Validation helper to check if the test/test sizes are meaningful wrt to the
    size of the data (n_samples)
    """
    if test_size is None and train_size is None:
        test_size = default_test_size

    test_size_type = np.asarray(test_size).dtype.kind
    train_size_type = np.asarray(train_size).dtype.kind

    if (test_size_type == 'i' and (test_size >= n_samples or test_size <= 0)
       or test_size_type == 'f' and (test_size <= 0 or test_size >= 1)):
        raise ValueError('test_size={0} should be either positive and smaller'
                         ' than the number of samples {1} or a float in the '
                         '(0, 1) range'.format(test_size, n_samples))

    if (train_size_type == 'i' and (train_size >= n_samples or train_size <= 0)
       or train_size_type == 'f' and (train_size <= 0 or train_size >= 1)):
        raise ValueError('train_size={0} should be either positive and smaller'
                         ' than the number of samples {1} or a float in the '
                         '(0, 1) range'.format(train_size, n_samples))

    if train_size is not None and train_size_type not in ('i', 'f'):
        raise ValueError("Invalid value for train_size: {}".format(train_size))
    if test_size is not None and test_size_type not in ('i', 'f'):
        raise ValueError("Invalid value for test_size: {}".format(test_size))

    if (train_size_type == 'f' and test_size_type == 'f' and
            train_size + test_size > 1):
        raise ValueError(
            'The sum of test_size and train_size = {}, should be in the (0, 1)'
            ' range. Reduce test_size and/or train_size.'
            .format(train_size + test_size))

    if test_size_type == 'f':
        n_test = ceil(test_size * n_samples)
    elif test_size_type == 'i':
        n_test = float(test_size)

    if train_size_type == 'f':
        n_train = floor(train_size * n_samples)
    elif train_size_type == 'i':
        n_train = float(train_size)

    if train_size is None:
        n_train = n_samples - n_test
    elif test_size is None:
        n_test = n_samples - n_train

    if n_train + n_test > n_samples:
        raise ValueError('The sum of train_size and test_size = %d, '
                         'should be smaller than the number of '
                         'samples %d. Reduce test_size and/or '
                         'train_size.' % (n_train + n_test, n_samples))

    n_train, n_test = int(n_train), int(n_test)

    if n_train == 0:
        raise ValueError(
            'With n_samples={}, test_size={} and train_size={}, the '
            'resulting train set will be empty. Adjust any of the '
            'aforementioned parameters.'.format(n_samples, test_size,
                                                train_size)
        )

    return n_train, n_test


class PredefinedSplit(BaseCrossValidator):
    """Predefined split cross-validator

    Provides train/test indices to split data into train/test sets using a
    predefined scheme specified by the user with the ``test_fold`` parameter.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    test_fold : array-like, shape (n_samples,)
        The entry ``test_fold[i]`` represents the index of the test set that
        sample ``i`` belongs to. It is possible to exclude sample ``i`` from
        any test set (i.e. include sample ``i`` in every training set) by
        setting ``test_fold[i]`` equal to -1.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import PredefinedSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> test_fold = [0, 1, -1, 1]
    >>> ps = PredefinedSplit(test_fold)
    >>> ps.get_n_splits()
    2
    >>> print(ps)
    PredefinedSplit(test_fold=array([ 0,  1, -1,  1]))
    >>> for train_index, test_index in ps.split():
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 2 3] TEST: [0]
    TRAIN: [0 2] TEST: [1 3]
    """

    def __init__(self, test_fold):
        self.test_fold = np.array(test_fold, dtype=np.int)
        self.test_fold = column_or_1d(self.test_fold)
        self.unique_folds = np.unique(self.test_fold)
        self.unique_folds = self.unique_folds[self.unique_folds != -1]

    def split(self, X=None, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        ind = np.arange(len(self.test_fold))
        for test_index in self._iter_test_masks():
            train_index = ind[np.logical_not(test_index)]
            test_index = ind[test_index]
            yield train_index, test_index

    def _iter_test_masks(self):
        """Generates boolean masks corresponding to test sets."""
        for f in self.unique_folds:
            test_index = np.where(self.test_fold == f)[0]
            test_mask = np.zeros(len(self.test_fold), dtype=np.bool)
            test_mask[test_index] = True
            yield test_mask

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return len(self.unique_folds)


class _CVIterableWrapper(BaseCrossValidator):
    """Wrapper class for old style cv objects and iterables."""
    def __init__(self, cv):
        self.cv = list(cv)

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return len(self.cv)

    def split(self, X=None, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        groups : object
            Always ignored, exists for compatibility.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        for train, test in self.cv:
            yield train, test


def check_cv(cv=5, y=None, classifier=False):
    """Input checker utility for building a cross-validator

    Parameters
    ----------
    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 5-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.

        For integer/None inputs, if classifier is True and ``y`` is either
        binary or multiclass, :class:`StratifiedKFold` is used. In all other
        cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

        .. versionchanged:: 0.22
            ``cv`` default value changed from 3-fold to 5-fold.

    y : array-like, optional
        The target variable for supervised learning problems.

    classifier : boolean, optional, default False
        Whether the task is a classification task, in which case
        stratified KFold will be used.

    Returns
    -------
    checked_cv : a cross-validator instance.
        The return value is a cross-validator which generates the train/test
        splits via the ``split`` method.
    """
    cv = 5 if cv is None else cv
    if isinstance(cv, numbers.Integral):
        if (classifier and (y is not None) and
                (type_of_target(y) in ('binary', 'multiclass'))):
            return StratifiedKFold(cv)
        else:
            return KFold(cv)

    if not hasattr(cv, 'split') or isinstance(cv, str):
        if not isinstance(cv, Iterable) or isinstance(cv, str):
            raise ValueError("Expected cv as an integer, cross-validation "
                             "object (from sklearn.model_selection) "
                             "or an iterable. Got %s." % cv)
        return _CVIterableWrapper(cv)

    return cv  # New style cv objects are passed without any modification


def train_test_split(*arrays, **options):
    """Split arrays or matrices into random train and test subsets

    Quick utility that wraps input validation and
    ``next(ShuffleSplit().split(X, y))`` and application to input data
    into a single call for splitting (and optionally subsampling) data in a
    oneliner.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    *arrays : sequence of indexables with same length / shape[0]
        Allowed inputs are lists, numpy arrays, scipy-sparse
        matrices or pandas dataframes.

    test_size : float, int or None, optional (default=None)
        If float, should be between 0.0 and 1.0 and represent the proportion
        of the dataset to include in the test split. If int, represents the
        absolute number of test samples. If None, the value is set to the
        complement of the train size. If ``train_size`` is also None, it will
        be set to 0.25.

    train_size : float, int, or None, (default=None)
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.

    random_state : int, RandomState instance or None, optional (default=None)
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    shuffle : boolean, optional (default=True)
        Whether or not to shuffle the data before splitting. If shuffle=False
        then stratify must be None.

    stratify : array-like or None (default=None)
        If not None, data is split in a stratified fashion, using this as
        the class labels.

    Returns
    -------
    splitting : list, length=2 * len(arrays)
        List containing train-test split of inputs.

        .. versionadded:: 0.16
            If the input is sparse, the output will be a
            ``scipy.sparse.csr_matrix``. Else, output type is the same as the
            input type.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.model_selection import train_test_split
    >>> X, y = np.arange(10).reshape((5, 2)), range(5)
    >>> X
    array([[0, 1],
           [2, 3],
           [4, 5],
           [6, 7],
           [8, 9]])
    >>> list(y)
    [0, 1, 2, 3, 4]

    >>> X_train, X_test, y_train, y_test = train_test_split(
    ...     X, y, test_size=0.33, random_state=42)
    ...
    >>> X_train
    array([[4, 5],
           [0, 1],
           [6, 7]])
    >>> y_train
    [2, 0, 3]
    >>> X_test
    array([[2, 3],
           [8, 9]])
    >>> y_test
    [1, 4]

    >>> train_test_split(y, shuffle=False)
    [[0, 1, 2], [3, 4]]

    """
    n_arrays = len(arrays)
    if n_arrays == 0:
        raise ValueError("At least one array required as input")
    test_size = options.pop('test_size', None)
    train_size = options.pop('train_size', None)
    random_state = options.pop('random_state', None)
    stratify = options.pop('stratify', None)
    shuffle = options.pop('shuffle', True)

    if options:
        raise TypeError("Invalid parameters passed: %s" % str(options))

    arrays = indexable(*arrays)

    n_samples = _num_samples(arrays[0])
    n_train, n_test = _validate_shuffle_split(n_samples, test_size, train_size,
                                              default_test_size=0.25)

    if shuffle is False:
        if stratify is not None:
            raise ValueError(
                "Stratified train/test split is not implemented for "
                "shuffle=False")

        train = np.arange(n_train)
        test = np.arange(n_train, n_train + n_test)

    else:
        if stratify is not None:
            CVClass = StratifiedShuffleSplit
        else:
            CVClass = ShuffleSplit

        cv = CVClass(test_size=n_test,
                     train_size=n_train,
                     random_state=random_state)

        train, test = next(cv.split(X=arrays[0], y=stratify))

    return list(chain.from_iterable((safe_indexing(a, train),
                                     safe_indexing(a, test)) for a in arrays))


# Tell nose that train_test_split is not a test.
# (Needed for external libraries that may use nose.)
train_test_split.__test__ = False


def _build_repr(self):
    # XXX This is copied from BaseEstimator's get_params
    cls = self.__class__
    init = getattr(cls.__init__, 'deprecated_original', cls.__init__)
    # Ignore varargs, kw and default values and pop self
    init_signature = signature(init)
    # Consider the constructor parameters excluding 'self'
    if init is object.__init__:
        args = []
    else:
        args = sorted([p.name for p in init_signature.parameters.values()
                       if p.name != 'self' and p.kind != p.VAR_KEYWORD])
    class_name = self.__class__.__name__
    params = dict()
    for key in args:
        # We need deprecation warnings to always be on in order to
        # catch deprecated param values.
        # This is set in utils/__init__.py but it gets overwritten
        # when running under python3 somehow.
        warnings.simplefilter("always", DeprecationWarning)
        try:
            with warnings.catch_warnings(record=True) as w:
                value = getattr(self, key, None)
                if value is None and hasattr(self, 'cvargs'):
                    value = self.cvargs.get(key, None)
            if len(w) and w[0].category == DeprecationWarning:
                # if the parameter is deprecated, don't show it
                continue
        finally:
            warnings.filters.pop(0)
        params[key] = value

    return '%s(%s)' % (class_name, _pprint(params, offset=len(class_name)))
