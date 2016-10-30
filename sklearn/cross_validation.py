
"""
The :mod:`sklearn.cross_validation` module includes utilities for cross-
validation and performance evaluation.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>,
#         Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

from __future__ import print_function
from __future__ import division

import warnings
from itertools import chain, combinations
from math import ceil, floor, factorial
import numbers
import time
from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp

from .base import is_classifier, clone
from .utils import indexable, check_random_state, safe_indexing
from .utils.validation import (_is_arraylike, _num_samples,
                               column_or_1d)
from .utils.multiclass import type_of_target
from .utils.random import choice
from .externals.joblib import Parallel, delayed, logger
from .externals.six import with_metaclass
from .externals.six.moves import zip
from .metrics.scorer import check_scoring
from .utils.fixes import bincount
from .gaussian_process.kernels import Kernel as GPKernel
from .exceptions import FitFailedWarning


warnings.warn("This module was deprecated in version 0.18 in favor of the "
              "model_selection module into which all the refactored classes "
              "and functions are moved. Also note that the interface of the "
              "new CV iterators are different from that of this module. "
              "This module will be removed in 0.20.", DeprecationWarning)


__all__ = ['KFold',
           'LabelKFold',
           'LeaveOneLabelOut',
           'LeaveOneOut',
           'LeavePLabelOut',
           'LeavePOut',
           'ShuffleSplit',
           'StratifiedKFold',
           'StratifiedShuffleSplit',
           'PredefinedSplit',
           'LabelShuffleSplit',
           'check_cv',
           'cross_val_score',
           'cross_val_predict',
           'permutation_test_score',
           'train_test_split']


class _PartitionIterator(with_metaclass(ABCMeta)):
    """Base class for CV iterators where train_mask = ~test_mask

    Implementations must define `_iter_test_masks` or `_iter_test_indices`.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.
    """

    def __init__(self, n):
        if abs(n - int(n)) >= np.finfo('f').eps:
            raise ValueError("n must be an integer")
        self.n = int(n)

    def __iter__(self):
        ind = np.arange(self.n)
        for test_index in self._iter_test_masks():
            train_index = np.logical_not(test_index)
            train_index = ind[train_index]
            test_index = ind[test_index]
            yield train_index, test_index

    # Since subclasses must implement either _iter_test_masks or
    # _iter_test_indices, neither can be abstract.
    def _iter_test_masks(self):
        """Generates boolean masks corresponding to test sets.

        By default, delegates to _iter_test_indices()
        """
        for test_index in self._iter_test_indices():
            test_mask = self._empty_mask()
            test_mask[test_index] = True
            yield test_mask

    def _iter_test_indices(self):
        """Generates integer indices corresponding to test sets."""
        raise NotImplementedError

    def _empty_mask(self):
        return np.zeros(self.n, dtype=np.bool)


class LeaveOneOut(_PartitionIterator):
    """Leave-One-Out cross validation iterator.

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.LeaveOneOut` instead.

    Provides train/test indices to split data in train test sets. Each
    sample is used once as a test set (singleton) while the remaining
    samples form the training set.

    Note: ``LeaveOneOut(n)`` is equivalent to ``KFold(n, n_folds=n)`` and
    ``LeavePOut(n, p=1)``.

    Due to the high number of test sets (which is the same as the
    number of samples) this cross validation method can be very costly.
    For large datasets one should favor KFold, StratifiedKFold or
    ShuffleSplit.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4]])
    >>> y = np.array([1, 2])
    >>> loo = cross_validation.LeaveOneOut(2)
    >>> len(loo)
    2
    >>> print(loo)
    sklearn.cross_validation.LeaveOneOut(n=2)
    >>> for train_index, test_index in loo:
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
    LeaveOneLabelOut for splitting the data according to explicit,
    domain-specific stratification of the dataset.
    """

    def _iter_test_indices(self):
        return range(self.n)

    def __repr__(self):
        return '%s.%s(n=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
        )

    def __len__(self):
        return self.n


class LeavePOut(_PartitionIterator):
    """Leave-P-Out cross validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.LeavePOut` instead.

    Provides train/test indices to split data in train test sets. This results
    in testing on all distinct samples of size p, while the remaining n - p
    samples form the training set in each iteration.

    Note: ``LeavePOut(n, p)`` is NOT equivalent to ``KFold(n, n_folds=n // p)``
    which creates non-overlapping test sets.

    Due to the high number of iterations which grows combinatorically with the
    number of samples this cross validation method can be very costly. For
    large datasets one should favor KFold, StratifiedKFold or ShuffleSplit.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.

    p : int
        Size of the test sets.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 3, 4])
    >>> lpo = cross_validation.LeavePOut(4, 2)
    >>> len(lpo)
    6
    >>> print(lpo)
    sklearn.cross_validation.LeavePOut(n=4, p=2)
    >>> for train_index, test_index in lpo:
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

    def __init__(self, n, p):
        super(LeavePOut, self).__init__(n)
        self.p = p

    def _iter_test_indices(self):
        for comb in combinations(range(self.n), self.p):
            yield np.array(comb)

    def __repr__(self):
        return '%s.%s(n=%i, p=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
            self.p,
        )

    def __len__(self):
        return int(factorial(self.n) / factorial(self.n - self.p)
                   / factorial(self.p))


class _BaseKFold(with_metaclass(ABCMeta, _PartitionIterator)):
    """Base class to validate KFold approaches"""

    @abstractmethod
    def __init__(self, n, n_folds, shuffle, random_state):
        super(_BaseKFold, self).__init__(n)

        if abs(n_folds - int(n_folds)) >= np.finfo('f').eps:
            raise ValueError("n_folds must be an integer")
        self.n_folds = n_folds = int(n_folds)

        if n_folds <= 1:
            raise ValueError(
                "k-fold cross validation requires at least one"
                " train / test split by setting n_folds=2 or more,"
                " got n_folds={0}.".format(n_folds))
        if n_folds > self.n:
            raise ValueError(
                ("Cannot have number of folds n_folds={0} greater"
                 " than the number of samples: {1}.").format(n_folds, n))

        if not isinstance(shuffle, bool):
            raise TypeError("shuffle must be True or False;"
                            " got {0}".format(shuffle))
        self.shuffle = shuffle
        self.random_state = random_state


class KFold(_BaseKFold):
    """K-Folds cross validation iterator.

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.KFold` instead.

    Provides train/test indices to split data in train test sets. Split
    dataset into k consecutive folds (without shuffling by default).

    Each fold is then used as a validation set once while the k - 1 remaining
    fold(s) form the training set.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n : int
        Total number of elements.

    n_folds : int, default=3
        Number of folds. Must be at least 2.

    shuffle : boolean, optional
        Whether to shuffle the data before splitting into batches.

    random_state : None, int or RandomState
        When shuffle=True, pseudo-random number generator state used for
        shuffling. If None, use default numpy RNG for shuffling.

    Examples
    --------
    >>> from sklearn.cross_validation import KFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4])
    >>> kf = KFold(4, n_folds=2)
    >>> len(kf)
    2
    >>> print(kf)  # doctest: +NORMALIZE_WHITESPACE
    sklearn.cross_validation.KFold(n=4, n_folds=2, shuffle=False,
                                   random_state=None)
    >>> for train_index, test_index in kf:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [2 3] TEST: [0 1]
    TRAIN: [0 1] TEST: [2 3]

    Notes
    -----
    The first n % n_folds folds have size n // n_folds + 1, other folds have
    size n // n_folds.

    See also
    --------
    StratifiedKFold take label information into account to avoid building
    folds with imbalanced class distributions (for binary or multiclass
    classification tasks).

    LabelKFold: K-fold iterator variant with non-overlapping labels.
    """

    def __init__(self, n, n_folds=3, shuffle=False,
                 random_state=None):
        super(KFold, self).__init__(n, n_folds, shuffle, random_state)
        self.idxs = np.arange(n)
        if shuffle:
            rng = check_random_state(self.random_state)
            rng.shuffle(self.idxs)

    def _iter_test_indices(self):
        n = self.n
        n_folds = self.n_folds
        fold_sizes = (n // n_folds) * np.ones(n_folds, dtype=np.int)
        fold_sizes[:n % n_folds] += 1
        current = 0
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            yield self.idxs[start:stop]
            current = stop

    def __repr__(self):
        return '%s.%s(n=%i, n_folds=%i, shuffle=%s, random_state=%s)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
            self.n_folds,
            self.shuffle,
            self.random_state,
        )

    def __len__(self):
        return self.n_folds


class LabelKFold(_BaseKFold):
    """K-fold iterator variant with non-overlapping labels.

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.GroupKFold` instead.

    The same label will not appear in two different folds (the number of
    distinct labels has to be at least equal to the number of folds).

    The folds are approximately balanced in the sense that the number of
    distinct labels is approximately the same in each fold.

    .. versionadded:: 0.17

    Parameters
    ----------
    labels : array-like with shape (n_samples, )
        Contains a label for each sample.
        The folds are built so that the same label does not appear in two
        different folds.

    n_folds : int, default=3
        Number of folds. Must be at least 2.

    Examples
    --------
    >>> from sklearn.cross_validation import LabelKFold
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 3, 4])
    >>> labels = np.array([0, 0, 2, 2])
    >>> label_kfold = LabelKFold(labels, n_folds=2)
    >>> len(label_kfold)
    2
    >>> print(label_kfold)
    sklearn.cross_validation.LabelKFold(n_labels=4, n_folds=2)
    >>> for train_index, test_index in label_kfold:
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
    LeaveOneLabelOut for splitting the data according to explicit,
    domain-specific stratification of the dataset.
    """
    def __init__(self, labels, n_folds=3):
        super(LabelKFold, self).__init__(len(labels), n_folds,
                                         shuffle=False, random_state=None)

        unique_labels, labels = np.unique(labels, return_inverse=True)
        n_labels = len(unique_labels)

        if n_folds > n_labels:
            raise ValueError(
                ("Cannot have number of folds n_folds={0} greater"
                 " than the number of labels: {1}.").format(n_folds,
                                                            n_labels))

        # Weight labels by their number of occurrences
        n_samples_per_label = np.bincount(labels)

        # Distribute the most frequent labels first
        indices = np.argsort(n_samples_per_label)[::-1]
        n_samples_per_label = n_samples_per_label[indices]

        # Total weight of each fold
        n_samples_per_fold = np.zeros(n_folds)

        # Mapping from label index to fold index
        label_to_fold = np.zeros(len(unique_labels))

        # Distribute samples by adding the largest weight to the lightest fold
        for label_index, weight in enumerate(n_samples_per_label):
            lightest_fold = np.argmin(n_samples_per_fold)
            n_samples_per_fold[lightest_fold] += weight
            label_to_fold[indices[label_index]] = lightest_fold

        self.idxs = label_to_fold[labels]

    def _iter_test_indices(self):
        for f in range(self.n_folds):
            yield np.where(self.idxs == f)[0]

    def __repr__(self):
        return '{0}.{1}(n_labels={2}, n_folds={3})'.format(
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
            self.n_folds,
        )

    def __len__(self):
        return self.n_folds


class StratifiedKFold(_BaseKFold):
    """Stratified K-Folds cross validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.StratifiedKFold` instead.

    Provides train/test indices to split data in train test sets.

    This cross-validation object is a variation of KFold that
    returns stratified folds. The folds are made by preserving
    the percentage of samples for each class.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    y : array-like, [n_samples]
        Samples to split in K folds.

    n_folds : int, default=3
        Number of folds. Must be at least 2.

    shuffle : boolean, optional
        Whether to shuffle each stratification of the data before splitting
        into batches.

    random_state : None, int or RandomState
        When shuffle=True, pseudo-random number generator state used for
        shuffling. If None, use default numpy RNG for shuffling.

    Examples
    --------
    >>> from sklearn.cross_validation import StratifiedKFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> skf = StratifiedKFold(y, n_folds=2)
    >>> len(skf)
    2
    >>> print(skf)  # doctest: +NORMALIZE_WHITESPACE
    sklearn.cross_validation.StratifiedKFold(labels=[0 0 1 1], n_folds=2,
                                             shuffle=False, random_state=None)
    >>> for train_index, test_index in skf:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 3] TEST: [0 2]
    TRAIN: [0 2] TEST: [1 3]

    Notes
    -----
    All the folds have size trunc(n_samples / n_folds), the last one has the
    complementary.

    See also
    --------
    LabelKFold: K-fold iterator variant with non-overlapping labels.
    """

    def __init__(self, y, n_folds=3, shuffle=False,
                 random_state=None):
        super(StratifiedKFold, self).__init__(
            len(y), n_folds, shuffle, random_state)
        y = np.asarray(y)
        n_samples = y.shape[0]
        unique_labels, y_inversed = np.unique(y, return_inverse=True)
        label_counts = bincount(y_inversed)
        min_labels = np.min(label_counts)
        if np.all(self.n_folds > label_counts):
            raise ValueError("All the n_labels for individual classes"
                             " are less than %d folds."
                             % (self.n_folds))
        if self.n_folds > min_labels:
            warnings.warn(("The least populated class in y has only %d"
                           " members, which is too few. The minimum"
                           " number of labels for any class cannot"
                           " be less than n_folds=%d."
                           % (min_labels, self.n_folds)), Warning)

        # don't want to use the same seed in each label's shuffle
        if self.shuffle:
            rng = check_random_state(self.random_state)
        else:
            rng = self.random_state

        # pre-assign each sample to a test fold index using individual KFold
        # splitting strategies for each label so as to respect the
        # balance of labels
        per_label_cvs = [
            KFold(max(c, self.n_folds), self.n_folds, shuffle=self.shuffle,
                  random_state=rng) for c in label_counts]
        test_folds = np.zeros(n_samples, dtype=np.int)
        for test_fold_idx, per_label_splits in enumerate(zip(*per_label_cvs)):
            for label, (_, test_split) in zip(unique_labels, per_label_splits):
                label_test_folds = test_folds[y == label]
                # the test split can be too big because we used
                # KFold(max(c, self.n_folds), self.n_folds) instead of
                # KFold(c, self.n_folds) to make it possible to not crash even
                # if the data is not 100% stratifiable for all the labels
                # (we use a warning instead of raising an exception)
                # If this is the case, let's trim it:
                test_split = test_split[test_split < len(label_test_folds)]
                label_test_folds[test_split] = test_fold_idx
                test_folds[y == label] = label_test_folds

        self.test_folds = test_folds
        self.y = y

    def _iter_test_masks(self):
        for i in range(self.n_folds):
            yield self.test_folds == i

    def __repr__(self):
        return '%s.%s(labels=%s, n_folds=%i, shuffle=%s, random_state=%s)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.y,
            self.n_folds,
            self.shuffle,
            self.random_state,
        )

    def __len__(self):
        return self.n_folds


class LeaveOneLabelOut(_PartitionIterator):
    """Leave-One-Label_Out cross-validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.LeaveOneGroupOut` instead.

    Provides train/test indices to split data according to a third-party
    provided label. This label information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the labels could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    labels : array-like of int with shape (n_samples,)
        Arbitrary domain-specific stratification of the data to be used
        to draw the splits.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 1, 2])
    >>> labels = np.array([1, 1, 2, 2])
    >>> lol = cross_validation.LeaveOneLabelOut(labels)
    >>> len(lol)
    2
    >>> print(lol)
    sklearn.cross_validation.LeaveOneLabelOut(labels=[1 1 2 2])
    >>> for train_index, test_index in lol:
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

    See also
    --------
    LabelKFold: K-fold iterator variant with non-overlapping labels.
    """

    def __init__(self, labels):
        super(LeaveOneLabelOut, self).__init__(len(labels))
        # We make a copy of labels to avoid side-effects during iteration
        self.labels = np.array(labels, copy=True)
        self.unique_labels = np.unique(labels)
        self.n_unique_labels = len(self.unique_labels)

    def _iter_test_masks(self):
        for i in self.unique_labels:
            yield self.labels == i

    def __repr__(self):
        return '%s.%s(labels=%s)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.labels,
        )

    def __len__(self):
        return self.n_unique_labels


class LeavePLabelOut(_PartitionIterator):
    """Leave-P-Label_Out cross-validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.LeavePGroupsOut` instead.

    Provides train/test indices to split data according to a third-party
    provided label. This label information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the labels could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    The difference between LeavePLabelOut and LeaveOneLabelOut is that
    the former builds the test sets with all the samples assigned to
    ``p`` different values of the labels while the latter uses samples
    all assigned the same labels.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    labels : array-like of int with shape (n_samples,)
        Arbitrary domain-specific stratification of the data to be used
        to draw the splits.

    p : int
        Number of samples to leave out in the test split.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [5, 6]])
    >>> y = np.array([1, 2, 1])
    >>> labels = np.array([1, 2, 3])
    >>> lpl = cross_validation.LeavePLabelOut(labels, p=2)
    >>> len(lpl)
    3
    >>> print(lpl)
    sklearn.cross_validation.LeavePLabelOut(labels=[1 2 3], p=2)
    >>> for train_index, test_index in lpl:
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
    LabelKFold: K-fold iterator variant with non-overlapping labels.
    """

    def __init__(self, labels, p):
        # We make a copy of labels to avoid side-effects during iteration
        super(LeavePLabelOut, self).__init__(len(labels))
        self.labels = np.array(labels, copy=True)
        self.unique_labels = np.unique(labels)
        self.n_unique_labels = len(self.unique_labels)
        self.p = p

    def _iter_test_masks(self):
        comb = combinations(range(self.n_unique_labels), self.p)
        for idx in comb:
            test_index = self._empty_mask()
            idx = np.array(idx)
            for l in self.unique_labels[idx]:
                test_index[self.labels == l] = True
            yield test_index

    def __repr__(self):
        return '%s.%s(labels=%s, p=%s)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.labels,
            self.p,
        )

    def __len__(self):
        return int(factorial(self.n_unique_labels) /
                   factorial(self.n_unique_labels - self.p) /
                   factorial(self.p))


class BaseShuffleSplit(with_metaclass(ABCMeta)):
    """Base class for ShuffleSplit and StratifiedShuffleSplit"""

    def __init__(self, n, n_iter=10, test_size=0.1, train_size=None,
                 random_state=None):
        self.n = n
        self.n_iter = n_iter
        self.test_size = test_size
        self.train_size = train_size
        self.random_state = random_state
        self.n_train, self.n_test = _validate_shuffle_split(n, test_size,
                                                            train_size)

    def __iter__(self):
        for train, test in self._iter_indices():
            yield train, test
        return

    @abstractmethod
    def _iter_indices(self):
        """Generate (train, test) indices"""


class ShuffleSplit(BaseShuffleSplit):
    """Random permutation cross-validation iterator.

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.ShuffleSplit` instead.

    Yields indices to split data into training and test sets.

    Note: contrary to other cross-validation strategies, random splits
    do not guarantee that all folds will be different, although this is
    still very likely for sizeable datasets.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n : int
        Total number of elements in the dataset.

    n_iter : int (default 10)
        Number of re-shuffling & splitting iterations.

    test_size : float (default 0.1), int, or None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the test split. If
        int, represents the absolute number of test samples. If None,
        the value is automatically set to the complement of the train size.

    train_size : float, int, or None (default is None)
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.

    random_state : int or RandomState
        Pseudo-random number generator state used for random sampling.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> rs = cross_validation.ShuffleSplit(4, n_iter=3,
    ...     test_size=.25, random_state=0)
    >>> len(rs)
    3
    >>> print(rs)
    ... # doctest: +ELLIPSIS
    ShuffleSplit(4, n_iter=3, test_size=0.25, ...)
    >>> for train_index, test_index in rs:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...
    TRAIN: [3 1 0] TEST: [2]
    TRAIN: [2 1 3] TEST: [0]
    TRAIN: [0 2 1] TEST: [3]

    >>> rs = cross_validation.ShuffleSplit(4, n_iter=3,
    ...     train_size=0.5, test_size=.25, random_state=0)
    >>> for train_index, test_index in rs:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...
    TRAIN: [3 1] TEST: [2]
    TRAIN: [2 1] TEST: [0]
    TRAIN: [0 2] TEST: [3]

    """

    def _iter_indices(self):
        rng = check_random_state(self.random_state)
        for i in range(self.n_iter):
            # random partition
            permutation = rng.permutation(self.n)
            ind_test = permutation[:self.n_test]
            ind_train = permutation[self.n_test:self.n_test + self.n_train]
            yield ind_train, ind_test

    def __repr__(self):
        return ('%s(%d, n_iter=%d, test_size=%s, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.n,
                    self.n_iter,
                    str(self.test_size),
                    self.random_state,
                ))

    def __len__(self):
        return self.n_iter


def _validate_shuffle_split(n, test_size, train_size):
    if test_size is None and train_size is None:
        raise ValueError(
            'test_size and train_size can not both be None')

    if test_size is not None:
        if np.asarray(test_size).dtype.kind == 'f':
            if test_size >= 1.:
                raise ValueError(
                    'test_size=%f should be smaller '
                    'than 1.0 or be an integer' % test_size)
        elif np.asarray(test_size).dtype.kind == 'i':
            if test_size >= n:
                raise ValueError(
                    'test_size=%d should be smaller '
                    'than the number of samples %d' % (test_size, n))
        else:
            raise ValueError("Invalid value for test_size: %r" % test_size)

    if train_size is not None:
        if np.asarray(train_size).dtype.kind == 'f':
            if train_size >= 1.:
                raise ValueError("train_size=%f should be smaller "
                                 "than 1.0 or be an integer" % train_size)
            elif np.asarray(test_size).dtype.kind == 'f' and \
                    train_size + test_size > 1.:
                raise ValueError('The sum of test_size and train_size = %f, '
                                 'should be smaller than 1.0. Reduce '
                                 'test_size and/or train_size.' %
                                 (train_size + test_size))
        elif np.asarray(train_size).dtype.kind == 'i':
            if train_size >= n:
                raise ValueError("train_size=%d should be smaller "
                                 "than the number of samples %d" %
                                 (train_size, n))
        else:
            raise ValueError("Invalid value for train_size: %r" % train_size)

    if np.asarray(test_size).dtype.kind == 'f':
        n_test = ceil(test_size * n)
    elif np.asarray(test_size).dtype.kind == 'i':
        n_test = float(test_size)

    if train_size is None:
        n_train = n - n_test
    else:
        if np.asarray(train_size).dtype.kind == 'f':
            n_train = floor(train_size * n)
        else:
            n_train = float(train_size)

    if test_size is None:
        n_test = n - n_train

    if n_train + n_test > n:
        raise ValueError('The sum of train_size and test_size = %d, '
                         'should be smaller than the number of '
                         'samples %d. Reduce test_size and/or '
                         'train_size.' % (n_train + n_test, n))

    return int(n_train), int(n_test)


def _approximate_mode(class_counts, n_draws, rng):
    """Computes approximate mode of multivariate hypergeometric.

    This is an approximation to the mode of the multivariate
    hypergeometric given by class_counts and n_draws.
    It shouldn't be off by more than one.

    It is the mostly likely outcome of drawing n_draws many
    samples from the population given by class_counts.

    Parameters
    ----------
    class_counts : ndarray of int
        Population per class.
    n_draws : int
        Number of draws (samples to draw) from the overall population.
    rng : random state
        Used to break ties.

    Returns
    -------
    sampled_classes : ndarray of int
        Number of samples drawn from each class.
        np.sum(sampled_classes) == n_draws
    """
    # this computes a bad approximation to the mode of the
    # multivariate hypergeometric given by class_counts and n_draws
    continuous = n_draws * class_counts / class_counts.sum()
    # floored means we don't overshoot n_samples, but probably undershoot
    floored = np.floor(continuous)
    # we add samples according to how much "left over" probability
    # they had, until we arrive at n_samples
    need_to_add = int(n_draws - floored.sum())
    if need_to_add > 0:
        remainder = continuous - floored
        values = np.sort(np.unique(remainder))[::-1]
        # add according to remainder, but break ties
        # randomly to avoid biases
        for value in values:
            inds, = np.where(remainder == value)
            # if we need_to_add less than what's in inds
            # we draw randomly from them.
            # if we need to add more, we add them all and
            # go to the next value
            add_now = min(len(inds), need_to_add)
            inds = choice(inds, size=add_now, replace=False, random_state=rng)
            floored[inds] += 1
            need_to_add -= add_now
            if need_to_add == 0:
                    break
    return floored.astype(np.int)


class StratifiedShuffleSplit(BaseShuffleSplit):
    """Stratified ShuffleSplit cross validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.StratifiedShuffleSplit` instead.

    Provides train/test indices to split data in train test sets.

    This cross-validation object is a merge of StratifiedKFold and
    ShuffleSplit, which returns stratified randomized folds. The folds
    are made by preserving the percentage of samples for each class.

    Note: like the ShuffleSplit strategy, stratified random splits
    do not guarantee that all folds will be different, although this is
    still very likely for sizeable datasets.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    y : array, [n_samples]
        Labels of samples.

    n_iter : int (default 10)
        Number of re-shuffling & splitting iterations.

    test_size : float (default 0.1), int, or None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the test split. If
        int, represents the absolute number of test samples. If None,
        the value is automatically set to the complement of the train size.

    train_size : float, int, or None (default is None)
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.

    random_state : int or RandomState
        Pseudo-random number generator state used for random sampling.

    Examples
    --------
    >>> from sklearn.cross_validation import StratifiedShuffleSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> sss = StratifiedShuffleSplit(y, 3, test_size=0.5, random_state=0)
    >>> len(sss)
    3
    >>> print(sss)       # doctest: +ELLIPSIS
    StratifiedShuffleSplit(labels=[0 0 1 1], n_iter=3, ...)
    >>> for train_index, test_index in sss:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 2] TEST: [3 0]
    TRAIN: [0 2] TEST: [1 3]
    TRAIN: [0 2] TEST: [3 1]
    """

    def __init__(self, y, n_iter=10, test_size=0.1, train_size=None,
                 random_state=None):

        super(StratifiedShuffleSplit, self).__init__(
            len(y), n_iter, test_size, train_size, random_state)

        self.y = np.array(y)
        self.classes, self.y_indices = np.unique(y, return_inverse=True)
        n_cls = self.classes.shape[0]

        if np.min(bincount(self.y_indices)) < 2:
            raise ValueError("The least populated class in y has only 1"
                             " member, which is too few. The minimum"
                             " number of labels for any class cannot"
                             " be less than 2.")

        if self.n_train < n_cls:
            raise ValueError('The train_size = %d should be greater or '
                             'equal to the number of classes = %d' %
                             (self.n_train, n_cls))
        if self.n_test < n_cls:
            raise ValueError('The test_size = %d should be greater or '
                             'equal to the number of classes = %d' %
                             (self.n_test, n_cls))

    def _iter_indices(self):
        rng = check_random_state(self.random_state)
        cls_count = bincount(self.y_indices)

        for n in range(self.n_iter):
            # if there are ties in the class-counts, we want
            # to make sure to break them anew in each iteration
            n_i = _approximate_mode(cls_count, self.n_train, rng)
            class_counts_remaining = cls_count - n_i
            t_i = _approximate_mode(class_counts_remaining, self.n_test, rng)

            train = []
            test = []

            for i, _ in enumerate(self.classes):
                permutation = rng.permutation(cls_count[i])
                perm_indices_class_i = np.where(
                    (i == self.y_indices))[0][permutation]

                train.extend(perm_indices_class_i[:n_i[i]])
                test.extend(perm_indices_class_i[n_i[i]:n_i[i] + t_i[i]])
            train = rng.permutation(train)
            test = rng.permutation(test)

            yield train, test

    def __repr__(self):
        return ('%s(labels=%s, n_iter=%d, test_size=%s, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.y,
                    self.n_iter,
                    str(self.test_size),
                    self.random_state,
                ))

    def __len__(self):
        return self.n_iter


class PredefinedSplit(_PartitionIterator):
    """Predefined split cross validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.PredefinedSplit` instead.

    Splits the data into training/test set folds according to a predefined
    scheme. Each sample can be assigned to at most one test set fold, as
    specified by the user through the ``test_fold`` parameter.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    test_fold : "array-like, shape (n_samples,)
        test_fold[i] gives the test set fold of sample i. A value of -1
        indicates that the corresponding sample is not part of any test set
        folds, but will instead always be put into the training fold.

    Examples
    --------
    >>> from sklearn.cross_validation import PredefinedSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> ps = PredefinedSplit(test_fold=[0, 1, -1, 1])
    >>> len(ps)
    2
    >>> print(ps)       # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
    sklearn.cross_validation.PredefinedSplit(test_fold=[ 0  1 -1  1])
    >>> for train_index, test_index in ps:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 2 3] TEST: [0]
    TRAIN: [0 2] TEST: [1 3]
    """

    def __init__(self, test_fold):
        super(PredefinedSplit, self).__init__(len(test_fold))
        self.test_fold = np.array(test_fold, dtype=np.int)
        self.test_fold = column_or_1d(self.test_fold)
        self.unique_folds = np.unique(self.test_fold)
        self.unique_folds = self.unique_folds[self.unique_folds != -1]

    def _iter_test_indices(self):
        for f in self.unique_folds:
            yield np.where(self.test_fold == f)[0]

    def __repr__(self):
        return '%s.%s(test_fold=%s)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.test_fold)

    def __len__(self):
        return len(self.unique_folds)


class LabelShuffleSplit(ShuffleSplit):
    """Shuffle-Labels-Out cross-validation iterator

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :class:`sklearn.model_selection.GroupShuffleSplit` instead.

    Provides randomized train/test indices to split data according to a
    third-party provided label. This label information can be used to encode
    arbitrary domain specific stratifications of the samples as integers.

    For instance the labels could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    The difference between LeavePLabelOut and LabelShuffleSplit is that
    the former generates splits using all subsets of size ``p`` unique labels,
    whereas LabelShuffleSplit generates a user-determined number of random
    test splits, each with a user-determined fraction of unique labels.

    For example, a less computationally intensive alternative to
    ``LeavePLabelOut(labels, p=10)`` would be
    ``LabelShuffleSplit(labels, test_size=10, n_iter=100)``.

    Note: The parameters ``test_size`` and ``train_size`` refer to labels, and
    not to samples, as in ShuffleSplit.

    .. versionadded:: 0.17

    Parameters
    ----------
    labels :  array, [n_samples]
        Labels of samples

    n_iter : int (default 5)
        Number of re-shuffling and splitting iterations.

    test_size : float (default 0.2), int, or None
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the labels to include in the test split. If
        int, represents the absolute number of test labels. If None,
        the value is automatically set to the complement of the train size.

    train_size : float, int, or None (default is None)
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the labels to include in the train split. If
        int, represents the absolute number of train labels. If None,
        the value is automatically set to the complement of the test size.

    random_state : int or RandomState
        Pseudo-random number generator state used for random sampling.

    """
    def __init__(self, labels, n_iter=5, test_size=0.2, train_size=None,
                 random_state=None):

        classes, label_indices = np.unique(labels, return_inverse=True)

        super(LabelShuffleSplit, self).__init__(
            len(classes),
            n_iter=n_iter,
            test_size=test_size,
            train_size=train_size,
            random_state=random_state)

        self.labels = labels
        self.classes = classes
        self.label_indices = label_indices

    def __repr__(self):
        return ('%s(labels=%s, n_iter=%d, test_size=%s, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.labels,
                    self.n_iter,
                    str(self.test_size),
                    self.random_state,
                ))

    def __len__(self):
        return self.n_iter

    def _iter_indices(self):
        for label_train, label_test in super(LabelShuffleSplit,
                                             self)._iter_indices():
            # these are the indices of classes in the partition
            # invert them into data indices

            train = np.flatnonzero(np.in1d(self.label_indices, label_train))
            test = np.flatnonzero(np.in1d(self.label_indices, label_test))

            yield train, test


##############################################################################
def _index_param_value(X, v, indices):
    """Private helper function for parameter value indexing."""
    if not _is_arraylike(v) or _num_samples(v) != _num_samples(X):
        # pass through: skip indexing
        return v
    if sp.issparse(v):
        v = v.tocsr()
    return safe_indexing(v, indices)


def cross_val_predict(estimator, X, y=None, cv=None, n_jobs=1,
                      verbose=0, fit_params=None, pre_dispatch='2*n_jobs'):
    """Generate cross-validated estimates for each input data point

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :func:`sklearn.model_selection.cross_val_predict` instead.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    estimator : estimator object implementing 'fit' and 'predict'
        The object to use to fit the data.

    X : array-like
        The data to fit. Can be, for example a list, or an array at least 2d.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - An object to be used as a cross-validation generator.
        - An iterable yielding train/test splits.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    verbose : integer, optional
        The verbosity level.

    fit_params : dict, optional
        Parameters to pass to the fit method of the estimator.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    Returns
    -------
    preds : ndarray
        This is the result of calling 'predict'

    Examples
    --------
    >>> from sklearn import datasets, linear_model
    >>> from sklearn.cross_validation import cross_val_predict
    >>> diabetes = datasets.load_diabetes()
    >>> X = diabetes.data[:150]
    >>> y = diabetes.target[:150]
    >>> lasso = linear_model.Lasso()
    >>> y_pred = cross_val_predict(lasso, X, y)
    """
    X, y = indexable(X, y)

    cv = check_cv(cv, X, y, classifier=is_classifier(estimator))
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    preds_blocks = parallel(delayed(_fit_and_predict)(clone(estimator), X, y,
                                                      train, test, verbose,
                                                      fit_params)
                            for train, test in cv)

    preds = [p for p, _ in preds_blocks]
    locs = np.concatenate([loc for _, loc in preds_blocks])
    if not _check_is_partition(locs, _num_samples(X)):
        raise ValueError('cross_val_predict only works for partitions')
    inv_locs = np.empty(len(locs), dtype=int)
    inv_locs[locs] = np.arange(len(locs))

    # Check for sparse predictions
    if sp.issparse(preds[0]):
        preds = sp.vstack(preds, format=preds[0].format)
    else:
        preds = np.concatenate(preds)
    return preds[inv_locs]


def _fit_and_predict(estimator, X, y, train, test, verbose, fit_params):
    """Fit estimator and predict values for a given dataset split.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    estimator : estimator object implementing 'fit' and 'predict'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    train : array-like, shape (n_train_samples,)
        Indices of training samples.

    test : array-like, shape (n_test_samples,)
        Indices of test samples.

    verbose : integer
        The verbosity level.

    fit_params : dict or None
        Parameters that will be passed to ``estimator.fit``.

    Returns
    -------
    preds : sequence
        Result of calling 'estimator.predict'

    test : array-like
        This is the value of the test parameter
    """
    # Adjust length of sample weights
    fit_params = fit_params if fit_params is not None else {}
    fit_params = dict([(k, _index_param_value(X, v, train))
                      for k, v in fit_params.items()])

    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, _ = _safe_split(estimator, X, y, test, train)

    if y_train is None:
        estimator.fit(X_train, **fit_params)
    else:
        estimator.fit(X_train, y_train, **fit_params)
    preds = estimator.predict(X_test)
    return preds, test


def _check_is_partition(locs, n):
    """Check whether locs is a reordering of the array np.arange(n)

    Parameters
    ----------
    locs : ndarray
        integer array to test
    n : int
        number of expected elements

    Returns
    -------
    is_partition : bool
        True iff sorted(locs) is range(n)
    """
    if len(locs) != n:
        return False
    hit = np.zeros(n, bool)
    hit[locs] = True
    if not np.all(hit):
        return False
    return True


def cross_val_score(estimator, X, y=None, scoring=None, cv=None, n_jobs=1,
                    verbose=0, fit_params=None, pre_dispatch='2*n_jobs'):
    """Evaluate a score by cross-validation

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :func:`sklearn.model_selection.cross_val_score` instead.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like
        The data to fit. Can be, for example a list, or an array at least 2d.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - An object to be used as a cross-validation generator.
        - An iterable yielding train/test splits.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    verbose : integer, optional
        The verbosity level.

    fit_params : dict, optional
        Parameters to pass to the fit method of the estimator.

    pre_dispatch : int, or string, optional
        Controls the number of jobs that get dispatched during parallel
        execution. Reducing this number can be useful to avoid an
        explosion of memory consumption when more jobs get dispatched
        than CPUs can process. This parameter can be:

            - None, in which case all the jobs are immediately
              created and spawned. Use this for lightweight and
              fast-running jobs, to avoid delays due to on-demand
              spawning of the jobs

            - An int, giving the exact number of total jobs that are
              spawned

            - A string, giving an expression as a function of n_jobs,
              as in '2*n_jobs'

    Returns
    -------
    scores : array of float, shape=(len(list(cv)),)
        Array of scores of the estimator for each run of the cross validation.

    Examples
    --------
    >>> from sklearn import datasets, linear_model
    >>> from sklearn.cross_validation import cross_val_score
    >>> diabetes = datasets.load_diabetes()
    >>> X = diabetes.data[:150]
    >>> y = diabetes.target[:150]
    >>> lasso = linear_model.Lasso()
    >>> print(cross_val_score(lasso, X, y))  # doctest:  +ELLIPSIS
    [ 0.33150734  0.08022311  0.03531764]

    See Also
    ---------
    :func:`sklearn.metrics.make_scorer`:
        Make a scorer from a performance metric or loss function.

    """
    X, y = indexable(X, y)

    cv = check_cv(cv, X, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, scoring=scoring)
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    scores = parallel(delayed(_fit_and_score)(clone(estimator), X, y, scorer,
                                              train, test, verbose, None,
                                              fit_params)
                      for train, test in cv)
    return np.array(scores)[:, 0]


def _fit_and_score(estimator, X, y, scorer, train, test, verbose,
                   parameters, fit_params, return_train_score=False,
                   return_parameters=False, error_score='raise'):
    """Fit estimator and compute scores for a given dataset split.

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    scorer : callable
        A scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    train : array-like, shape (n_train_samples,)
        Indices of training samples.

    test : array-like, shape (n_test_samples,)
        Indices of test samples.

    verbose : integer
        The verbosity level.

    error_score : 'raise' (default) or numeric
        Value to assign to the score if an error occurs in estimator fitting.
        If set to 'raise', the error is raised. If a numeric value is given,
        FitFailedWarning is raised. This parameter does not affect the refit
        step, which will always raise the error.

    parameters : dict or None
        Parameters to be set on the estimator.

    fit_params : dict or None
        Parameters that will be passed to ``estimator.fit``.

    return_train_score : boolean, optional, default: False
        Compute and return score on training set.

    return_parameters : boolean, optional, default: False
        Return parameters that has been used for the estimator.

    Returns
    -------
    train_score : float, optional
        Score on training set, returned only if `return_train_score` is `True`.

    test_score : float
        Score on test set.

    n_test_samples : int
        Number of test samples.

    scoring_time : float
        Time spent for fitting and scoring in seconds.

    parameters : dict or None, optional
        The parameters that have been evaluated.
    """
    if verbose > 1:
        if parameters is None:
            msg = ''
        else:
            msg = '%s' % (', '.join('%s=%s' % (k, v)
                          for k, v in parameters.items()))
        print("[CV] %s %s" % (msg, (64 - len(msg)) * '.'))

    # Adjust length of sample weights
    fit_params = fit_params if fit_params is not None else {}
    fit_params = dict([(k, _index_param_value(X, v, train))
                      for k, v in fit_params.items()])

    if parameters is not None:
        estimator.set_params(**parameters)

    start_time = time.time()

    X_train, y_train = _safe_split(estimator, X, y, train)
    X_test, y_test = _safe_split(estimator, X, y, test, train)

    try:
        if y_train is None:
            estimator.fit(X_train, **fit_params)
        else:
            estimator.fit(X_train, y_train, **fit_params)

    except Exception as e:
        if error_score == 'raise':
            raise
        elif isinstance(error_score, numbers.Number):
            test_score = error_score
            if return_train_score:
                train_score = error_score
            warnings.warn("Classifier fit failed. The score on this train-test"
                          " partition for these parameters will be set to %f. "
                          "Details: \n%r" % (error_score, e), FitFailedWarning)
        else:
            raise ValueError("error_score must be the string 'raise' or a"
                             " numeric value. (Hint: if using 'raise', please"
                             " make sure that it has been spelled correctly.)"
                             )

    else:
        test_score = _score(estimator, X_test, y_test, scorer)
        if return_train_score:
            train_score = _score(estimator, X_train, y_train, scorer)

    scoring_time = time.time() - start_time

    if verbose > 2:
        msg += ", score=%f" % test_score
    if verbose > 1:
        end_msg = "%s -%s" % (msg, logger.short_format_time(scoring_time))
        print("[CV] %s %s" % ((64 - len(end_msg)) * '.', end_msg))

    ret = [train_score] if return_train_score else []
    ret.extend([test_score, _num_samples(X_test), scoring_time])
    if return_parameters:
        ret.append(parameters)
    return ret


def _safe_split(estimator, X, y, indices, train_indices=None):
    """Create subset of dataset and properly handle kernels."""
    if hasattr(estimator, 'kernel') and callable(estimator.kernel) \
       and not isinstance(estimator.kernel, GPKernel):
        # cannot compute the kernel values with custom function
        raise ValueError("Cannot use a custom kernel function. "
                         "Precompute the kernel matrix instead.")

    if not hasattr(X, "shape"):
        if getattr(estimator, "_pairwise", False):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        X_subset = [X[idx] for idx in indices]
    else:
        if getattr(estimator, "_pairwise", False):
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


def _score(estimator, X_test, y_test, scorer):
    """Compute the score of an estimator on a given test set."""
    if y_test is None:
        score = scorer(estimator, X_test)
    else:
        score = scorer(estimator, X_test, y_test)
    if hasattr(score, 'item'):
        try:
            # e.g. unwrap memmapped scalars
            score = score.item()
        except ValueError:
            # non-scalar?
            pass
    if not isinstance(score, numbers.Number):
        raise ValueError("scoring must return a number, got %s (%s) instead."
                         % (str(score), type(score)))
    return score


def _permutation_test_score(estimator, X, y, cv, scorer):
    """Auxiliary function for permutation_test_score"""
    avg_score = []
    for train, test in cv:
        estimator.fit(X[train], y[train])
        avg_score.append(scorer(estimator, X[test], y[test]))
    return np.mean(avg_score)


def _shuffle(y, labels, random_state):
    """Return a shuffled copy of y eventually shuffle among same labels."""
    if labels is None:
        ind = random_state.permutation(len(y))
    else:
        ind = np.arange(len(labels))
        for label in np.unique(labels):
            this_mask = (labels == label)
            ind[this_mask] = random_state.permutation(ind[this_mask])
    return y[ind]


def check_cv(cv, X=None, y=None, classifier=False):
    """Input checker utility for building a CV in a user friendly way.

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :func:`sklearn.model_selection.check_cv` instead.

    Parameters
    ----------
    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - An object to be used as a cross-validation generator.
        - An iterable yielding train/test splits.

        For integer/None inputs, if classifier is True and ``y`` is binary or
        multiclass, :class:`StratifiedKFold` is used. In all other cases,
        :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    X : array-like
        The data the cross-val object will be applied on.

    y : array-like
        The target variable for a supervised learning problem.

    classifier : boolean optional
        Whether the task is a classification task, in which case
        stratified KFold will be used.

    Returns
    -------
    checked_cv: a cross-validation generator instance.
        The return value is guaranteed to be a cv generator instance, whatever
        the input type.
    """
    is_sparse = sp.issparse(X)
    if cv is None:
        cv = 3
    if isinstance(cv, numbers.Integral):
        if classifier:
            if type_of_target(y) in ['binary', 'multiclass']:
                cv = StratifiedKFold(y, cv)
            else:
                cv = KFold(_num_samples(y), cv)
        else:
            if not is_sparse:
                n_samples = len(X)
            else:
                n_samples = X.shape[0]
            cv = KFold(n_samples, cv)
    return cv


def permutation_test_score(estimator, X, y, cv=None,
                           n_permutations=100, n_jobs=1, labels=None,
                           random_state=0, verbose=0, scoring=None):
    """Evaluate the significance of a cross-validated score with permutations

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :func:`sklearn.model_selection.permutation_test_score` instead.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like
        The target variable to try to predict in the case of
        supervised learning.

    scoring : string, callable or None, optional, default: None
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.

    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:

        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - An object to be used as a cross-validation generator.
        - An iterable yielding train/test splits.

        For integer/None inputs, if the estimator is a classifier and ``y`` is
        either binary or multiclass, :class:`StratifiedKFold` is used. In all
        other cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

    n_permutations : integer, optional
        Number of times to permute ``y``.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    labels : array-like of shape [n_samples] (optional)
        Labels constrain the permutation among groups of samples with
        a same label.

    random_state : RandomState or an int seed (0 by default)
        A random number generator instance to define the state of the
        random permutations generator.

    verbose : integer, optional
        The verbosity level.

    Returns
    -------
    score : float
        The true score without permuting targets.

    permutation_scores : array, shape (n_permutations,)
        The scores obtained for each permutations.

    pvalue : float
        The returned value equals p-value if `scoring` returns bigger
        numbers for better scores (e.g., accuracy_score). If `scoring` is
        rather a loss function (i.e. when lower is better such as with
        `mean_squared_error`) then this is actually the complement of the
        p-value:  1 - p-value.

    Notes
    -----
    This function implements Test 1 in:

        Ojala and Garriga. Permutation Tests for Studying Classifier
        Performance.  The Journal of Machine Learning Research (2010)
        vol. 11

    """
    X, y = indexable(X, y)
    cv = check_cv(cv, X, y, classifier=is_classifier(estimator))
    scorer = check_scoring(estimator, scoring=scoring)
    random_state = check_random_state(random_state)

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    score = _permutation_test_score(clone(estimator), X, y, cv, scorer)
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_permutation_test_score)(
            clone(estimator), X, _shuffle(y, labels, random_state), cv,
            scorer)
        for _ in range(n_permutations))
    permutation_scores = np.array(permutation_scores)
    pvalue = (np.sum(permutation_scores >= score) + 1.0) / (n_permutations + 1)
    return score, permutation_scores, pvalue


permutation_test_score.__test__ = False  # to avoid a pb with nosetests


def train_test_split(*arrays, **options):
    """Split arrays or matrices into random train and test subsets

    .. deprecated:: 0.18
        This module will be removed in 0.20.
        Use :func:`sklearn.model_selection.train_test_split` instead.

    Quick utility that wraps input validation and
    ``next(iter(ShuffleSplit(n_samples)))`` and application to input
    data into a single call for splitting (and optionally subsampling)
    data in a oneliner.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    *arrays : sequence of indexables with same length / shape[0]
        Allowed inputs are lists, numpy arrays, scipy-sparse
        matrices or pandas dataframes.

    test_size : float, int, or None (default is None)
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the test split. If
        int, represents the absolute number of test samples. If None,
        the value is automatically set to the complement of the train size.
        If train size is also None, test size is set to 0.25.

    train_size : float, int, or None (default is None)
        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split. If
        int, represents the absolute number of train samples. If None,
        the value is automatically set to the complement of the test size.

    random_state : int or RandomState
        Pseudo-random number generator state used for random sampling.

    stratify : array-like or None (default is None)
        If not None, data is split in a stratified fashion, using this as
        the labels array.

        .. versionadded:: 0.17
           *stratify* splitting

    Returns
    -------
    splitting : list, length = 2 * len(arrays),
        List containing train-test split of inputs.

        .. versionadded:: 0.16
            If the input is sparse, the output will be a
            ``scipy.sparse.csr_matrix``. Else, output type is the same as the
            input type.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.cross_validation import train_test_split
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

    """
    n_arrays = len(arrays)
    if n_arrays == 0:
        raise ValueError("At least one array required as input")

    test_size = options.pop('test_size', None)
    train_size = options.pop('train_size', None)
    random_state = options.pop('random_state', None)
    stratify = options.pop('stratify', None)

    if options:
        raise TypeError("Invalid parameters passed: %s" % str(options))

    if test_size is None and train_size is None:
        test_size = 0.25
    arrays = indexable(*arrays)
    if stratify is not None:
        cv = StratifiedShuffleSplit(stratify, test_size=test_size,
                                    train_size=train_size,
                                    random_state=random_state)
    else:
        n_samples = _num_samples(arrays[0])
        cv = ShuffleSplit(n_samples, test_size=test_size,
                          train_size=train_size,
                          random_state=random_state)

    train, test = next(iter(cv))
    return list(chain.from_iterable((safe_indexing(a, train),
                                     safe_indexing(a, test)) for a in arrays))


train_test_split.__test__ = False  # to avoid a pb with nosetests
