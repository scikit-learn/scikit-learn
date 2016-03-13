"""
The :mod:`sklearn.model_selection._split` module includes classes and
functions to split the data based on a preset strategy.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>,
#         Olivier Girsel <olivier.grisel@ensta.org>
#         Raghav R V <rvraghav93@gmail.com>
# License: BSD 3 clause


from __future__ import print_function
from __future__ import division

import warnings
import inspect
from itertools import chain, combinations
from collections import Iterable
from math import ceil, floor
import numbers
from abc import ABCMeta, abstractmethod

import numpy as np

from scipy.misc import comb
from ..utils import indexable, check_random_state, safe_indexing
from ..utils.validation import _num_samples, column_or_1d
from ..utils.multiclass import type_of_target
from ..externals.six import with_metaclass
from ..externals.six.moves import zip
from ..utils.fixes import bincount
from ..utils.fixes import signature
from ..base import _pprint
from ..gaussian_process.kernels import Kernel as GPKernel

__all__ = ['BaseCrossValidator',
           'KFold',
           'LabelKFold',
           'LeaveOneLabelOut',
           'LeaveOneOut',
           'LeavePLabelOut',
           'LeavePOut',
           'ShuffleSplit',
           'LabelShuffleSplit',
           'StratifiedKFold',
           'StratifiedShuffleSplit',
           'PredefinedSplit',
           'train_test_split',
           'check_cv']


class BaseCrossValidator(with_metaclass(ABCMeta)):
    """Base class for all cross-validators

    Implementations must define `_iter_test_masks` or `_iter_test_indices`.
    """

    def __init__(self):
        # We need this for the build_repr to work properly in py2.7
        # see #6304
        pass

    def split(self, X, y=None, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, of length n_samples
            The target variable for supervised learning problems.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        X, y, labels = indexable(X, y, labels)
        indices = np.arange(_num_samples(X))
        for test_index in self._iter_test_masks(X, y, labels):
            train_index = indices[np.logical_not(test_index)]
            test_index = indices[test_index]
            yield train_index, test_index

    # Since subclasses must implement either _iter_test_masks or
    # _iter_test_indices, neither can be abstract.
    def _iter_test_masks(self, X=None, y=None, labels=None):
        """Generates boolean masks corresponding to test sets.

        By default, delegates to _iter_test_indices(X, y, labels)
        """
        for test_index in self._iter_test_indices(X, y, labels):
            test_mask = np.zeros(_num_samples(X), dtype=np.bool)
            test_mask[test_index] = True
            yield test_mask

    def _iter_test_indices(self, X=None, y=None, labels=None):
        """Generates integer indices corresponding to test sets."""
        raise NotImplementedError

    @abstractmethod
    def get_n_splits(self, X=None, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator"""

    def __repr__(self):
        return _build_repr(self)


class LeaveOneOut(BaseCrossValidator):
    """Leave-One-Out cross-validator

    Provides train/test indices to split data in train/test sets. Each
    sample is used once as a test set (singleton) while the remaining
    samples form the training set.

    Note: ``LeaveOneOut()`` is equivalent to ``KFold(n_folds=n)`` and
    ``LeavePOut(p=1)`` where ``n`` is the number of samples.

    Due to the high number of test sets (which is the same as the
    number of samples) this cross-validation method can be very costly.
    For large datasets one should favor :class:`KFold`, :class:`ShuffleSplit`
    or :class:`StratifiedKFold`.

    Read more in the :ref:`User Guide <cross_validation>`.

    Examples
    --------
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
    LeaveOneLabelOut
        For splitting the data according to explicit, domain-specific
        stratification of the dataset.

    LabelKFold: K-fold iterator variant with non-overlapping labels.
    """

    def _iter_test_indices(self, X, y=None, labels=None):
        return range(_num_samples(X))

    def get_n_splits(self, X, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if X is None:
            raise ValueError("The X parameter should not be None")
        return _num_samples(X)


class LeavePOut(BaseCrossValidator):
    """Leave-P-Out cross-validator

    Provides train/test indices to split data in train/test sets. This results
    in testing on all distinct samples of size p, while the remaining n - p
    samples form the training set in each iteration.

    Note: ``LeavePOut(p)`` is NOT equivalent to
    ``KFold(n_folds=n_samples // p)`` which creates non-overlapping test sets.

    Due to the high number of iterations which grows combinatorically with the
    number of samples this cross-validation method can be very costly. For
    large datasets one should favor :class:`KFold`, :class:`StratifiedKFold`
    or :class:`ShuffleSplit`.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    p : int
        Size of the test sets.

    Examples
    --------
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

    def _iter_test_indices(self, X, y=None, labels=None):
        for combination in combinations(range(_num_samples(X)), self.p):
            yield np.array(combination)

    def get_n_splits(self, X, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.
        """
        if X is None:
            raise ValueError("The X parameter should not be None")
        return int(comb(_num_samples(X), self.p, exact=True))


class _BaseKFold(with_metaclass(ABCMeta, BaseCrossValidator)):
    """Base class for KFold and StratifiedKFold"""

    @abstractmethod
    def __init__(self, n_folds, shuffle, random_state):
        if not isinstance(n_folds, numbers.Integral):
            raise ValueError('The number of folds must be of Integral type. '
                             '%s of type %s was passed.'
                             % (n_folds, type(n_folds)))
        n_folds = int(n_folds)

        if n_folds <= 1:
            raise ValueError(
                "k-fold cross-validation requires at least one"
                " train/test split by setting n_folds=2 or more,"
                " got n_folds={0}.".format(n_folds))

        if not isinstance(shuffle, bool):
            raise TypeError("shuffle must be True or False;"
                            " got {0}".format(shuffle))

        self.n_folds = n_folds
        self.shuffle = shuffle
        self.random_state = random_state

    def split(self, X, y=None, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        X, y, labels = indexable(X, y, labels)
        n_samples = _num_samples(X)
        if self.n_folds > n_samples:
            raise ValueError(
                ("Cannot have number of folds n_folds={0} greater"
                 " than the number of samples: {1}.").format(self.n_folds,
                                                             n_samples))

        for train, test in super(_BaseKFold, self).split(X, y, labels):
            yield train, test

    def get_n_splits(self, X=None, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return self.n_folds


class KFold(_BaseKFold):
    """K-Folds cross-validator

    Provides train/test indices to split data in train/test sets. Split
    dataset into k consecutive folds (without shuffling by default).

    Each fold is then used once as a validation while the k - 1 remaining
    folds form the training set.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
    n_folds : int, default=3
        Number of folds. Must be at least 2.

    shuffle : boolean, optional
        Whether to shuffle the data before splitting into batches.

    random_state : None, int or RandomState
        When shuffle=True, pseudo-random number generator state used for
        shuffling. If None, use default numpy RNG for shuffling.

    Examples
    --------
    >>> from sklearn.model_selection import KFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4])
    >>> kf = KFold(n_folds=2)
    >>> kf.get_n_splits(X)
    2
    >>> print(kf)  # doctest: +NORMALIZE_WHITESPACE
    KFold(n_folds=2, random_state=None, shuffle=False)
    >>> for train_index, test_index in kf.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [2 3] TEST: [0 1]
    TRAIN: [0 1] TEST: [2 3]

    Notes
    -----
    The first ``n_samples % n_folds`` folds have size
    ``n_samples // n_folds + 1``, other folds have size
    ``n_samples // n_folds``, where ``n_samples`` is the number of samples.

    See also
    --------
    StratifiedKFold
        Takes label information into account to avoid building folds with
        imbalanced class distributions (for binary or multiclass
        classification tasks).

    LabelKFold: K-fold iterator variant with non-overlapping labels.
    """

    def __init__(self, n_folds=3, shuffle=False,
                 random_state=None):
        super(KFold, self).__init__(n_folds, shuffle, random_state)

    def _iter_test_indices(self, X, y=None, labels=None):
        n_samples = _num_samples(X)
        indices = np.arange(n_samples)
        if self.shuffle:
            check_random_state(self.random_state).shuffle(indices)

        n_folds = self.n_folds
        fold_sizes = (n_samples // n_folds) * np.ones(n_folds, dtype=np.int)
        fold_sizes[:n_samples % n_folds] += 1
        current = 0
        for fold_size in fold_sizes:
            start, stop = current, current + fold_size
            yield indices[start:stop]
            current = stop


class LabelKFold(_BaseKFold):
    """K-fold iterator variant with non-overlapping labels.

    The same label will not appear in two different folds (the number of
    distinct labels has to be at least equal to the number of folds).

    The folds are approximately balanced in the sense that the number of
    distinct labels is approximately the same in each fold.

    Parameters
    ----------
    n_folds : int, default=3
        Number of folds. Must be at least 2.

    Examples
    --------
    >>> from sklearn.model_selection import LabelKFold
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 3, 4])
    >>> labels = np.array([0, 0, 2, 2])
    >>> label_kfold = LabelKFold(n_folds=2)
    >>> label_kfold.get_n_splits(X, y, labels)
    2
    >>> print(label_kfold)
    LabelKFold(n_folds=2)
    >>> for train_index, test_index in label_kfold.split(X, y, labels):
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
    LeaveOneLabelOut
        For splitting the data according to explicit domain-specific
        stratification of the dataset.
    """
    def __init__(self, n_folds=3):
        super(LabelKFold, self).__init__(n_folds, shuffle=False,
                                         random_state=None)

    def _iter_test_indices(self, X, y, labels):
        if labels is None:
            raise ValueError("The labels parameter should not be None")

        unique_labels, labels = np.unique(labels, return_inverse=True)
        n_labels = len(unique_labels)

        if self.n_folds > n_labels:
            raise ValueError("Cannot have number of folds n_folds=%d greater"
                             " than the number of labels: %d."
                             % (self.n_folds, n_labels))

        # Weight labels by their number of occurrences
        n_samples_per_label = np.bincount(labels)

        # Distribute the most frequent labels first
        indices = np.argsort(n_samples_per_label)[::-1]
        n_samples_per_label = n_samples_per_label[indices]

        # Total weight of each fold
        n_samples_per_fold = np.zeros(self.n_folds)

        # Mapping from label index to fold index
        label_to_fold = np.zeros(len(unique_labels))

        # Distribute samples by adding the largest weight to the lightest fold
        for label_index, weight in enumerate(n_samples_per_label):
            lightest_fold = np.argmin(n_samples_per_fold)
            n_samples_per_fold[lightest_fold] += weight
            label_to_fold[indices[label_index]] = lightest_fold

        indices = label_to_fold[labels]

        for f in range(self.n_folds):
            yield np.where(indices == f)[0]


class StratifiedKFold(_BaseKFold):
    """Stratified K-Folds cross-validator

    Provides train/test indices to split data in train/test sets.

    This cross-validation object is a variation of KFold that returns
    stratified folds. The folds are made by preserving the percentage of
    samples for each class.

    Read more in the :ref:`User Guide <cross_validation>`.

    Parameters
    ----------
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
    >>> from sklearn.model_selection import StratifiedKFold
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> skf = StratifiedKFold(n_folds=2)
    >>> skf.get_n_splits(X, y)
    2
    >>> print(skf)  # doctest: +NORMALIZE_WHITESPACE
    StratifiedKFold(n_folds=2, random_state=None, shuffle=False)
    >>> for train_index, test_index in skf.split(X, y):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 3] TEST: [0 2]
    TRAIN: [0 2] TEST: [1 3]

    Notes
    -----
    All the folds have size ``trunc(n_samples / n_folds)``, the last one has
    the complementary.

    """

    def __init__(self, n_folds=3, shuffle=False, random_state=None):
        super(StratifiedKFold, self).__init__(n_folds, shuffle, random_state)

    def _make_test_folds(self, X, y=None, labels=None):
        if self.shuffle:
            rng = check_random_state(self.random_state)
        else:
            rng = self.random_state
        y = np.asarray(y)
        n_samples = y.shape[0]
        unique_y, y_inversed = np.unique(y, return_inverse=True)
        y_counts = bincount(y_inversed)
        min_labels = np.min(y_counts)
        if np.all(self.n_folds > y_counts):
            raise ValueError("All the n_labels for individual classes"
                             " are less than %d folds."
                             % (self.n_folds))
        if self.n_folds > min_labels:
            warnings.warn(("The least populated class in y has only %d"
                           " members, which is too few. The minimum"
                           " number of labels for any class cannot"
                           " be less than n_folds=%d."
                           % (min_labels, self.n_folds)), Warning)

        # pre-assign each sample to a test fold index using individual KFold
        # splitting strategies for each class so as to respect the balance of
        # classes
        # NOTE: Passing the data corresponding to ith class say X[y==class_i]
        # will break when the data is not 100% stratifiable for all classes.
        # So we pass np.zeroes(max(c, n_folds)) as data to the KFold
        per_cls_cvs = [
            KFold(self.n_folds, shuffle=self.shuffle,
                  random_state=rng).split(np.zeros(max(count, self.n_folds)))
            for count in y_counts]

        test_folds = np.zeros(n_samples, dtype=np.int)
        for test_fold_indices, per_cls_splits in enumerate(zip(*per_cls_cvs)):
            for cls, (_, test_split) in zip(unique_y, per_cls_splits):
                cls_test_folds = test_folds[y == cls]
                # the test split can be too big because we used
                # KFold(...).split(X[:max(c, n_folds)]) when data is not 100%
                # stratifiable for all the classes
                # (we use a warning instead of raising an exception)
                # If this is the case, let's trim it:
                test_split = test_split[test_split < len(cls_test_folds)]
                cls_test_folds[test_split] = test_fold_indices
                test_folds[y == cls] = cls_test_folds

        return test_folds

    def _iter_test_masks(self, X, y=None, labels=None):
        test_folds = self._make_test_folds(X, y)
        for i in range(self.n_folds):
            yield test_folds == i

    def split(self, X, y, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        return super(StratifiedKFold, self).split(X, y, labels)

class LeaveOneLabelOut(BaseCrossValidator):
    """Leave One Label Out cross-validator

    Provides train/test indices to split data according to a third-party
    provided label. This label information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the labels could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    Read more in the :ref:`User Guide <cross_validation>`.

    Examples
    --------
    >>> from sklearn.model_selection import LeaveOneLabelOut
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 1, 2])
    >>> labels = np.array([1, 1, 2, 2])
    >>> lol = LeaveOneLabelOut()
    >>> lol.get_n_splits(X, y, labels)
    2
    >>> print(lol)
    LeaveOneLabelOut()
    >>> for train_index, test_index in lol.split(X, y, labels):
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

    def _iter_test_masks(self, X, y, labels):
        if labels is None:
            raise ValueError("The labels parameter should not be None")
        # We make a copy of labels to avoid side-effects during iteration
        labels = np.array(labels, copy=True)
        unique_labels = np.unique(labels)
        for i in unique_labels:
            yield labels == i

    def get_n_splits(self, X, y, labels):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if labels is None:
            raise ValueError("The labels parameter should not be None")
        return len(np.unique(labels))


class LeavePLabelOut(BaseCrossValidator):
    """Leave P Labels Out cross-validator

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
    n_labels : int
        Number of labels (``p``) to leave out in the test split.

    Examples
    --------
    >>> from sklearn.model_selection import LeavePLabelOut
    >>> X = np.array([[1, 2], [3, 4], [5, 6]])
    >>> y = np.array([1, 2, 1])
    >>> labels = np.array([1, 2, 3])
    >>> lpl = LeavePLabelOut(n_labels=2)
    >>> lpl.get_n_splits(X, y, labels)
    3
    >>> print(lpl)
    LeavePLabelOut(n_labels=2)
    >>> for train_index, test_index in lpl.split(X, y, labels):
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

    def __init__(self, n_labels):
        self.n_labels = n_labels

    def _iter_test_masks(self, X, y, labels):
        if labels is None:
            raise ValueError("The labels parameter should not be None")
        labels = np.array(labels, copy=True)
        unique_labels = np.unique(labels)
        combi = combinations(range(len(unique_labels)), self.n_labels)
        for indices in combi:
            test_index = np.zeros(_num_samples(X), dtype=np.bool)
            for l in unique_labels[np.array(indices)]:
                test_index[labels == l] = True
            yield test_index

    def get_n_splits(self, X, y, labels):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        if labels is None:
            raise ValueError("The labels parameter should not be None")
        return int(comb(len(np.unique(labels)), self.n_labels, exact=True))


class BaseShuffleSplit(with_metaclass(ABCMeta)):
    """Base class for ShuffleSplit and StratifiedShuffleSplit"""

    def __init__(self, n_iter=10, test_size=0.1, train_size=None,
                 random_state=None):
        _validate_shuffle_split_init(test_size, train_size)
        self.n_iter = n_iter
        self.test_size = test_size
        self.train_size = train_size
        self.random_state = random_state

    def split(self, X, y=None, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        X, y, labels = indexable(X, y, labels)
        for train, test in self._iter_indices(X, y, labels):
            yield train, test

    @abstractmethod
    def _iter_indices(self, X, y=None, labels=None):
        """Generate (train, test) indices"""

    def get_n_splits(self, X=None, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return self.n_iter

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
    n_iter : int (default 10)
        Number of re-shuffling & splitting iterations.

    test_size : float, int, or None, default 0.1
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
    >>> from sklearn.model_selection import ShuffleSplit
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 1, 2])
    >>> rs = ShuffleSplit(n_iter=3, test_size=.25, random_state=0)
    >>> rs.get_n_splits(X)
    3
    >>> print(rs)
    ShuffleSplit(n_iter=3, random_state=0, test_size=0.25, train_size=None)
    >>> for train_index, test_index in rs.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...  # doctest: +ELLIPSIS
    TRAIN: [3 1 0] TEST: [2]
    TRAIN: [2 1 3] TEST: [0]
    TRAIN: [0 2 1] TEST: [3]
    >>> rs = ShuffleSplit(n_iter=3, train_size=0.5, test_size=.25,
    ...                   random_state=0)
    >>> for train_index, test_index in rs.split(X):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...  # doctest: +ELLIPSIS
    TRAIN: [3 1] TEST: [2]
    TRAIN: [2 1] TEST: [0]
    TRAIN: [0 2] TEST: [3]
    """

    def _iter_indices(self, X, y=None, labels=None):
        n_samples = _num_samples(X)
        n_train, n_test = _validate_shuffle_split(n_samples, self.test_size,
                                                  self.train_size)
        rng = check_random_state(self.random_state)
        for i in range(self.n_iter):
            # random partition
            permutation = rng.permutation(n_samples)
            ind_test = permutation[:n_test]
            ind_train = permutation[n_test:(n_test + n_train)]
            yield ind_train, ind_test


class LabelShuffleSplit(ShuffleSplit):
    '''Shuffle-Labels-Out cross-validation iterator

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
    ``LeavePLabelOut(p=10)`` would be
    ``LabelShuffleSplit(test_size=10, n_iter=100)``.

    Note: The parameters ``test_size`` and ``train_size`` refer to labels, and
    not to samples, as in ShuffleSplit.


    Parameters
    ----------
    n_iter : int (default 5)
        Number of re-shuffling & splitting iterations.

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
    '''

    def __init__(self, n_iter=5, test_size=0.2, train_size=None,
                 random_state=None):
        super(LabelShuffleSplit, self).__init__(
            n_iter=n_iter,
            test_size=test_size,
            train_size=train_size,
            random_state=random_state)

    def _iter_indices(self, X, y, labels):
        if labels is None:
            raise ValueError("The labels parameter should not be None")
        classes, label_indices = np.unique(labels, return_inverse=True)
        for label_train, label_test in super(
                LabelShuffleSplit, self)._iter_indices(X=classes):
            # these are the indices of classes in the partition
            # invert them into data indices

            train = np.flatnonzero(np.in1d(label_indices, label_train))
            test = np.flatnonzero(np.in1d(label_indices, label_test))

            yield train, test


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
    >>> from sklearn.model_selection import StratifiedShuffleSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> sss = StratifiedShuffleSplit(n_iter=3, test_size=0.5, random_state=0)
    >>> sss.get_n_splits(X, y)
    3
    >>> print(sss)       # doctest: +ELLIPSIS
    StratifiedShuffleSplit(n_iter=3, random_state=0, ...)
    >>> for train_index, test_index in sss.split(X, y):
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 2] TEST: [3 0]
    TRAIN: [0 2] TEST: [1 3]
    TRAIN: [0 2] TEST: [3 1]
    """

    def __init__(self, n_iter=10, test_size=0.1, train_size=None,
                 random_state=None):
        super(StratifiedShuffleSplit, self).__init__(
            n_iter, test_size, train_size, random_state)

    def _iter_indices(self, X, y, labels=None):
        n_samples = _num_samples(X)
        n_train, n_test = _validate_shuffle_split(n_samples, self.test_size,
                                                  self.train_size)
        classes, y_indices = np.unique(y, return_inverse=True)
        n_classes = classes.shape[0]

        class_counts = bincount(y_indices)
        if np.min(class_counts) < 2:
            raise ValueError("The least populated class in y has only 1"
                             " member, which is too few. The minimum"
                             " number of labels for any class cannot"
                             " be less than 2.")

        if n_train < n_classes:
            raise ValueError('The train_size = %d should be greater or '
                             'equal to the number of classes = %d' %
                             (n_train, n_classes))
        if n_test < n_classes:
            raise ValueError('The test_size = %d should be greater or '
                             'equal to the number of classes = %d' %
                             (n_test, n_classes))

        rng = check_random_state(self.random_state)
        p_i = class_counts / float(n_samples)
        n_i = np.round(n_train * p_i).astype(int)
        t_i = np.minimum(class_counts - n_i,
                         np.round(n_test * p_i).astype(int))

        for _ in range(self.n_iter):
            train = []
            test = []

            for i, class_i in enumerate(classes):
                permutation = rng.permutation(class_counts[i])
                perm_indices_class_i = np.where((y == class_i))[0][permutation]

                train.extend(perm_indices_class_i[:n_i[i]])
                test.extend(perm_indices_class_i[n_i[i]:n_i[i] + t_i[i]])

            # Because of rounding issues (as n_train and n_test are not
            # dividers of the number of elements per class), we may end
            # up here with less samples in train and test than asked for.
            if len(train) + len(test) < n_train + n_test:
                # We complete by affecting randomly the missing indexes
                missing_indices = np.where(bincount(train + test,
                                                    minlength=len(y)) == 0)[0]
                missing_indices = rng.permutation(missing_indices)
                n_missing_train = n_train - len(train)
                n_missing_test = n_test - len(test)

                if n_missing_train > 0:
                    train.extend(missing_indices[:n_missing_train])
                if n_missing_test > 0:
                    test.extend(missing_indices[-n_missing_test:])

            train = rng.permutation(train)
            test = rng.permutation(test)

            yield train, test

    def split(self, X, y, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            The target variable for supervised learning problems.

        labels : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Returns
        -------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        return super(StratifiedShuffleSplit, self).split(X, y, labels)


def _validate_shuffle_split_init(test_size, train_size):
    """Validation helper to check the test_size and train_size at init

    NOTE This does not take into account the number of samples which is known
    only at split
    """
    if test_size is None and train_size is None:
        raise ValueError('test_size and train_size can not both be None')

    if test_size is not None:
        if np.asarray(test_size).dtype.kind == 'f':
            if test_size >= 1.:
                raise ValueError(
                    'test_size=%f should be smaller '
                    'than 1.0 or be an integer' % test_size)
        elif np.asarray(test_size).dtype.kind != 'i':
            # int values are checked during split based on the input
            raise ValueError("Invalid value for test_size: %r" % test_size)

    if train_size is not None:
        if np.asarray(train_size).dtype.kind == 'f':
            if train_size >= 1.:
                raise ValueError("train_size=%f should be smaller "
                                 "than 1.0 or be an integer" % train_size)
            elif (np.asarray(test_size).dtype.kind == 'f' and
                    (train_size + test_size) > 1.):
                raise ValueError('The sum of test_size and train_size = %f, '
                                 'should be smaller than 1.0. Reduce '
                                 'test_size and/or train_size.' %
                                 (train_size + test_size))
        elif np.asarray(train_size).dtype.kind != 'i':
            # int values are checked during split based on the input
            raise ValueError("Invalid value for train_size: %r" % train_size)


def _validate_shuffle_split(n_samples, test_size, train_size):
    """
    Validation helper to check if the test/test sizes are meaningful wrt to the
    size of the data (n_samples)
    """
    if (test_size is not None and np.asarray(test_size).dtype.kind == 'i'
            and test_size >= n_samples):
        raise ValueError('test_size=%d should be smaller than the number of '
                         'samples %d' % (test_size, n_samples))

    if (train_size is not None and np.asarray(train_size).dtype.kind == 'i'
            and train_size >= n_samples):
        raise ValueError("train_size=%d should be smaller than the number of"
                         " samples %d" % (train_size, n_samples))

    if np.asarray(test_size).dtype.kind == 'f':
        n_test = ceil(test_size * n_samples)
    elif np.asarray(test_size).dtype.kind == 'i':
        n_test = float(test_size)

    if train_size is None:
        n_train = n_samples - n_test
    elif np.asarray(train_size).dtype.kind == 'f':
        n_train = floor(train_size * n_samples)
    else:
        n_train = float(train_size)

    if test_size is None:
        n_test = n_samples - n_train

    if n_train + n_test > n_samples:
        raise ValueError('The sum of train_size and test_size = %d, '
                         'should be smaller than the number of '
                         'samples %d. Reduce test_size and/or '
                         'train_size.' % (n_train + n_test, n_samples))

    return int(n_train), int(n_test)


class PredefinedSplit(BaseCrossValidator):
    """Predefined split cross-validator

    Splits the data into training/test set folds according to a predefined
    scheme. Each sample can be assigned to at most one test set fold, as
    specified by the user through the ``test_fold`` parameter.

    Read more in the :ref:`User Guide <cross_validation>`.

    Examples
    --------
    >>> from sklearn.model_selection import PredefinedSplit
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> test_fold = [0, 1, -1, 1]
    >>> ps = PredefinedSplit(test_fold)
    >>> ps.get_n_splits()
    2
    >>> print(ps)       # doctest: +NORMALIZE_WHITESPACE +ELLIPSIS
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

    def split(self, X=None, y=None, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.

        Returns
        -------
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

    def get_n_splits(self, X=None, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : object
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
        self.cv = cv

    def get_n_splits(self, X=None, y=None, labels=None):
        """Returns the number of splitting iterations in the cross-validator

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.

        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return len(self.cv)  # Both iterables and old-cv objects support len

    def split(self, X=None, y=None, labels=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.

        y : object
            Always ignored, exists for compatibility.

        labels : object
            Always ignored, exists for compatibility.

        Returns
        -------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        for train, test in self.cv:
            yield train, test


def check_cv(cv=3, y=None, classifier=False):
    """Input checker utility for building a cross-validator

    Parameters
    ----------
    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
          - None, to use the default 3-fold cross-validation,
          - integer, to specify the number of folds.
          - An object to be used as a cross-validation generator.
          - An iterable yielding train/test splits.

        For integer/None inputs, if classifier is True and ``y`` is either
        binary or multiclass, :class:`StratifiedKFold` used. In all other
        cases, :class:`KFold` is used.

        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.

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
    if cv is None:
        cv = 3

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

        allowed inputs are lists, numpy arrays, scipy-sparse
        matrices or pandas dataframes.

        .. versionadded:: 0.16
            preserves input type instead of always casting to numpy array.

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

    Returns
    -------
    splitting : list, length=2 * len(arrays)
        List containing train-test split of inputs.

        .. versionadded:: 0.16
            Output type is the same as the input type.

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
        CVClass = StratifiedShuffleSplit
    else:
        CVClass = ShuffleSplit

    cv = CVClass(test_size=test_size,
                 train_size=train_size,
                 random_state=random_state)

    train, test = next(cv.split(X=arrays[0], y=stratify))
    return list(chain.from_iterable((safe_indexing(a, train),
                                     safe_indexing(a, test)) for a in arrays))


train_test_split.__test__ = False  # to avoid a pb with nosetests


def _safe_split(estimator, X, y, indices, train_indices=None):
    """Create subset of dataset and properly handle kernels."""
    if (hasattr(estimator, 'kernel') and callable(estimator.kernel) and
            not isinstance(estimator.kernel, GPKernel)):
        # cannot compute the kernel values with custom function
        raise ValueError("Cannot use a custom kernel function. "
                         "Precompute the kernel matrix instead.")

    if not hasattr(X, "shape"):
        if getattr(estimator, "_pairwise", False):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        X_subset = [X[index] for index in indices]
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
            if len(w) and w[0].category == DeprecationWarning:
                # if the parameter is deprecated, don't show it
                continue
        finally:
            warnings.filters.pop(0)
        params[key] = value

    return '%s(%s)' % (class_name, _pprint(params, offset=len(class_name)))
