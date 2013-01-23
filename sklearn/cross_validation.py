"""
The :mod:`sklearn.cross_validation` module includes utilities for cross-
validation and performance evaluation.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>,
#         Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD Style.

from __future__ import print_function

import warnings
from itertools import combinations
from math import ceil, floor, factorial
import numbers

import numpy as np
import scipy.sparse as sp

from .base import is_classifier, clone
from .utils import check_arrays, check_random_state, safe_mask
from .utils.fixes import unique
from .externals.joblib import Parallel, delayed

__all__ = ['Bootstrap',
           'KFold',
           'LeaveOneLabelOut',
           'LeaveOneOut',
           'LeavePLabelOut',
           'LeavePOut',
           'ShuffleSplit',
           'StratifiedKFold',
           'StratifiedShuffleSplit',
           'check_cv',
           'cross_val_score',
           'permutation_test_score',
           'train_test_split']


class LeaveOneOut(object):
    """Leave-One-Out cross validation iterator.

    Provides train/test indices to split data in train test sets. Each
    sample is used once as a test set (singleton) while the remaining
    samples form the training set.

    Due to the high number of test sets (which is the same as the
    number of samples) this cross validation method can be very costly.
    For large datasets one should favor KFold, StratifiedKFold or
    ShuffleSplit.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

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

    def __init__(self, n, indices=True):
        self.n = n
        self.indices = indices

    def __iter__(self):
        n = self.n
        for i in xrange(n):
            test_index = np.zeros(n, dtype=np.bool)
            test_index[i] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(n)
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

    def __repr__(self):
        return '%s.%s(n=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
        )

    def __len__(self):
        return self.n


class LeavePOut(object):
    """Leave-P-Out cross validation iterator

    Provides train/test indices to split data in train test sets. The
    test set is built using p samples while the remaining samples form
    the training set.

    Due to the high number of iterations which grows with the number of
    samples this cross validation method can be very costly. For large
    datasets one should favor KFold, StratifiedKFold or ShuffleSplit.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.

    p : int
        Size of the test sets.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

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

    def __init__(self, n, p, indices=True):
        self.n = n
        self.p = p
        self.indices = indices

    def __iter__(self):
        n = self.n
        p = self.p
        comb = combinations(range(n), p)
        for idx in comb:
            test_index = np.zeros(n, dtype=np.bool)
            test_index[np.array(idx)] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(n)
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

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


def _validate_kfold(k, n_samples):
    if k <= 0:
        raise ValueError("Cannot have number of folds k below 1.")
    if k > n_samples:
        raise ValueError("Cannot have number of folds k=%d greater than"
                         " the number of samples: %d." % (k, n_samples))


class KFold(object):
    """K-Folds cross validation iterator.

    Provides train/test indices to split data in train test sets. Split
    dataset into k consecutive folds (without shuffling).

    Each fold is then used a validation set once while the k - 1 remaining
    fold form the training set.

    Parameters
    ----------
    n : int
        Total number of elements.

    n_folds : int, default=3
        Number of folds.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    shuffle : boolean, optional
        Whether to shuffle the data before splitting into batches.

    random_state : int or RandomState
            Pseudo number generator state used for random sampling.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4])
    >>> kf = cross_validation.KFold(4, n_folds=2)
    >>> len(kf)
    2
    >>> print(kf)
    sklearn.cross_validation.KFold(n=4, n_folds=2)
    >>> for train_index, test_index in kf:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [2 3] TEST: [0 1]
    TRAIN: [0 1] TEST: [2 3]

    Notes
    -----
    All the folds have size trunc(n_samples / n_folds), the last one has the
    complementary.

    See also
    --------
    StratifiedKFold: take label information into account to avoid building
    folds with imbalanced class distributions (for binary or multiclass
    classification tasks).
    """

    def __init__(self, n, n_folds=3, indices=True, shuffle=False,
                 random_state=None, k=None):
        if k is not None:  # pragma: no cover
            warnings.warn("The parameter k was renamed to n_folds and will be"
                          " removed in 0.15.", DeprecationWarning)
            n_folds = k
        _validate_kfold(n_folds, n)
        random_state = check_random_state(random_state)

        if abs(n - int(n)) >= np.finfo('f').eps:
            raise ValueError("n must be an integer")
        self.n = int(n)
        if abs(n_folds - int(n_folds)) >= np.finfo('f').eps:
            raise ValueError("n_folds must be an integer")
        self.n_folds = int(n_folds)
        self.indices = indices
        self.idxs = np.arange(n)
        if shuffle:
            random_state.shuffle(self.idxs)

    def __iter__(self):
        n = self.n
        n_folds = self.n_folds
        fold_size = n // n_folds

        for i in xrange(n_folds):
            test_index = np.zeros(n, dtype=np.bool)
            if i < n_folds - 1:
                test_index[self.idxs[i * fold_size:(i + 1) * fold_size]] = True
            else:
                test_index[self.idxs[i * fold_size:]] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                train_index = self.idxs[train_index]
                test_index = self.idxs[test_index]
            yield train_index, test_index

    def __repr__(self):
        return '%s.%s(n=%i, n_folds=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
            self.n_folds,
        )

    def __len__(self):
        return self.n_folds


class StratifiedKFold(object):
    """Stratified K-Folds cross validation iterator

    Provides train/test indices to split data in train test sets.

    This cross-validation object is a variation of KFold, which
    returns stratified folds. The folds are made by preserving
    the percentage of samples for each class.

    Parameters
    ----------
    y : array-like, [n_samples]
        Samples to split in K folds.

    n_folds : int, default=3
        Number of folds.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> skf = cross_validation.StratifiedKFold(y, n_folds=2)
    >>> len(skf)
    2
    >>> print(skf)
    sklearn.cross_validation.StratifiedKFold(labels=[0 0 1 1], n_folds=2)
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
    """

    def __init__(self, y, n_folds=3, indices=True, k=None):
        if k is not None:  # pragma: no cover
            warnings.warn("The parameter k was renamed to n_folds and will be"
                          " removed in 0.15.", DeprecationWarning)
            n_folds = k
        y = np.asarray(y)
        n = y.shape[0]
        _validate_kfold(n_folds, n)
        _, y_sorted = unique(y, return_inverse=True)
        min_labels = np.min(np.bincount(y_sorted))
        if n_folds > min_labels:
            warnings.warn(("The least populated class in y has only %d"
                          " members, which is too few. The minimum"
                          " number of labels for any class cannot"
                          " be less than n_folds=%d."
                          % (min_labels, n_folds)), Warning)
        self.y = y
        self.n_folds = n_folds
        self.indices = indices

    def __iter__(self):
        y = self.y.copy()
        n_folds = self.n_folds
        n = y.size
        idx = np.argsort(y)

        for i in xrange(n_folds):
            test_index = np.zeros(n, dtype=np.bool)
            test_index[idx[i::n_folds]] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(n)
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

    def __repr__(self):
        return '%s.%s(labels=%s, n_folds=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.y,
            self.n_folds,
        )

    def __len__(self):
        return self.n_folds


class LeaveOneLabelOut(object):
    """Leave-One-Label_Out cross-validation iterator

    Provides train/test indices to split data according to a third-party
    provided label. This label information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the labels could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    Parameters
    ----------
    labels : array-like of int with shape (n_samples,)
        Arbitrary domain-specific stratification of the data to be used
        to draw the splits.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

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

    """

    def __init__(self, labels, indices=True):
        self.labels = labels
        self.n_unique_labels = unique(labels).size
        self.indices = indices

    def __iter__(self):
        # We make a copy here to avoid side-effects during iteration
        labels = np.array(self.labels, copy=True)
        for i in unique(labels):
            test_index = np.zeros(len(labels), dtype=np.bool)
            test_index[labels == i] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(len(labels))
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

    def __repr__(self):
        return '%s.%s(labels=%s)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.labels,
        )

    def __len__(self):
        return self.n_unique_labels


class LeavePLabelOut(object):
    """Leave-P-Label_Out cross-validation iterator

    Provides train/test indices to split data according to a third-party
    provided label. This label information can be used to encode arbitrary
    domain specific stratifications of the samples as integers.

    For instance the labels could be the year of collection of the samples
    and thus allow for cross-validation against time-based splits.

    The difference between LeavePLabelOut and LeaveOneLabelOut is that
    the former builds the test sets with all the samples assigned to
    ``p`` different values of the labels while the latter uses samples
    all assigned the same labels.

    Parameters
    ----------
    labels : array-like of int with shape (n_samples,)
        Arbitrary domain-specific stratification of the data to be used
        to draw the splits.

    p : int
        Number of samples to leave out in the test split.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

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
    """

    def __init__(self, labels, p, indices=True):
        self.labels = labels
        self.unique_labels = unique(self.labels)
        self.n_unique_labels = self.unique_labels.size
        self.p = p
        self.indices = indices

    def __iter__(self):
        # We make a copy here to avoid side-effects during iteration
        labels = np.array(self.labels, copy=True)
        unique_labels = unique(labels)
        comb = combinations(range(self.n_unique_labels), self.p)

        for idx in comb:
            test_index = np.zeros(labels.size, dtype=np.bool)
            idx = np.array(idx)
            for l in unique_labels[idx]:
                test_index[labels == l] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(labels.size)
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

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


class Bootstrap(object):
    """Random sampling with replacement cross-validation iterator

    Provides train/test indices to split data in train test sets
    while resampling the input n_iter times: each time a new
    random split of the data is performed and then samples are drawn
    (with replacement) on each side of the split to build the training
    and test sets.

    Note: contrary to other cross-validation strategies, bootstrapping
    will allow some samples to occur several times in each splits. However
    a sample that occurs in the train split will never occur in the test
    split and vice-versa.

    If you want each sample to occur at most once you should probably
    use ShuffleSplit cross validation instead.

    Parameters
    ----------
    n : int
        Total number of elements in the dataset.

    n_iter : int (default is 3)
        Number of bootstrapping iterations

    train_size : int or float (default is 0.5)
        If int, number of samples to include in the training split
        (should be smaller than the total number of samples passed
        in the dataset).

        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split.

    test_size : int or float or None (default is None)
        If int, number of samples to include in the training set
        (should be smaller than the total number of samples passed
        in the dataset).

        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the test split.

        If None, n_test is set as the complement of n_train.

    random_state : int or RandomState
        Pseudo number generator state used for random sampling.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> bs = cross_validation.Bootstrap(9, random_state=0)
    >>> len(bs)
    3
    >>> print(bs)
    Bootstrap(9, n_iter=3, train_size=5, test_size=4, random_state=0)
    >>> for train_index, test_index in bs:
    ...    print("TRAIN:", train_index, "TEST:", test_index)
    ...
    TRAIN: [1 8 7 7 8] TEST: [0 3 0 5]
    TRAIN: [5 4 2 4 2] TEST: [6 7 1 0]
    TRAIN: [4 7 0 1 1] TEST: [5 3 6 5]

    See also
    --------
    ShuffleSplit: cross validation using random permutations.
    """

    # Static marker to be able to introspect the CV type
    indices = True

    def __init__(self, n, n_iter=3, train_size=.5, test_size=None,
                 random_state=None, n_bootstraps=None):
        self.n = n
        if n_bootstraps is not None:  # pragma: no cover
            warnings.warn("n_bootstraps was renamed to n_iter and will "
                          "be removed in 0.16.", DeprecationWarning)
            n_iter = n_bootstraps
        self.n_iter = n_iter
        if (isinstance(train_size, numbers.Real) and train_size >= 0.0
                and train_size <= 1.0):
            self.train_size = int(ceil(train_size * n))
        elif isinstance(train_size, numbers.Integral):
            self.train_size = train_size
        else:
            raise ValueError("Invalid value for train_size: %r" %
                             train_size)
        if self.train_size > n:
            raise ValueError("train_size=%d should not be larger than n=%d" %
                             (self.train_size, n))

        if isinstance(test_size, numbers.Real) and 0.0 <= test_size <= 1.0:
            self.test_size = int(ceil(test_size * n))
        elif isinstance(test_size, numbers.Integral):
            self.test_size = test_size
        elif test_size is None:
            self.test_size = self.n - self.train_size
        else:
            raise ValueError("Invalid value for test_size: %r" % test_size)
        if self.test_size > n:
            raise ValueError("test_size=%d should not be larger than n=%d" %
                             (self.test_size, n))

        self.random_state = random_state

    def __iter__(self):
        rng = check_random_state(self.random_state)
        for i in range(self.n_iter):
            # random partition
            permutation = rng.permutation(self.n)
            ind_train = permutation[:self.train_size]
            ind_test = permutation[self.train_size:self.train_size
                                   + self.test_size]

            # bootstrap in each split individually
            train = rng.randint(0, self.train_size,
                                size=(self.train_size,))
            test = rng.randint(0, self.test_size,
                               size=(self.test_size,))
            yield ind_train[train], ind_test[test]

    def __repr__(self):
        return ('%s(%d, n_iter=%d, train_size=%d, test_size=%d, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.n,
                    self.n_iter,
                    self.train_size,
                    self.test_size,
                    self.random_state,
                ))

    def __len__(self):
        return self.n_iter


class ShuffleSplit(object):
    """Random permutation cross-validation iterator.

    Yields indices to split data into training and test sets.

    Note: contrary to other cross-validation strategies, random splits
    do not guarantee that all folds will be different, although this is
    still very likely for sizeable datasets.

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

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

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
    ShuffleSplit(4, n_iter=3, test_size=0.25, indices=True, ...)
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

    See also
    --------
    Bootstrap: cross-validation using re-sampling with replacement.
    """

    def __init__(self, n, n_iter=10, test_size=0.1, train_size=None,
                 indices=True, random_state=None, n_iterations=None):
        self.n = n
        self.n_iter = n_iter
        if n_iterations is not None:  # pragma: no cover
            warnings.warn("n_iterations was renamed to n_iter for consistency "
                          " and will be removed in 0.16.")
            self.n_iter = n_iterations
        self.test_size = test_size
        self.train_size = train_size
        self.random_state = random_state
        self.indices = indices

        self.n_train, self.n_test = _validate_shuffle_split(n,
                                                            test_size,
                                                            train_size)

    def __iter__(self):
        rng = check_random_state(self.random_state)
        for i in range(self.n_iter):
            # random partition
            permutation = rng.permutation(self.n)
            ind_test = permutation[:self.n_test]
            ind_train = permutation[self.n_test:self.n_test + self.n_train]

            if self.indices:
                yield ind_train, ind_test
            else:
                train_mask = np.zeros(self.n, dtype=np.bool)
                train_mask[ind_train] = True
                test_mask = np.zeros(self.n, dtype=np.bool)
                test_mask[ind_test] = True
                yield train_mask, test_mask

    def __repr__(self):
        return ('%s(%d, n_iter=%d, test_size=%s, indices=%s, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.n,
                    self.n_iter,
                    str(self.test_size),
                    self.indices,
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


def _validate_stratified_shuffle_split(y, test_size, train_size):
    classes, y = unique(y, return_inverse=True)
    n_cls = classes.shape[0]

    if np.min(np.bincount(y)) < 2:
        raise ValueError("The least populated class in y has only 1"
                         " member, which is too few. The minimum"
                         " number of labels for any class cannot"
                         " be less than 2.")

    n_train, n_test = _validate_shuffle_split(y.size, test_size, train_size)

    if n_train < n_cls:
        raise ValueError('The train_size = %d should be greater or '
                         'equal to the number of classes = %d' %
                         (n_train, n_cls))
    if n_test < n_cls:
        raise ValueError('The test_size = %d should be greater or '
                         'equal to the number of classes = %d' %
                         (n_test, n_cls))

    return n_train, n_test, classes, y


class StratifiedShuffleSplit(object):
    """Stratified ShuffleSplit cross validation iterator

    Provides train/test indices to split data in train test sets.

    This cross-validation object is a merge of StratifiedKFold and
    ShuffleSplit, which returns stratified randomized folds. The folds
    are made by preserving the percentage of samples for each class.

    Note: like the ShuffleSplit strategy, stratified random splits
    do not guarantee that all folds will be different, although this is
    still very likely for sizeable datasets.

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

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

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
                 indices=True, random_state=None, n_iterations=None):

        self.y = np.array(y)
        self.n = self.y.size
        self.n_iter = n_iter
        if n_iterations is not None:  # pragma: no cover
            warnings.warn("n_iterations was renamed to n_iter for consistency"
                          " and will be removed in 0.16.")
            self.n_iter = n_iterations
        self.test_size = test_size
        self.train_size = train_size
        self.random_state = random_state
        self.indices = indices
        self.n_train, self.n_test, self.classes, self.y_indices = \
            _validate_stratified_shuffle_split(y, test_size, train_size)

    def __iter__(self):
        rng = check_random_state(self.random_state)
        cls_count = np.bincount(self.y_indices)
        p_i = cls_count / float(self.n)
        n_i = np.round(self.n_train * p_i).astype(int)
        t_i = np.minimum(cls_count - n_i,
                         np.round(self.n_test * p_i).astype(int))

        for n in range(self.n_iter):
            train = []
            test = []

            for i, cls in enumerate(self.classes):
                permutation = rng.permutation(n_i[i] + t_i[i])
                cls_i = np.where((self.y == cls))[0][permutation]

                train.extend(cls_i[:n_i[i]])
                test.extend(cls_i[n_i[i]:n_i[i] + t_i[i]])

            train = rng.permutation(train)
            test = rng.permutation(test)

            if self.indices:
                yield train, test
            else:
                train_m = np.zeros(self.n, dtype=bool)
                test_m = np.zeros(self.n, dtype=bool)
                train_m[train] = True
                test_m[test] = True

                yield train_m, test_m

    def __repr__(self):
        return ('%s(labels=%s, n_iter=%d, test_size=%s, indices=%s, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.y,
                    self.n_iter,
                    str(self.test_size),
                    self.indices,
                    self.random_state,
                ))

    def __len__(self):
        return self.n_iter


##############################################################################

def _cross_val_score(estimator, X, y, score_func, train, test, verbose,
                     fit_params):
    """Inner loop for cross validation"""
    n_samples = X.shape[0] if sp.issparse(X) else len(X)
    fit_params = dict([(k, np.asarray(v)[train]
                        if hasattr(v, '__len__')
                        and len(v) == n_samples else v)
                       for k, v in fit_params.items()])
    if not hasattr(X, "shape"):
        if getattr(estimator, "_pairwise", False):
            raise ValueError("Precomputed kernels or affinity matrices have "
                             "to be passed as arrays or sparse matrices.")
        X_train = [X[idx] for idx in train]
        X_test = [X[idx] for idx in test]
    else:
        if getattr(estimator, "_pairwise", False):
            # X is a precomputed square kernel matrix
            if X.shape[0] != X.shape[1]:
                raise ValueError("X should be a square kernel matrix")
            X_train = X[np.ix_(train, train)]
            X_test = X[np.ix_(test, train)]
        else:
            X_train = X[safe_mask(X, train)]
            X_test = X[safe_mask(X, test)]

    if y is None:
        estimator.fit(X_train, **fit_params)
        if score_func is None:
            score = estimator.score(X_test)
        else:
            score = score_func(X_test)
    else:
        estimator.fit(X_train, y[train], **fit_params)
        if score_func is None:
            score = estimator.score(X_test, y[test])
        else:
            score = score_func(y[test], estimator.predict(X_test))
    if verbose > 1:
        print("score: %f" % score)
    return score


def cross_val_score(estimator, X, y=None, score_func=None, cv=None, n_jobs=1,
                    verbose=0, fit_params=None):
    """Evaluate a score by cross-validation

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like, optional
        The target variable to try to predict in the case of
        supervised learning.

    score_func : callable, optional
        Score function to use for evaluation.
        Has priority over the score function in the estimator.
        In a non-supervised setting, where y is None, it takes the test
        data (X_test) as its only argument. In a supervised setting it takes
        the test target (y_true) and the test prediction (y_pred) as arguments.

    cv : cross-validation generator, optional
        A cross-validation generator. If None, a 3-fold cross
        validation is used or 3-fold stratified cross-validation
        when y is supplied and estimator is a classifier.

    n_jobs : integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    verbose : integer, optional
        The verbosity level.

    fit_params : dict, optional
        Parameters to pass to the fit method of the estimator.

    Returns
    -------
    scores : array of float, shape=(len(list(cv)),)
        Array of scores of the estimator for each run of the cross validation.
    """
    X, y = check_arrays(X, y, sparse_format='csr', allow_lists=True)
    cv = check_cv(cv, X, y, classifier=is_classifier(estimator))
    if score_func is None:
        if not hasattr(estimator, 'score'):
            raise TypeError(
                "If no score_func is specified, the estimator passed "
                "should have a 'score' method. The estimator %s "
                "does not." % estimator)
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    fit_params = fit_params if fit_params is not None else {}
    scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_cross_val_score)(
            clone(estimator), X, y, score_func,
            train, test, verbose, fit_params)
        for train, test in cv)
    return np.array(scores)


def _permutation_test_score(estimator, X, y, cv, score_func):
    """Auxilary function for permutation_test_score"""
    avg_score = []
    for train, test in cv:
        avg_score.append(score_func(y[test],
                                    estimator.fit(X[train],
                                                  y[train]).predict(X[test])))
    return np.mean(avg_score)


def _shuffle(y, labels, random_state):
    """Return a shuffled copy of y eventually shuffle among same labels."""
    if labels is None:
        ind = random_state.permutation(y.size)
    else:
        ind = np.arange(labels.size)
        for label in unique(labels):
            this_mask = (labels == label)
            ind[this_mask] = random_state.permutation(ind[this_mask])
    return y[ind]


def check_cv(cv, X=None, y=None, classifier=False):
    """Input checker utility for building a CV in a user friendly way.

    Parameters
    ----------
    cv : int, a cv generator instance, or None
        The input specifying which cv generator to use. It can be an
        integer, in which case it is the number of folds in a KFold,
        None, in which case 3 fold is used, or another object, that
        will then be used as a cv generator.

    X : array-like
        The data the cross-val object will be applied on.

    y : array-like
        The target variable for a supervised learning problem.

    classifier : boolean optional
        Whether the task is a classification task, in which case
        stratified KFold will be used.
    """
    is_sparse = sp.issparse(X)
    needs_indices = is_sparse or not hasattr(X, "shape")
    if cv is None:
        cv = 3
    if isinstance(cv, numbers.Integral):
        if classifier:
            cv = StratifiedKFold(y, cv, indices=needs_indices)
        else:
            if not is_sparse:
                n_samples = len(X)
            else:
                n_samples = X.shape[0]
            cv = KFold(n_samples, cv, indices=needs_indices)
    if needs_indices and not getattr(cv, "indices", True):
        raise ValueError("Sparse data and lists require indices-based cross"
                         " validation generator, got: %r", cv)
    return cv


def permutation_test_score(estimator, X, y, score_func, cv=None,
                           n_permutations=100, n_jobs=1, labels=None,
                           random_state=0, verbose=0):
    """Evaluate the significance of a cross-validated score with permutations

    Parameters
    ----------
    estimator : estimator object implementing 'fit'
        The object to use to fit the data.

    X : array-like of shape at least 2D
        The data to fit.

    y : array-like
        The target variable to try to predict in the case of
        supervised learning.

    score_func : callable
        Callable taking as arguments the test targets (y_test) and
        the predicted targets (y_pred) and returns a float. The score
        functions are expected to return a bigger value for a better result
        otherwise the returned value does not correspond to a p-value (see
        Returns below for further details).

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects.

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

    permutation_scores : array, shape = [n_permutations]
        The scores obtained for each permutations.

    pvalue : float
        The returned value equals p-value if `score_func` returns bigger
        numbers for better scores (e.g., accuracy_score). If `score_func` is
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
    X, y = check_arrays(X, y, sparse_format='csr')
    cv = check_cv(cv, X, y, classifier=is_classifier(estimator))

    random_state = check_random_state(random_state)

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    score = _permutation_test_score(clone(estimator), X, y, cv, score_func)
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(_permutation_test_score)(clone(estimator), X,
                                         _shuffle(y, labels, random_state),
                                         cv, score_func)
        for _ in range(n_permutations))
    permutation_scores = np.array(permutation_scores)
    pvalue = (np.sum(permutation_scores >= score) + 1.0) / (n_permutations + 1)
    return score, permutation_scores, pvalue


permutation_test_score.__test__ = False  # to avoid a pb with nosetests


def train_test_split(*arrays, **options):
    """Split arrays or matrices into random train and test subsets

    Quick utility that wraps calls to ``check_arrays`` and
    ``next(iter(ShuffleSplit(n_samples)))`` and application to input
    data into a single call for splitting (and optionally subsampling)
    data in a oneliner.

    Parameters
    ----------
    *arrays : sequence of arrays or scipy.sparse matrices with same shape[0]
        Python lists or tuples occurring in arrays are converted to 1D numpy
        arrays.

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

    dtype : a numpy dtype instance, None by default
        Enforce a specific dtype.

    Returns
    -------
    splitting : list of arrays, length=2 * len(arrays)
        List containing train-test split of input array.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.cross_validation import train_test_split
    >>> a, b = np.arange(10).reshape((5, 2)), range(5)
    >>> a
    array([[0, 1],
           [2, 3],
           [4, 5],
           [6, 7],
           [8, 9]])
    >>> list(b)
    [0, 1, 2, 3, 4]

    >>> a_train, a_test, b_train, b_test = train_test_split(
    ...     a, b, test_size=0.33, random_state=42)
    ...
    >>> a_train
    array([[4, 5],
           [0, 1],
           [6, 7]])
    >>> b_train
    array([2, 0, 3])
    >>> a_test
    array([[2, 3],
           [8, 9]])
    >>> b_test
    array([1, 4])

    """
    n_arrays = len(arrays)
    if n_arrays == 0:
        raise ValueError("At least one array required as input")

    test_size = options.pop('test_size', None)
    train_size = options.pop('train_size', None)
    random_state = options.pop('random_state', None)
    options['sparse_format'] = 'csr'

    if test_size is None and train_size is None:
        test_size = 0.25

    arrays = check_arrays(*arrays, **options)
    n_samples = arrays[0].shape[0]
    cv = ShuffleSplit(n_samples, test_size=test_size,
                      train_size=train_size,
                      random_state=random_state,
                      indices=True)
    train, test = next(iter(cv))
    splitted = []
    for a in arrays:
        splitted.append(a[train])
        splitted.append(a[test])
    return splitted

train_test_split.__test__ = False  # to avoid a pb with nosetests
