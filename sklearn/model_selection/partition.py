"""
The :mod:`sklearn.model_selection.partition` module includes
"""
#TODO Complete docstring

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
from abc import ABCMeta, abstractmethod

import numpy as np
import scipy.sparse as sp

from sklearn.utils import check_arrays, check_random_state
from sklearn.externals.six import with_metaclass
from sklearn.externals.six.moves import zip


class _PartitionIterator(with_metaclass(ABCMeta)):
    """Base class for CV iterators where train_mask = ~test_mask

    Implementations must define `_iter_test_masks` or `_iter_test_indices`.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.
    """

    def __init__(self, n, indices=None):
        if indices is None:
            indices = True
        else:
            warnings.warn("The indices parameter is deprecated and will be "
                          "removed (assumed True) in 0.17", DeprecationWarning,
                          stacklevel=1)
        if abs(n - int(n)) >= np.finfo('f').eps:
            raise ValueError("n must be an integer")
        self.n = int(n)
        self._indices = indices

    @property
    def indices(self):
        warnings.warn("The indices attribute is deprecated and will be "
                      "removed (assumed True) in 0.17", DeprecationWarning,
                      stacklevel=1)
        return self._indices

    def __iter__(self):
        indices = self._indices
        if indices:
            ind = np.arange(self.n)
        for test_index in self._iter_test_masks():
            train_index = np.logical_not(test_index)
            if indices:
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

    Provides train/test indices to split data in train test sets. Each
    sample is used once as a test set (singleton) while the remaining
    samples form the training set.

    Note: ``LeaveOneOut(n)`` is equivalent to ``KFold(n, n_folds=n)`` and
    ``LeavePOut(n, p=1)``.

    Due to the high number of test sets (which is the same as the
    number of samples) this cross validation method can be very costly.
    For large datasets one should favor KFold, StratifiedKFold or
    ShuffleSplit.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.

    Examples
    --------
    >>> from sklearn import model_selection
    >>> X = np.array([[1, 2], [3, 4]])
    >>> y = np.array([1, 2])
    >>> loo = model_selection.LeaveOneOut(2)
    >>> len(loo)
    2
    >>> print(loo)
    sklearn.model_selection.partition.LeaveOneOut(n=2)
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

    Provides train/test indices to split data in train test sets. This results
    in testing on all distinct samples of size p, while the remaining n - p
    samples form the training set in each iteration.

    Note: ``LeavePOut(n, p)`` is NOT equivalent to ``KFold(n, n_folds=n // p)``
    which creates non-overlapping test sets.

    Due to the high number of iterations which grows combinatorically with the
    number of samples this cross validation method can be very costly. For
    large datasets one should favor KFold, StratifiedKFold or ShuffleSplit.

    Parameters
    ----------
    n : int
        Total number of elements in dataset.

    p : int
        Size of the test sets.

    Examples
    --------
    >>> from sklearn import model_selection
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 3, 4])
    >>> lpo = model_selection.LeavePOut(4, 2)
    >>> len(lpo)
    6
    >>> print(lpo)
    sklearn.model_selection.partition.LeavePOut(n=4, p=2)
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

    def __init__(self, n, p, indices=None):
        super(LeavePOut, self).__init__(n, indices)
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
    def __init__(self, n, n_folds, indices, shuffle, random_state):
        super(_BaseKFold, self).__init__(n, indices)

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

    Provides train/test indices to split data in train test sets. Split
    dataset into k consecutive folds (without shuffling).

    Each fold is then used a validation set once while the k - 1 remaining
    fold form the training set.

    Parameters
    ----------
    n : int
        Total number of elements.

    n_folds : int, default=3
        Number of folds. Must be at least 2.

    shuffle : boolean, optional
        Whether to shuffle the data before splitting into batches.

    random_state : None, int or RandomState
        Pseudo-random number generator state used for random
        sampling. If None, use default numpy RNG for shuffling

    Examples
    --------
    >>> from sklearn import model_selection
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4])
    >>> kf = model_selection.KFold(4, n_folds=2)
    >>> len(kf)
    2
    >>> print(kf)  # doctest: +NORMALIZE_WHITESPACE
    sklearn.model_selection.partition.KFold(n=4, n_folds=2, shuffle=False,
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
    StratifiedKFold: take label information into account to avoid building
    folds with imbalanced class distributions (for binary or multiclass
    classification tasks).
    """

    def __init__(self, n, n_folds=3, indices=None, shuffle=False,
                 random_state=None):
        super(KFold, self).__init__(n, n_folds, indices, shuffle, random_state)
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


class StratifiedKFold(_BaseKFold):
    """Stratified K-Folds cross validation iterator

    Provides train/test indices to split data in train test sets.

    This cross-validation object is a variation of KFold that
    returns stratified folds. The folds are made by preserving
    the percentage of samples for each class.

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
        Pseudo-random number generator state used for random
        sampling. If None, use default numpy RNG for shuffling

    Examples
    --------
    >>> from sklearn import model_selection
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> skf = model_selection.StratifiedKFold(y, n_folds=2)
    >>> len(skf)
    2
    >>> print(skf)  # doctest: +NORMALIZE_WHITESPACE
    sklearn.model_selection.partition.StratifiedKFold(labels=[0 0 1 1], n_folds=2,
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

    """

    def __init__(self, y, n_folds=3, indices=None, shuffle=False,
                 random_state=None):
        super(StratifiedKFold, self).__init__(
            len(y), n_folds, indices, shuffle, random_state)
        y = np.asarray(y)
        n_samples = y.shape[0]
        unique_labels, y_inversed = np.unique(y, return_inverse=True)
        label_counts = np.bincount(y_inversed)
        min_labels = np.min(label_counts)
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

    Examples
    --------
    >>> from sklearn import model_selection
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 1, 2])
    >>> labels = np.array([1, 1, 2, 2])
    >>> lol = model_selection.LeaveOneLabelOut(labels)
    >>> len(lol)
    2
    >>> print(lol)
    sklearn.model_selection.partition.LeaveOneLabelOut(labels=[1 1 2 2])
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

    def __init__(self, labels, indices=None):
        super(LeaveOneLabelOut, self).__init__(len(labels), indices)
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

    Examples
    --------
    >>> from sklearn import model_selection
    >>> X = np.array([[1, 2], [3, 4], [5, 6]])
    >>> y = np.array([1, 2, 1])
    >>> labels = np.array([1, 2, 3])
    >>> lpl = model_selection.LeavePLabelOut(labels, p=2)
    >>> len(lpl)
    3
    >>> print(lpl)
    sklearn.model_selection.partition.LeavePLabelOut(labels=[1 2 3], p=2)
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

    def __init__(self, labels, p, indices=None):
        # We make a copy of labels to avoid side-effects during iteration
        super(LeavePLabelOut, self).__init__(len(labels), indices)
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
    >>> from sklearn.model_selection.partition import train_test_split
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
    options['allow_nans'] = True

    if test_size is None and train_size is None:
        test_size = 0.25

    arrays = check_arrays(*arrays, **options)
    n_samples = arrays[0].shape[0]
    cv = ShuffleSplit(n_samples, test_size=test_size,
                      train_size=train_size,
                      random_state=random_state)

    train, test = next(iter(cv))
    return list(chain.from_iterable((a[train], a[test]) for a in arrays))


train_test_split.__test__ = False  # to avoid a pb with nosetests


class BaseShuffleSplit(with_metaclass(ABCMeta)):
    """Base class for ShuffleSplit and StratifiedShuffleSplit"""

    def __init__(self, n, n_iter=10, test_size=0.1, train_size=None,
                 indices=None, random_state=None, n_iterations=None):
        if indices is None:
            indices = True
        else:
            warnings.warn("The indices parameter is deprecated and will be "
                          "removed (assumed True) in 0.17", DeprecationWarning)
        self.n = n
        self.n_iter = n_iter
        if n_iterations is not None:  # pragma: no cover
            warnings.warn("n_iterations was renamed to n_iter for consistency "
                          " and will be removed in 0.16.")
            self.n_iter = n_iterations
        self.test_size = test_size
        self.train_size = train_size
        self.random_state = random_state
        self._indices = indices
        self.n_train, self.n_test = _validate_shuffle_split(n,
                                                            test_size,
                                                            train_size)

    @property
    def indices(self):
        warnings.warn("The indices attribute is deprecated and will be "
                      "removed (assumed True) in 0.17", DeprecationWarning,
                      stacklevel=1)
        return self._indices

    def __iter__(self):
        if self._indices:
            for train, test in self._iter_indices():
                yield train, test
            return
        for train, test in self._iter_indices():
            train_m = np.zeros(self.n, dtype=bool)
            test_m = np.zeros(self.n, dtype=bool)
            train_m[train] = True
            test_m[test] = True
            yield train_m, test_m

    @abstractmethod
    def _iter_indices(self):
        """Generate (train, test) indices"""


class ShuffleSplit(BaseShuffleSplit):
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

    See also
    --------
    Bootstrap: cross-validation using re-sampling with replacement.
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


class StratifiedShuffleSplit(BaseShuffleSplit):
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
                 indices=None, random_state=None, n_iterations=None):

        super(StratifiedShuffleSplit, self).__init__(
            len(y), n_iter, test_size, train_size, indices, random_state,
            n_iterations)
        self.y = np.array(y)
        self.classes, self.y_indices = np.unique(y, return_inverse=True)
        n_cls = self.classes.shape[0]

        if np.min(np.bincount(self.y_indices)) < 2:
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
        cls_count = np.bincount(self.y_indices)
        p_i = cls_count / float(self.n)
        n_i = np.round(self.n_train * p_i).astype(int)
        t_i = np.minimum(cls_count - n_i,
                         np.round(self.n_test * p_i).astype(int))

        for n in range(self.n_iter):
            train = []
            test = []

            for i, cls in enumerate(self.classes):
                permutation = rng.permutation(cls_count[i])
                cls_i = np.where((self.y == cls))[0][permutation]

                train.extend(cls_i[:n_i[i]])
                test.extend(cls_i[n_i[i]:n_i[i] + t_i[i]])

            # Because of rounding issues (as n_train and n_test are not
            # dividers of the number of elements per class), we may end
            # up here with less samples in train and test than asked for.
            if len(train) < self.n_train or len(test) < self.n_test:
                # We complete by affecting randomly the missing indexes
                missing_idx = np.where(np.bincount(train + test,
                                                   minlength=len(self.y)) == 0,
                                       )[0]
                missing_idx = rng.permutation(missing_idx)
                train.extend(missing_idx[:(self.n_train - len(train))])
                test.extend(missing_idx[-(self.n_test - len(test)):])

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

    Returns
    -------
    checked_cv: a cross-validation generator instance.
        The return value is guaranteed to be a cv generator instance, whatever
        the input type.
    """
    return _check_cv(cv, X=X, y=y, classifier=classifier, warn_mask=True)


def _check_cv(cv, X=None, y=None, classifier=False, warn_mask=False):
    # This exists for internal use while indices is being deprecated.
    is_sparse = sp.issparse(X)
    needs_indices = is_sparse or not hasattr(X, "shape")
    if cv is None:
        cv = 3
    if isinstance(cv, numbers.Integral):
        if warn_mask and not needs_indices:
            warnings.warn('check_cv will return indices instead of boolean '
                          'masks from 0.17', DeprecationWarning)
        else:
            needs_indices = None
        if classifier:
            cv = StratifiedKFold(y, cv, indices=needs_indices)
        else:
            if not is_sparse:
                n_samples = len(X)
            else:
                n_samples = X.shape[0]
            cv = KFold(n_samples, cv, indices=needs_indices)
    if needs_indices and not getattr(cv, "_indices", True):
        raise ValueError("Sparse data and lists require indices-based cross"
                         " validation generator, got: %r", cv)
    return cv
