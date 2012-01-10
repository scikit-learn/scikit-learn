"""
The :mod:`sklearn.cross_validation` module includes utilities for cross-
validation and performance evaluation.
"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux <gael.varoquaux@normalesup.org>,
#         Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD Style.

from itertools import combinations
from math import ceil, floor, factorial
import operator

import numpy as np
import scipy.sparse as sp

from .base import is_classifier, clone
from .utils import check_arrays, check_random_state
from .utils.fixes import unique
from .externals.joblib import Parallel, delayed


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
    n: int
        Total number of elements

    indices: boolean, optional (default True)
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
    >>> print loo
    sklearn.cross_validation.LeaveOneOut(n=2)
    >>> for train_index, test_index in loo:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    ...    print X_train, X_test, y_train, y_test
    TRAIN: [1] TEST: [0]
    [[3 4]] [[1 2]] [2] [1]
    TRAIN: [0] TEST: [1]
    [[1 2]] [[3 4]] [1] [2]

    See also
    ========
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
    n: int
        Total number of elements

    p: int
        Size of the test sets

    indices: boolean, optional (default True)
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
    >>> print lpo
    sklearn.cross_validation.LeavePOut(n=4, p=2)
    >>> for train_index, test_index in lpo:
    ...    print "TRAIN:", train_index, "TEST:", test_index
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
        return (factorial(self.n) / factorial(self.n - self.p)
                / factorial(self.p))


class KFold(object):
    """K-Folds cross validation iterator

    Provides train/test indices to split data in train test sets. Split
    dataset into k consecutive folds (without shuffling).

    Each fold is then used a validation set once while the k - 1 remaining
    fold form the training set.

    Parameters
    ----------
    n: int
        Total number of elements

    k: int
        Number of folds

    indices: boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([1, 2, 3, 4])
    >>> kf = cross_validation.KFold(4, k=2)
    >>> len(kf)
    2
    >>> print kf
    sklearn.cross_validation.KFold(n=4, k=2)
    >>> for train_index, test_index in kf:
    ...    print "TRAIN:", train_index, "TEST:", test_index
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

    def __init__(self, n, k, indices=True):
        assert k > 0, ValueError('Cannot have number of folds k below 1.')
        assert k <= n, ValueError('Cannot have number of folds k=%d, '
                                  'greater than the number '
                                  'of samples: %d.' % (k, n))
        self.n = n
        self.k = k
        self.indices = indices

    def __iter__(self):
        n = self.n
        k = self.k
        j = ceil(n / k)

        for i in xrange(k):
            test_index = np.zeros(n, dtype=np.bool)
            if i < k - 1:
                test_index[i * j:(i + 1) * j] = True
            else:
                test_index[i * j:] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(n)
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

    def __repr__(self):
        return '%s.%s(n=%i, k=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.n,
            self.k,
        )

    def __len__(self):
        return self.k


class StratifiedKFold(object):
    """Stratified K-Folds cross validation iterator

    Provides train/test indices to split data in train test sets.

    This cross-validation object is a variation of KFold, which
    returns stratified folds. The folds are made by preserving
    the percentage of samples for each class.

    Parameters
    ----------
    y: array, [n_samples]
        Samples to split in K folds

    k: int
        Number of folds

    indices: boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
    >>> y = np.array([0, 0, 1, 1])
    >>> skf = cross_validation.StratifiedKFold(y, k=2)
    >>> len(skf)
    2
    >>> print skf
    sklearn.cross_validation.StratifiedKFold(labels=[0 0 1 1], k=2)
    >>> for train_index, test_index in skf:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    TRAIN: [1 3] TEST: [0 2]
    TRAIN: [0 2] TEST: [1 3]

    Notes
    -----
    All the folds have size trunc(n_samples / n_folds), the last one has the
    complementary.
    """

    def __init__(self, y, k, indices=True):
        y = np.asarray(y)
        n = y.shape[0]
        assert k > 0, ValueError('Cannot have number of folds k below 1.')
        assert k <= n, ValueError('Cannot have number of folds k=%d, '
                                  'greater than the number '
                                  'of samples: %d.' % (k, n))
        _, y_sorted = unique(y, return_inverse=True)
        min_labels = np.min(np.bincount(y_sorted))
        assert k <= min_labels, ValueError(
            'Cannot have number of folds k=%d, smaller than %d, the minimum '
            'number of labels for any class.' % (k, min_labels))
        self.y = y
        self.k = k
        self.indices = indices

    def __iter__(self):
        y = self.y.copy()
        k = self.k
        n = y.size
        idx = np.argsort(y)

        for i in xrange(k):
            test_index = np.zeros(n, dtype=np.bool)
            test_index[idx[i::k]] = True
            train_index = np.logical_not(test_index)
            if self.indices:
                ind = np.arange(n)
                train_index = ind[train_index]
                test_index = ind[test_index]
            yield train_index, test_index

    def __repr__(self):
        return '%s.%s(labels=%s, k=%i)' % (
            self.__class__.__module__,
            self.__class__.__name__,
            self.y,
            self.k,
        )

    def __len__(self):
        return self.k


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

    indices: boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    Examples
    ----------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
    >>> y = np.array([1, 2, 1, 2])
    >>> labels = np.array([1, 1, 2, 2])
    >>> lol = cross_validation.LeaveOneLabelOut(labels)
    >>> len(lol)
    2
    >>> print lol
    sklearn.cross_validation.LeaveOneLabelOut(labels=[1 1 2 2])
    >>> for train_index, test_index in lol:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    ...    print X_train, X_test, y_train, y_test
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

    indices: boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    Examples
    ----------
    >>> from sklearn import cross_validation
    >>> X = np.array([[1, 2], [3, 4], [5, 6]])
    >>> y = np.array([1, 2, 1])
    >>> labels = np.array([1, 2, 3])
    >>> lpl = cross_validation.LeavePLabelOut(labels, p=2)
    >>> len(lpl)
    3
    >>> print lpl
    sklearn.cross_validation.LeavePLabelOut(labels=[1 2 3], p=2)
    >>> for train_index, test_index in lpl:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test = X[train_index], X[test_index]
    ...    y_train, y_test = y[train_index], y[test_index]
    ...    print X_train, X_test, y_train, y_test
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
        return (factorial(self.n_unique_labels) /
                factorial(self.n_unique_labels - self.p) /
                factorial(self.p))


class Bootstrap(object):
    """Random sampling with replacement cross-validation iterator

    Provides train/test indices to split data in train test sets
    while resampling the input n_bootstraps times: each time a new
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

    n_bootstraps : int (default is 3)
        Number of bootstrapping iterations

    n_train : int or float (default is 0.5)
        If int, number of samples to include in the training split
        (should be smaller than the total number of samples passed
        in the dataset).

        If float, should be between 0.0 and 1.0 and represent the
        proportion of the dataset to include in the train split.

    n_test : int or float or None (default is None)
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
    >>> print bs
    Bootstrap(9, n_bootstraps=3, n_train=5, n_test=4, random_state=0)
    >>> for train_index, test_index in bs:
    ...    print "TRAIN:", train_index, "TEST:", test_index
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

    def __init__(self, n, n_bootstraps=3, n_train=0.5, n_test=None,
                 random_state=None):
        self.n = n
        self.n_bootstraps = n_bootstraps

        if isinstance(n_train, float) and n_train >= 0.0 and n_train <= 1.0:
            self.n_train = ceil(n_train * n)
        elif isinstance(n_train, int):
            self.n_train = n_train
        else:
            raise ValueError("Invalid value for n_train: %r" % n_train)
        if self.n_train > n:
            raise ValueError("n_train=%d should not be larger than n=%d" %
                             (self.n_train, n))

        if isinstance(n_test, float) and n_test >= 0.0 and n_test <= 1.0:
            self.n_test = ceil(n_test * n)
        elif isinstance(n_test, int):
            self.n_test = n_test
        elif n_test is None:
            self.n_test = self.n - self.n_train
        else:
            raise ValueError("Invalid value for n_test: %r" % n_test)
        if self.n_test > n:
            raise ValueError("n_test=%d should not be larger than n=%d" %
                             (self.n_test, n))

        self.random_state = random_state

    def __iter__(self):
        rng = self.random_state = check_random_state(self.random_state)
        for i in range(self.n_bootstraps):
            # random partition
            permutation = rng.permutation(self.n)
            ind_train = permutation[:self.n_train]
            ind_test = permutation[self.n_train:self.n_train + self.n_test]

            # bootstrap in each split individually
            train = rng.randint(0, self.n_train, size=(self.n_train,))
            test = rng.randint(0, self.n_test, size=(self.n_test,))
            yield ind_train[train], ind_test[test]

    def __repr__(self):
        return ('%s(%d, n_bootstraps=%d, n_train=%d, n_test=%d, '
                'random_state=%d)' % (
                    self.__class__.__name__,
                    self.n,
                    self.n_bootstraps,
                    self.n_train,
                    self.n_test,
                    self.random_state,
                ))

    def __len__(self):
        return self.n_bootstraps


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

    n_iterations : int (default 10)
        Number of re-shuffling & splitting iterations.

    test_fraction : float (default 0.1)
        Should be between 0.0 and 1.0 and represent the proportion of
        the dataset to include in the test split.

    train_fraction : float or None (default is None)
        Should be between 0.0 and 1.0 and represent the proportion of
        the dataset to include in the train split. If None, the value is
        automatically set to the complement of the test fraction.

    indices : boolean, optional (default True)
        Return train/test split as arrays of indices, rather than a boolean
        mask array. Integer indices are required when dealing with sparse
        matrices, since those cannot be indexed by boolean masks.

    random_state : int or RandomState
        Pseudo-random number generator state used for random sampling.

    Examples
    --------
    >>> from sklearn import cross_validation
    >>> rs = cross_validation.ShuffleSplit(4, n_iterations=3,
    ...     test_fraction=.25, random_state=0)
    >>> len(rs)
    3
    >>> print rs
    ... # doctest: +ELLIPSIS
    ShuffleSplit(4, n_iterations=3, test_fraction=0.25, indices=True, ...)
    >>> for train_index, test_index in rs:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...
    TRAIN: [3 1 0] TEST: [2]
    TRAIN: [2 1 3] TEST: [0]
    TRAIN: [0 2 1] TEST: [3]

    >>> rs = cross_validation.ShuffleSplit(4, n_iterations=3,
    ...     train_fraction=0.5, test_fraction=.25, random_state=0)
    >>> for train_index, test_index in rs:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...
    TRAIN: [3 1] TEST: [2]
    TRAIN: [2 1] TEST: [0]
    TRAIN: [0 2] TEST: [3]

    See also
    --------
    Bootstrap: cross-validation using re-sampling with replacement.
    """

    def __init__(self, n, n_iterations=10, test_fraction=0.1,
                 train_fraction=None, indices=True, random_state=None):
        self.n = n
        self.n_iterations = n_iterations
        self.test_fraction = test_fraction
        self.train_fraction = train_fraction
        self.random_state = random_state
        self.indices = indices
        if test_fraction >= 1.0:
            raise ValueError(
                "test_fraction=%f should be smaller than 1.0" % test_fraction)
        if (train_fraction is not None
            and train_fraction + test_fraction > 1.0):
            raise ValueError(
                'The sum of train_fraction=%f and test_fraction=%f '
                'should be smaller or equal than 1.0' %
                (train_fraction, test_fraction))

    def __iter__(self):
        rng = self.random_state = check_random_state(self.random_state)
        n_test = ceil(self.test_fraction * self.n)
        if self.train_fraction is None:
            n_train = self.n - n_test
        else:
            n_train = floor(self.train_fraction * self.n)
        for i in range(self.n_iterations):
            # random partition
            permutation = rng.permutation(self.n)
            ind_test = permutation[:n_test]
            ind_train = permutation[n_test:n_test + n_train]

            if self.indices:
                yield ind_train, ind_test
            else:
                train_mask = np.zeros(self.n, dtype=np.bool)
                train_mask[ind_train] = True
                test_mask = np.zeros(self.n, dtype=np.bool)
                test_mask[ind_test] = True
                yield train_mask, test_mask

    def __repr__(self):
        return ('%s(%d, n_iterations=%d, test_fraction=%s, indices=%s, '
                'random_state=%s)' % (
                    self.__class__.__name__,
                    self.n,
                    self.n_iterations,
                    str(self.test_fraction),
                    self.indices,
                    self.random_state,
                ))

    def __len__(self):
        return self.n_iterations


##############################################################################

def _cross_val_score(estimator, X, y, score_func, train, test):
    """Inner loop for cross validation"""
    if y is None:
        estimator.fit(X[train])
        if score_func is None:
            return estimator.score(X[test])
        else:
            return score_func(X[test])
    else:
        estimator.fit(X[train], y[train])
        if score_func is None:
            return estimator.score(X[test], y[test])
        else:
            return score_func(y[test], estimator.predict(X[test]))


def cross_val_score(estimator, X, y=None, score_func=None, cv=None, n_jobs=1,
                    verbose=0):
    """Evaluate a score by cross-validation

    Parameters
    ----------
    estimator: estimator object implementing 'fit'
        The object to use to fit the data

    X: array-like of shape at least 2D
        The data to fit.

    y: array-like, optional
        The target variable to try to predict in the case of
        supervised learning.

    score_func: callable, optional
        callable, has priority over the score function in the estimator.
        In a non-supervised setting, where y is None, it takes the test
        data (X_test) as its only argument. In a supervised setting it takes
        the test target (y_true) and the test prediction (y_pred) as arguments.

    cv: cross-validation generator, optional
        A cross-validation generator. If None, a 3-fold cross
        validation is used or 3-fold stratified cross-validation
        when y is supplied and estimator is a classifier.

    n_jobs: integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    verbose: integer, optional
        The verbosity level
    """
    X, y = check_arrays(X, y, sparse_format='csr')
    cv = check_cv(cv, X, y, classifier=is_classifier(estimator))
    if score_func is None:
        assert hasattr(estimator, 'score'), ValueError(
                "If no score_func is specified, the estimator passed "
                "should have a 'score' method. The estimator %s "
                "does not." % estimator)
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(_cross_val_score)(clone(estimator), X, y, score_func,
                                          train, test)
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
        for label in np.unique(labels):
            this_mask = (labels == label)
            ind[this_mask] = random_state.permutation(ind[this_mask])
    return y[ind]


def check_cv(cv, X=None, y=None, classifier=False):
    """Input checker utility for building a CV in a user friendly way.

    Parameters
    ----------
    cv: an integer, a cv generator instance, or None
        The input specifying which cv generator to use. It can be an
        integer, in which case it is the number of folds in a KFold,
        None, in which case 3 fold is used, or another object, that
        will then be used as a cv generator.

    X: 2D ndarray
        the data the cross-val object will be applied on

    y: 1D ndarray
        the target variable for a supervised learning problem

    classifier: boolean optional
        whether the task is a classification task, in which case
        stratified KFold will be used.
    """
    is_sparse = sp.issparse(X)
    if cv is None:
        cv = 3
    if operator.isNumberType(cv):
        if classifier:
            cv = StratifiedKFold(y, cv, indices=is_sparse)
        else:
            if not is_sparse:
                n_samples = len(X)
            else:
                n_samples = X.shape[0]
            cv = KFold(n_samples, cv, indices=is_sparse)
    if is_sparse and not getattr(cv, "indices", True):
        raise ValueError("Sparse data require indices-based cross validation"
                         " generator, got: %r", cv)
    return cv


def permutation_test_score(estimator, X, y, score_func, cv=None,
                      n_permutations=100, n_jobs=1, labels=None,
                      random_state=0, verbose=0):
    """Evaluate the significance of a cross-validated score with permutations

    Parameters
    ----------
    estimator: estimator object implementing 'fit'
        The object to use to fit the data

    X: array-like of shape at least 2D
        The data to fit.

    y: array-like
        The target variable to try to predict in the case of
        supervised learning.

    score_func: callable
        Callable taking as arguments the test targets (y_test) and
        the predicted targets (y_pred) and returns a float. The score
        functions are expected to return a bigger value for a better result
        otherwise the returned value does not correspond to a p-value (see
        Returns below for further details).

    cv : integer or crossvalidation generator, optional
        If an integer is passed, it is the number of fold (default 3).
        Specific crossvalidation objects can be passed, see
        sklearn.cross_validation module for the list of possible objects

    n_jobs: integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.

    labels: array-like of shape [n_samples] (optional)
        Labels constrain the permutation among groups of samples with
        a same label.

    random_state: RandomState or an int seed (0 by default)
        A random number generator instance to define the state of the
        random permutations generator.

    verbose: integer, optional
        The verbosity level

    Returns
    -------
    score: float
        The true score without permuting targets.

    permutation_scores : array, shape = [n_permutations]
        The scores obtained for each permutations.

    pvalue: float
        The returned value equals p-value if `score_func` returns bigger
        numbers for better scores (e.g., zero_one). If `score_func` is rather a
        loss function (i.e. when lower is better such as with
        `mean_square_error`) then this is actually the complement of the
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
