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
from math import ceil
import numbers

from .utils import check_random_state

from .model_selection.partition import LeaveOneOut, LeavePOut, KFold, \
        StratifiedKFold, LeaveOneLabelOut, LeavePLabelOut, \
        ShuffleSplit, StratifiedShuffleSplit, train_test_split, check_cv
from .model_selection.validate import cross_val_score, permutation_test_score

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


#TODO: move Boostrap somewhere else and import it from here
#TODO: issue a DeprecationWarning when this module is imported

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
        # See, e.g., http://youtu.be/BzHz0J9a6k0?t=9m38s for a motivation
        # behind this deprecation
        warnings.warn("Bootstrap will no longer be supported as a " +
                      "cross-validation method as of version 0.15 and " +
                      "will be removed in 0.17", DeprecationWarning)
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
        if self.test_size > n - self.train_size:
            raise ValueError(("test_size + train_size=%d, should not be " +
                              "larger than n=%d") %
                             (self.test_size + self.train_size, n))

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


##############################################################################


