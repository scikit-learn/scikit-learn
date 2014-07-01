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

from .model_selection.partition import LeaveOneOut, LeavePOut, KFold, check_cv
from .model_selection.partition import StratifiedKFold, LeaveOneLabelOut
from .model_selection.partition import ShuffleSplit, StratifiedShuffleSplit
from .model_selection.partition import LeavePLabelOut, train_test_split
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
           'PredefinedSplit',
           'check_cv',
           'cross_val_score',
           'cross_val_predict',
           'permutation_test_score',
           'train_test_split']


# TODO: move Boostrap somewhere else and import it from here
# TODO: issue a DeprecationWarning when this module is imported

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
                 random_state=None):
        # See, e.g., http://youtu.be/BzHz0J9a6k0?t=9m38s for a motivation
        # behind this deprecation
        warnings.warn("Bootstrap will no longer be supported as a " +
                      "cross-validation method as of version 0.15 and " +
                      "will be removed in 0.17", DeprecationWarning)
        self.n = n
        self.n_iter = n_iter
        if isinstance(train_size, numbers.Integral):
            self.train_size = train_size
        elif (isinstance(train_size, numbers.Real) and train_size >= 0.0
                and train_size <= 1.0):
            self.train_size = int(ceil(train_size * n))
        else:
            raise ValueError("Invalid value for train_size: %r" %
                             train_size)
        if self.train_size > n:
            raise ValueError("train_size=%d should not be larger than n=%d" %
                             (self.train_size, n))

        if isinstance(test_size, numbers.Integral):
            self.test_size = test_size
        elif isinstance(test_size, numbers.Real) and 0.0 <= test_size <= 1.0:
            self.test_size = int(ceil(test_size * n))
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

    Parameters
    ----------
    estimator : estimator object implementing 'fit' and 'predict'
        The object to use to fit the data.

    X : array-like
        The data to fit. Can be, for example a list, or an array at least 2d.

    y : array-like, optional, default: None
        The target variable to try to predict in the case of
        supervised learning.

    cv : cross-validation generator or int, optional, default: None
        A cross-validation generator to use. If int, determines
        the number of folds in StratifiedKFold if y is binary
        or multiclass and estimator is a classifier, or the number
        of folds in KFold otherwise. If None, it is equivalent to cv=3.
        This generator must include all elements in the test set exactly once.
        Otherwise, a ValueError is raised.

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
    """
    X, y = indexable(X, y)

    cv = _check_cv(cv, X, y, classifier=is_classifier(estimator))
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    parallel = Parallel(n_jobs=n_jobs, verbose=verbose,
                        pre_dispatch=pre_dispatch)
    preds_blocks = parallel(delayed(_fit_and_predict)(clone(estimator), X, y,
                                                      train, test, verbose,
                                                      fit_params)
                            for train, test in cv)
    p = np.concatenate([p for p, _ in preds_blocks])
    locs = np.concatenate([loc for _, loc in preds_blocks])
    if not _check_is_partition(locs, X.shape[0]):
        raise ValueError('cross_val_predict only works for partitions')
    preds = p.copy()
    preds[locs] = p
    return preds


def _fit_and_predict(estimator, X, y, train, test, verbose, fit_params):
    """Fit estimator and predict values for a given dataset split.

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
