"""Utilities for cross validation and performance evaluation"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>,
#         Gael Varoquaux    <gael.varoquaux@normalesup.org>
# License: BSD Style.

from math import ceil
import numpy as np

from .base import is_classifier, clone
from .utils.extmath import factorial, combinations
from .utils.fixes import unique
from .externals.joblib import Parallel, delayed


class LeaveOneOut(object):
    """Leave-One-Out cross validation iterator

    Provides train/test indices to split data in train test sets
    """

    def __init__(self, n, indices=False):
        """Leave-One-Out cross validation iterator

        Provides train/test indices to split data in train test sets

        Parameters
        ===========
        n: int
            Total number of elements
        indices: boolean, optional (default False)
            Return train/test split with integer indices or boolean mask.
            Integer indices are useful when dealing with sparse matrices
            that cannot be indexed by boolean masks.

        Examples
        ========
        >>> from scikits.learn import cross_val
        >>> X = np.array([[1, 2], [3, 4]])
        >>> y = np.array([1, 2])
        >>> loo = cross_val.LeaveOneOut(2)
        >>> len(loo)
        2
        >>> print loo
        scikits.learn.cross_val.LeaveOneOut(n=2)
        >>> for train_index, test_index in loo:
        ...    print "TRAIN:", train_index, "TEST:", test_index
        ...    X_train, X_test = X[train_index], X[test_index]
        ...    y_train, y_test = y[train_index], y[test_index]
        ...    print X_train, X_test, y_train, y_test
        TRAIN: [False  True] TEST: [ True False]
        [[3 4]] [[1 2]] [2] [1]
        TRAIN: [ True False] TEST: [False  True]
        [[1 2]] [[3 4]] [1] [2]
        """
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
        return '%s.%s(n=%i)' % (self.__class__.__module__,
                                self.__class__.__name__,
                                self.n,
                                )

    def __len__(self):
        return self.n


class LeavePOut(object):
    """Leave-P-Out cross validation iterator

    Provides train/test indices to split data in train test sets
    """

    def __init__(self, n, p, indices=False):
        """Leave-P-Out cross validation iterator

        Provides train/test indices to split data in train test sets

        Parameters
        ===========
        n: int
            Total number of elements
        p: int
            Size test sets
        indices: boolean, optional (default False)
            Return train/test split with integer indices or boolean mask.
            Integer indices are useful when dealing with sparse matrices
            that cannot be indexed by boolean masks.

        Examples
        ========
        >>> from scikits.learn import cross_val
        >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
        >>> y = np.array([1, 2, 3, 4])
        >>> lpo = cross_val.LeavePOut(4, 2)
        >>> len(lpo)
        6
        >>> print lpo
        scikits.learn.cross_val.LeavePOut(n=4, p=2)
        >>> for train_index, test_index in lpo:
        ...    print "TRAIN:", train_index, "TEST:", test_index
        ...    X_train, X_test = X[train_index], X[test_index]
        ...    y_train, y_test = y[train_index], y[test_index]
        TRAIN: [False False  True  True] TEST: [ True  True False False]
        TRAIN: [False  True False  True] TEST: [ True False  True False]
        TRAIN: [False  True  True False] TEST: [ True False False  True]
        TRAIN: [ True False False  True] TEST: [False  True  True False]
        TRAIN: [ True False  True False] TEST: [False  True False  True]
        TRAIN: [ True  True False False] TEST: [False False  True  True]
        """
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
        return factorial(self.n) / factorial(self.n - self.p) \
               / factorial(self.p)


class KFold(object):
    """K-Folds cross validation iterator

    Provides train/test indices to split data in train test sets
    """

    def __init__(self, n, k, indices=False):
        """K-Folds cross validation iterator

        Provides train/test indices to split data in train test sets

        Parameters
        ----------
        n: int
            Total number of elements
        k: int
            number of folds
        indices: boolean, optional (default False)
            Return train/test split with integer indices or boolean mask.
            Integer indices are useful when dealing with sparse matrices
            that cannot be indexed by boolean masks.

        Examples
        --------
        >>> from scikits.learn import cross_val
        >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
        >>> y = np.array([1, 2, 3, 4])
        >>> kf = cross_val.KFold(4, k=2)
        >>> len(kf)
        2
        >>> print kf
        scikits.learn.cross_val.KFold(n=4, k=2)
        >>> for train_index, test_index in kf:
        ...    print "TRAIN:", train_index, "TEST:", test_index
        ...    X_train, X_test = X[train_index], X[test_index]
        ...    y_train, y_test = y[train_index], y[test_index]
        TRAIN: [False False  True  True] TEST: [ True  True False False]
        TRAIN: [ True  True False False] TEST: [False False  True  True]

        Notes
        -----
        All the folds have size trunc(n/k), the last one has the complementary
        """
        assert k > 0, ('cannot have k below 1')
        assert k <= n, ('cannot have k=%d greater than the number '
                        'of samples: %d' % (k, n))
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

    Provides train/test indices to split data in train test sets

    This cross-validation object is a variation of KFold, which
    returns stratified folds. The folds are made by preserving
    the percentage of samples for each class.
    """

    def __init__(self, y, k, indices=False):
        """K-Folds cross validation iterator

        Provides train/test indices to split data in train test sets

        Parameters
        ----------
        y: array, [n_samples]
            Samples to split in K folds
        k: int
            number of folds
        indices: boolean, optional (default False)
            Return train/test split with integer indices or boolean mask.
            Integer indices are useful when dealing with sparse matrices
            that cannot be indexed by boolean masks.

        Examples
        --------
        >>> from scikits.learn import cross_val
        >>> X = np.array([[1, 2], [3, 4], [1, 2], [3, 4]])
        >>> y = np.array([0, 0, 1, 1])
        >>> skf = cross_val.StratifiedKFold(y, k=2)
        >>> len(skf)
        2
        >>> print skf
        scikits.learn.cross_val.StratifiedKFold(labels=[0 0 1 1], k=2)
        >>> for train_index, test_index in skf:
        ...    print "TRAIN:", train_index, "TEST:", test_index
        ...    X_train, X_test = X[train_index], X[test_index]
        ...    y_train, y_test = y[train_index], y[test_index]
        TRAIN: [False  True False  True] TEST: [ True False  True False]
        TRAIN: [ True False  True False] TEST: [False  True False  True]

        Notes
        -----
        All the folds have size trunc(n/k), the last one has the complementary
        """
        y = np.asanyarray(y)
        n = y.shape[0]
        assert k > 0, ValueError('cannot have k below 1')
        assert k <= n, ValueError('cannot have k=%d greater than the number '
                                  'of samples %d' % (k, n))
        _, y_sorted = unique(y, return_inverse=True)
        assert k <= np.min(np.bincount(y_sorted))
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


##############################################################################

class LeaveOneLabelOut(object):
    """Leave-One-Label_Out cross-validation iterator

    Provides train/test indices to split data in train test sets
    """

    def __init__(self, labels, indices=False):
        """Leave-One-Label_Out cross validation

        Provides train/test indices to split data in train test sets

        Parameters
        ----------
        labels : list
                List of labels
        indices: boolean, optional (default False)
            Return train/test split with integer indices or boolean mask.
            Integer indices are useful when dealing with sparse matrices
            that cannot be indexed by boolean masks.

        Examples
        ----------
        >>> from scikits.learn import cross_val
        >>> X = np.array([[1, 2], [3, 4], [5, 6], [7, 8]])
        >>> y = np.array([1, 2, 1, 2])
        >>> labels = np.array([1, 1, 2, 2])
        >>> lol = cross_val.LeaveOneLabelOut(labels)
        >>> len(lol)
        2
        >>> print lol
        scikits.learn.cross_val.LeaveOneLabelOut(labels=[1 1 2 2])
        >>> for train_index, test_index in lol:
        ...    print "TRAIN:", train_index, "TEST:", test_index
        ...    X_train, X_test = X[train_index], X[test_index]
        ...    y_train, y_test = y[train_index], y[test_index]
        ...    print X_train, X_test, y_train, y_test
        TRAIN: [False False  True  True] TEST: [ True  True False False]
        [[5 6]
         [7 8]] [[1 2]
         [3 4]] [1 2] [1 2]
        TRAIN: [ True  True False False] TEST: [False False  True  True]
        [[1 2]
         [3 4]] [[5 6]
         [7 8]] [1 2] [1 2]

        """
        self.labels = labels
        self.n_labels = unique(labels).size
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
        return self.n_labels


class LeavePLabelOut(object):
    """Leave-P-Label_Out cross-validation iterator

    Provides train/test indices to split data in train test sets
    """

    def __init__(self, labels, p, indices=False):
        """Leave-P-Label_Out cross validation

        Provides train/test indices to split data in train test sets

        Parameters
        ----------
        labels : list
            List of labels
        indices: boolean, optional (default False)
            Return train/test split with integer indices or boolean mask.
            Integer indices are useful when dealing with sparse matrices
            that cannot be indexed by boolean masks.

        Examples
        ----------
        >>> from scikits.learn import cross_val
        >>> X = np.array([[1, 2], [3, 4], [5, 6]])
        >>> y = np.array([1, 2, 1])
        >>> labels = np.array([1, 2, 3])
        >>> lpl = cross_val.LeavePLabelOut(labels, p=2)
        >>> len(lpl)
        3
        >>> print lpl
        scikits.learn.cross_val.LeavePLabelOut(labels=[1 2 3], p=2)
        >>> for train_index, test_index in lpl:
        ...    print "TRAIN:", train_index, "TEST:", test_index
        ...    X_train, X_test = X[train_index], X[test_index]
        ...    y_train, y_test = y[train_index], y[test_index]
        ...    print X_train, X_test, y_train, y_test
        TRAIN: [False False  True] TEST: [ True  True False]
        [[5 6]] [[1 2]
         [3 4]] [1] [1 2]
        TRAIN: [False  True False] TEST: [ True False  True]
        [[3 4]] [[1 2]
         [5 6]] [2] [1 1]
        TRAIN: [ True False False] TEST: [False  True  True]
        [[1 2]] [[3 4]
         [5 6]] [1] [2 1]

        """
        self.labels = labels
        self.unique_labels = unique(self.labels)
        self.n_labels = self.unique_labels.size
        self.p = p
        self.indices = indices

    def __iter__(self):
        # We make a copy here to avoid side-effects during iteration
        labels = np.array(self.labels, copy=True)
        unique_labels = unique(labels)
        n_labels = unique_labels.size
        comb = combinations(range(n_labels), self.p)

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
        return factorial(self.n_labels) / factorial(self.n_labels - self.p) \
               / factorial(self.p)


def _cross_val_score(estimator, X, y, score_func, train, test, iid):
    """Inner loop for cross validation"""
    if score_func is None:
        score_func = lambda self, *args: self.score(*args)
    if y is None:
        score = score_func(estimator.fit(X[train]), X[test])
    else:
        score = score_func(estimator.fit(X[train], y[train]), X[test], y[test])
    if iid:
        if y is not None:
            score *= len(y[test])
        else:
            score *= len(X[test])
    return score


def cross_val_score(estimator, X, y=None, score_func=None, cv=None, iid=False,
                n_jobs=1, verbose=0):
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
        callable taking as arguments the fitted estimator, the
        test data (X_test) and the test target (y_test) if y is
        not None.
    cv: cross-validation generator, optional
        A cross-validation generator. If None, a 3-fold cross
        validation is used or 3-fold stratified cross-validation
        when y is supplied and estimator is a classifier.
    iid: boolean, optional
        If True, the data is assumed to be identically distributed across
        the folds, and the loss minimized is the total loss per sample,
        and not the mean loss across the folds.
    n_jobs: integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.
    verbose: integer, optional
        The verbosity level
    """
    n_samples = len(X)
    if cv is None:
        if y is not None and is_classifier(estimator):
            cv = StratifiedKFold(y, k=3)
        else:
            cv = KFold(n_samples, k=3)
    if score_func is None:
        assert hasattr(estimator, 'score'), ValueError(
                "If no score_func is specified, the estimator passed "
                "should have a 'score' method. The estimator %s "
                "does not." % estimator)
    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(_cross_val_score)(clone(estimator), X, y, score_func,
                                                        train, test, iid)
                for train, test in cv)
    return np.array(scores)


def _permutation_test_score(estimator, X, y, cv, score_func):
    """Auxilary function for permutation_test_score
    """
    y_test = list()
    y_pred = list()
    for train, test in cv:
        y_test.append(y[test])
        y_pred.append(estimator.fit(X[train], y[train]).predict(X[test]))
    return score_func(np.ravel(y_test), np.ravel(y_pred))


def _shuffle(y, labels, rng):
    """Return a shuffled copy of y eventually shuffle among same labels.
    """
    if labels is None:
        ind = rng.permutation(y.size)
    else:
        ind = np.arange(labels.size)
        for label in np.unique(labels):
            this_mask = (labels == label)
            ind[this_mask] = rng.permutation(ind[this_mask])
    return y[ind]


def permutation_test_score(estimator, X, y, score_func, cv=None,
                      n_permutations=100, n_jobs=1, labels=None,
                      rng=0, verbose=0):
    """Evaluate the significance of a cross-validated score with permutations

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
        callable taking as arguments the test targets (y_test) and
        the predicted targets (y_pred). Returns a float.
    cv: cross-validation generator, optional
        A cross-validation generator. If None, a 3-fold cross
        validation is used or 3-fold stratified cross-validation
        when the estimator is a classifier.
    n_jobs: integer, optional
        The number of CPUs to use to do the computation. -1 means
        'all CPUs'.
    labels: array-like of shape [n_samples] (optional)
        Labels constrain the permutation among groups of samples with
        a same label.
    rng: RandomState or an int seed (0 by default)
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
        The p-value.

    Notes
    -----
    In corresponds to Test 1 in :
    Ojala and Garriga. Permutation Tests for Studying Classifier Performance.
    The Journal of Machine Learning Research (2010) vol. 11
    """
    n_samples = len(X)
    if cv is None:
        if is_classifier(estimator):
            cv = StratifiedKFold(y, k=3)
        else:
            cv = KFold(n_samples, k=3)

    if rng is None:
        rng = np.random.RandomState()
    elif isinstance(rng, int):
        rng = np.random.RandomState(rng)

    # We clone the estimator to make sure that all the folds are
    # independent, and that it is pickle-able.
    score = _permutation_test_score(clone(estimator), X, y, cv, score_func)
    permutation_scores = Parallel(n_jobs=n_jobs, verbose=verbose)(
                delayed(_permutation_test_score)(clone(estimator), X,
                                            _shuffle(y, labels, rng),
                                            cv, score_func)
                for _ in range(n_permutations))
    permutation_scores = np.array(permutation_scores)
    pvalue = (np.sum(permutation_scores >= score) + 1.0) / (n_permutations + 1)
    return score, permutation_scores, pvalue


permutation_test_score.__test__ = False  # to avoid a pb with nosetests
