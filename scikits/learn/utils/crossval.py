# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np

try:
    from itertools import combinations
except: # Using Python < 2.6
    def combinations(seq, r=None):
        """Generator returning combinations of items from sequence <seq>
        taken <r> at a time. Order is not significant. If <r> is not given,
        the entire sequence is returned.
        """
        if r == None:
            r = len(seq)
        if r <= 0:
            yield []
        else:
            for i in xrange(len(seq)):
                for cc in combinations(seq[i+1:], r-1):
                    yield [seq[i]]+cc

def leave_one_out(n):
    """
    Leave-One-Out cross validation:
    Provides train/test indexes to split data in train test sets

    Parameters
    ===========
    n: int
        Total number of elements

    Examples
    ========
    >>> import numpy as np
    >>> from scikits.learn.utils import crossval
    >>> X = [[1, 2], [3, 4]]
    >>> y = [1, 2]
    >>> loo = crossval.leave_one_out(2)
    >>> for train_index, test_index in loo:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test, y_train, y_test = crossval.split(train_index, test_index, X, y)
    ...    print X_train, X_test, y_train, y_test
    TRAIN: [False  True] TEST: [ True False]
    [[3 4]] [[1 2]] [2] [1]
    TRAIN: [ True False] TEST: [False  True]
    [[1 2]] [[3 4]] [1] [2]
    """
    for i in xrange(n):
        test_index  = np.zeros(n, dtype=np.bool)
        test_index[i] = True
        train_index = np.logical_not(test_index)
        yield train_index, test_index

def leave_p_out(n, p):
    """
    Leave-P-Out cross validation:
    Provides train/test indexes to split data in train test sets

    Parameters
    ===========
    n: int
        Total number of elements
    p: int
        Size test sets

    Examples
    ========
    >>> import numpy as np
    >>> from scikits.learn.utils import crossval
    >>> X = [[1, 2], [3, 4], [5, 6], [7, 8]]
    >>> y = [1, 2, 3, 4]
    >>> lpo = crossval.leave_p_out(4, 2)
    >>> for train_index, test_index in lpo:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test, y_train, y_test = crossval.split(train_index, test_index, X, y)
    TRAIN: [False False  True  True] TEST: [ True  True False False]
    TRAIN: [False  True False  True] TEST: [ True False  True False]
    TRAIN: [False  True  True False] TEST: [ True False False  True]
    TRAIN: [ True False False  True] TEST: [False  True  True False]
    TRAIN: [ True False  True False] TEST: [False  True False  True]
    TRAIN: [ True  True False False] TEST: [False False  True  True]
    """
    comb = combinations(range(n), p)
    for idx in comb:
        test_index = np.zeros(n, dtype=np.bool)
        test_index[np.array(idx)] = True
        train_index = np.logical_not(test_index)
        yield train_index, test_index


def k_fold(n, k):
    """
    K-Folds cross validation:
    Provides train/test indexes to split data in train test sets

    Parameters
    ===========
    n: int
        Total number of elements
    k: int
        number of folds

    Examples
    ========
    >>> import numpy as np
    >>> from scikits.learn.utils import crossval
    >>> X = [[1, 2], [3, 4], [1, 2], [3, 4]]
    >>> y = [1, 2, 3, 4]
    >>> kf = crossval.k_fold(4, 2)
    >>> for train_index, test_index in kf:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test, y_train, y_test = crossval.split(train_index, test_index, X, y)
    TRAIN: [False False  True  True] TEST: [ True  True False False]
    TRAIN: [ True  True False False] TEST: [False False  True  True]

    Note
    ====
    All the folds have size trunc(n/k), the last one has the complementary
    """
    assert k>0, ValueError('cannot have k below 1')
    assert k<n, ValueError('cannot have k=%d greater than %d'% (k, n))
    j = np.ceil(n/k)

    for i in xrange(k):
        test_index  = np.zeros(n, dtype=np.bool)
        if i<k-1:
            test_index[i*j:(i+1)*j] = True
        else:
            test_index[i*j:] = True
        train_index = np.logical_not(test_index)
        yield train_index, test_index


def leave_one_label_out(labels):
    """
    Leave-One-Label_Out cross validation:
    Provides train/test indexes to split data in train test sets

    Parameters
    ----------
    labels : list
            List of labels

    Examples
    ----------
    >>> import numpy as np
    >>> from scikits.learn.utils import crossval
    >>> X = [[1, 2], [3, 4], [5, 6], [7, 8]]
    >>> y = [1, 2, 1, 2]
    >>> labels = [1, 1, 2, 2]
    >>> lol = crossval.leave_one_label_out(labels)
    >>> for train_index, test_index in lol:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    X_train, X_test, y_train, y_test = crossval.split(train_index, \
           test_index, X, y)
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
    for i in np.unique(labels):
        test_index  = np.zeros(len(labels), dtype=np.bool)
        test_index[labels==i] = True
        train_index = np.logical_not(test_index)
        yield train_index, test_index


def split(train_indexes, test_indexes, *args):
    """
    For each arg return a train and test subsets defined by indexes provided
    in train_indexes and test_indexes
    """
    ret = []
    for arg in args:
        arg = np.asanyarray(arg)
        arg_train = arg[train_indexes]
        arg_test  = arg[test_indexes]
        ret.append(arg_train)
        ret.append(arg_test)
    return ret

