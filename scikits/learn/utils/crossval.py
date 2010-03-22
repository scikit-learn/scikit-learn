# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import exceptions

import numpy as np

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


def leave_one_label_out(label):
    """
    Leave-One-Label_Out cross validation:
    Provides train/test indexes to split data in train test sets

    Parameters
    ----------
    label : list
            List of labels

    Examples
    ----------
    >>> import numpy as np
    >>> from scikits.learn.utils import crossval
    >>> X = [[1, 2], [3, 4], [5, 6], [7, 8]]
    >>> y = [1, 2, 1, 2]
    >>> label = [1,1,2,2]
    >>> lol = crossval.leave_one_label_out(label)
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
    for i in np.unique(label):
        test_index  = np.zeros(len(label), dtype=np.bool)
        test_index[label==i] = True
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
