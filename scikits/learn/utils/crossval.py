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
    >>> import scikits.learn.utils.crossval
    >>> n_samples, n_features = 5, 10
    >>> X = np.random.randn(n_samples, n_features)
    >>> loo = crossval.LOO(n_samples)
    >>> for train_index, test_index in loo:
    ...    print "TRAIN:", train_index, "TEST:", test_index
    ...    Xtrain, Xtest, Ytrain, Ytest = split(train_index, test_index, X, y)
    ...    print Xtrain, Xtest, Ytrain, Ytest
    """
    for i in xrange(n):
        test_index  = np.zeros(n, dtype=np.bool)
        test_index[i] = True
        train_index = np.logical_not(test_index)
        yield train_index, test_index


def split(train_indexes, test_indexes, *args):
    """
    For each arg return a train and test subsets defined by indexes provided
    in train_indexes and test_indexes
    """
    ret = []
    for arg in args:
        arg_train = arg[train_indexes,:]
        arg_test  = arg[test_indexes,:]
        ret.append(arg_train)
        ret.append(arg_test)
    return ret


if __name__ == "__main__":
    print "Leave One Out crossvalidation"
    n_samples, n_features = 3, 4
    X = np.random.randn(n_samples, n_features)
    y = np.random.randn(n_samples)
    print X
    loo = leave_one_out(n_samples)
    for train_index, test_index in loo:
        print "TRAIN:", train_index, "TEST:", test_index
        Xtrain, Xtest, Ytrain, Ytest = split(train_index, test_index, X, y)
        print Xtrain, Xtest, Ytrain, Ytest
