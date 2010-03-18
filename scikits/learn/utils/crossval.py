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
    # commented doctest, see issue #34
    # >>> import numpy as np
    # >>> from scikits.learn.utils import crossval
    # >>> n_samples, n_features = 5, 10
    # >>> X = np.random.randn(n_samples, n_features)
    # >>> loo = crossval.leave_one_out(n_samples)
    # >>> for train_index, test_index in loo:
    # ...    print "TRAIN:", train_index, "TEST:", test_index
    # ...    Xtrain, Xtest, Ytrain, Ytest = split(train_index, test_index, X, y)
    # ...    print Xtrain, Xtest, Ytrain, Ytest
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
    n_samples, n_features = 4, 2
    X = np.random.randn(n_samples, n_features)
    y = np.random.randn(n_samples)
    print X
    loo = leave_one_out(n_samples)
    for train_index, test_index in loo:
        print "TRAIN:", train_index, "TEST:", test_index
        Xtrain, Xtest, Ytrain, Ytest = split(train_index, test_index, X, y)
        print Xtrain, Xtest, Ytrain, Ytest

    print "K-Fold crossvalidation"
    k = 2
    kf = k_fold(n_samples, k)
    for train_index, test_index in kf:
        print "TRAIN:", train_index, "TEST:", test_index
