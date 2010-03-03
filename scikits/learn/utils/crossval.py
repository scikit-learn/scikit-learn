# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id: cd.py 473 2010-03-03 16:27:38Z twigster $

import numpy as np
import exceptions

class LOO:
    """
    Leave-One-Out cross validation:
    Provides train/test indexes to split data in train test sets

    Examples:
    import scikits.learn.utils.crossval
    import numpy as np
    n_samples, n_features = 5, 10
    X = np.random.randn(n_samples, n_features)
    print X
    loo = crossval.LOO(n_samples)
    print loo[1]
    for train_index, test_index in loo:
        print "TRAIN:", train_index, "TEST:", test_index
        Xtrain, Xtest, Ytrain, Ytest = split(train_index, test_index, X, y)
        print Xtrain, Xtest, Ytrain, Ytest
    """

    def __init__(self,n):
        """
        n : is the size of the dataset to split
        """
        self.n_folds = n
        self.iter = 0

    def __getitem__(self,item):
        test_index  = np.zeros(self.n_folds,dtype=np.bool)
        test_index[item] = True
        train_index = np.logical_not(test_index)
        return train_index, test_index

    def next(self):
        if self.iter < self.n_folds:
            self.iter += 1
            return self.__getitem__(self.iter-1)
        raise StopIteration

    def __iter__(self):
        return self

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
    loo = LOO(n_samples)
    print loo[1]
    for train_index, test_index in loo:
        print "TRAIN:", train_index, "TEST:", test_index
        Xtrain, Xtest, Ytrain, Ytest = split(train_index, test_index, X, y)
        print Xtrain, Xtest, Ytrain, Ytest
