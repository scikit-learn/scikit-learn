from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from . import libdecisiontree
import numpy as np

class DecisionTree(BaseEstimator, ClassifierMixin):

    def __init__(self, minleafsize = None, maxdepth = -1, maxnodes = -1, nbins = 20, sepcriterion = "gini"):

        self.minleafsize = minleafsize
        self.maxdepth = maxdepth
        self.maxnodes = maxnodes
        self.nbins = nbins
        self.sepcriterion = sepcriterion
        self.root = None
        self.nfeatures = 0

    def fit(self, X, y, sample_weight = [], **params):
        
        self._set_params(**params)

        X = np.atleast_2d(np.asanyarray(X, dtype = np.float64, order = 'C'))
        y = np.asanyarray(y, dtype = np.float64, order = 'C')
        sample_weight = np.asanyarray(sample_weight, dtype = np.float64, order = 'C')

        if X.shape[0] != y.shape[0]:
            raise ValueError("X and y have incompatible shapes.\n" +
                             "X has %s samples, but y has %s." % \
                             (X.shape[0], y.shape[0]))

        if len(sample_weight) == 0:
            sample_weight = np.ones(X.shape[0], dtype=np.float64)
        elif sample_weight.shape[0] != X.shape[0]:
            raise ValueError("sample_weight and X have incompatible shapes.\n" +
                             "sample_weight has %s samples while X has %s" % \
                             (sample_weight.shape[0], X.shape[0]))
       
        labels = np.unique(y)
        if len(labels) != 2:
            raise ValueError("Only binary classificiation is supported.\n" +
                             "y must contain two unique values.")
        sig_score = labels[1]
        bkg_score = labels[0]
        print sig_score
        print bkg_score

        self.root = libdecisiontree.fit(X, y, sample_weight, self.minleafsize, self.nbins, sig_score, bkg_score)
        self.nfeatures = X.shape[1]
        return self

    def predict(self, X):

        if self.root is None:
            return None

        X = np.atleast_2d(np.asanyarray(X, dtype = np.float64, order = 'C'))
        
        if X.shape[1] != self.nfeatures:
            raise ValueError("X.shape[1] should be equal to the number of "
                             "features at training time!")

        return libdecisiontree.predict(X, self.root)
        
"""
class DecisionTreeRegressor(BaseEstimator, RegressorMixin): pass

class AlternatingDecisionTree(BaseEstimator, ClassifierMixin): pass

class RandomProjectionTree(BaseEstimator, ClassifierMixin): pass
"""
