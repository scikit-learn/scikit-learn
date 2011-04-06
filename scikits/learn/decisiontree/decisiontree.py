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

    def fit(self, X, Y, sample_weight = np.empty(0)):

        self.root = libdecisiontree.fit(X, Y, sample_weight, self.minleafsize, self.nbins)
        return self

    def predict(self, X):

        return libdecisiontree.predict(X, self.root)
        
"""
class DecisionTreeRegressor(BaseEstimator, RegressorMixin): pass

class AlternatingDecisionTree(BaseEstimator, ClassifierMixin): pass

class RandomProjectionTree(BaseEstimator, ClassifierMixin): pass
"""
