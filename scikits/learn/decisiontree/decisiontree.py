from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from . import libdecisiontree

class DecisionTree(BaseEstimator, ClassifierMixin): pass

class DecisionTreeRegressor(BaseEstimator, RegressorMixin): pass

class AlternatingDecisionTree(BaseEstimator, ClassifierMixin): pass

class RandomProjectionTree(BaseEstimator, ClassifierMixin): pass
