
from ..base import BaseEstimator
"""
Base class for all boosting classes
"""
class BaseBoost(BaseEstimator):

    def __init__(self, boosts, estimator, **params):

        self.estimator = estimator
        self.params = params
        self.estimators = []
        self.boosts = boosts

    def __len__(self):

        return len(self.estimators)

    def __getitem__(self, index):

        return self.estimators[index]

    def __setitem__(self, index, thing):

        return self.estimators[index] = thing

    def __delitem__(self, index):

        return del self.estimators[index]
    
    def append(self, thing):

        return self.estimators.append(thing)
