
from ..base import BaseEstimator
"""
Base class for all ensemble classes
"""
class BaseEnsemble(BaseEstimator):

    def __init__(self, boosts, estimator, **params):

        self.estimator = estimator
        self.params = params
        self.estimators = []

    def __len__(self):

        return len(self.estimators)

    def __getitem__(self, index):

        return self.estimators[index]

    def __setitem__(self, index, thing):

        self.estimators[index] = thing

    def __delitem__(self, index):

        del self.estimators[index]
    
    def append(self, thing):

        return self.estimators.append(thing)
