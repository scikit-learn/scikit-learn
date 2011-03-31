import numpy as np
from base import BaseBoost, clone

class AdaBoost(BaseBoost):

    def fit(self, X, Y, W = None):
        """
        X: list of instance vectors
        Y: target values/classes
        W: weights
        """
        if W is None:
            W = np.ones(X.shape[0], dtype = np.float64)
        for boost in xrange(self.boosts):
            if boost > 0:
                self.append(clone(self[-1]))
            classification = self[-1].fit(X,Y,W)

    def predict(self, X):
        
        prediction = np.zeros(len(X))
        for weight, estimator in self:
            prediction += weight * estimator
        prediction /= len(self)
        return prediction
