import numpy as np
from ..base import BaseEnsemble
import math

class AdaBoost(BaseEnsemble):

    def fit(self, X, Y, sample_weight = [], boosts = 1, **params):
        """
        X: list of instance vectors
        Y: target values/classes
        sample_weight: sample (X) weights

        Notes: currently only binary classification is supported
        I am making the assumption that one class label is positive and the other is negative
        """
        if len(sample_weight) == 0:
            sample_weight = np.ones(X.shape[0], dtype = np.float64)
        else:
            sample_weight = np.copy(sample_weight)
        # remove any previous ensemble
        self[:] = []
        for boost in xrange(boosts+1):
            estimator = self.estimator(**self.params)
            estimator.fit(X,Y,sample_weight,**params)
            # TODO request that classifiers return classification of training sets when fitting
            # which would make the following line unnecessary 
            T = estimator.predict(X)
            incorrect = (T*Y)<0
            err = np.sum(sample_weight * incorrect) / np.sum(sample_weight)
            alpha = math.log((1 - err) / err)
            sample_weight *= np.exp(alpha * incorrect)
            self.append((alpha, estimator))
        return self

    def predict(self, X):
        
        prediction = np.zeros(X.shape[0], dtype = np.float64)
        norm = 0.
        for alpha, estimator in self:
            prediction += alpha * estimator.predict(X)
            norm += alpha
        if norm > 0:
            prediction /= norm
        return prediction
