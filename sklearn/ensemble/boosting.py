from ..base import clone
from .base import BaseEnsemble

import math
import numpy as np


class AdaBoost(BaseEnsemble):
    def __init__(self, estimator, boosts=1, beta=.5, **params):
        if boosts < 1:
            raise ValueError(
                "You must specify a number of boosts greater than 0")
        if beta <= 0:
            raise ValueError("Beta must be positive and non-zero")

        self.beta = beta
        self.boosts = boosts

        super(AdaBoost, self).__init__(self, estimator, boosts + 1)
        self.estimator_weights_ = np.ones(boosts + 1)

    def fit(self, X, Y, sample_weight=[], **fit_params):
        """
        X: list of instance vectors
        Y: target values/classes
        sample_weight: sample (X) weights

        Notes: currently only binary classification is supported
        I am making the assumption that one class label is
        positive and the other is negative
        """

        if len(sample_weight) == 0:
            # initialize weights to 1/N
            sample_weight = np.ones(X.shape[0], dtype=np.float64)\
                / X.shape[0]
        else:
            sample_weight = np.copy(sample_weight)

        for i in xrange(self.boosts + 1):
            estimator = self.estimators[i]
            estimator.fit(X, Y, sample_weight, **fit_params)
            # TODO request that classifiers return classification
            # of training sets when fitting
            # which would make the following line unnecessary
            T = estimator.predict(X)
            # instances incorrectly classified
            incorrect = ((T * Y) < 0).astype(np.int32)
            # error fraction
            err = np.sum(sample_weight * incorrect) / np.sum(sample_weight)
            # sanity check
            if err >= 0.5:
                break
            # boost weight
            alpha = self.beta * math.log((1 - err) / err)
            self.alpha_[i] = alpha
            if i < self.boosts:
                correct = incorrect ^ 1
                sample_weight *= np.exp(alpha * (incorrect - correct))

        return self

    def predict(self, X):
        X = np.atleast_2d(X)
        prediction = np.zeros(X.shape[0], dtype=np.float64)
        norm = 0.
        for i in xrange(len(self.estimators)):
            alpha = self.alpha_[i]
            prediction += alpha * self.estimators[i].predict(X)
        norm = self.alpha_.sum()
        if norm > 0:
            prediction /= norm
        return prediction
    

"""
YET TO BE IMPLEMENTED

class GradientBoost(BaseEnsemble): pass

class StochasticGradientBoost(BaseEnsemble): pass
"""
