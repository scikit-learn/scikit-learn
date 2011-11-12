from ..base import clone
from .base import BaseEnsemble
import numpy as np


class Bagged(BaseEnsemble):
    def __init__(self, estimator, sample_fraction=.5, baggs=10, **params):
        """
        sample_fraction: fraction of X and Y randomly sampled
        baggs: number of sampling/training iterations
        """
        if not 0 < sample_fraction < 1:
            raise ValueError(
                "You must specify sample_fraction between 0 and 1 (exclusive)")
        if baggs < 2:
            raise ValueError("baggs must be greater than 1")

        self.sample_fraction = sample_fraction
        self.baggs = baggs

        super(Bagged, self).__init__(estimator, 

    def fit(self, X, Y, sample_weight=[],
                        sample_fraction=.5, baggs=10, **fit_params):
        """
        X: list of instance vectors
        Y: target values/classes
        """
        if len(sample_weight) == 0:
            # initialize weights to 1/N
            sample_weight = np.ones(X.shape[0], dtype=np.float64)
        else:
            sample_weight = np.copy(sample_weight)

        for i in xrange(self.baggs):
            # weight a random subsample with 0
            random_sample_weight = sample_weight * \
                (np.random.random_sample(sample_weight.shape[0]) \
                    < self.sample_fraction)
            estimator = self.estimators[i]
            estimator.fit(X, Y, random_sample_weight, **fit_params)
        return self

    def predict(self, X):
        X = np.atleast_2d(X)
        if len(self) == 0:
            return None
        prediction = np.sum([est.predict(X) for est in self.estimators],
                            axis=0)
        prediction /= len(self)
        return prediction
