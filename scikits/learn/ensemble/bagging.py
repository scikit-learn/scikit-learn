from .base import BaseEnsemble
import numpy as np


class Bagged(BaseEnsemble):

    def fit(self, X, Y, sample_weight=[],
                        sample_fraction=.5, baggs=10, **params):
        """
        X: list of instance vectors
        Y: target values/classes
        sample_fraction: fraction of X and Y randomly sampled
        baggs: number of sampling/training iterations
        """
        if not 0 < sample_fraction < 1:
            raise ValueError(
                "You must specify sample_fraction between 0 and 1 (exclusive)")
        if baggs < 2:
            raise ValueError("baggs must be greater than 1")

        if len(sample_weight) == 0:
            # initialize weights to 1/N
            sample_weight = np.ones(X.shape[0], dtype=np.float64)
        else:
            sample_weight = np.copy(sample_weight)
        # remove any previous ensemble
        self[:] = []
        for bagg in xrange(baggs):
            # weight a random subsample with 0
            random_sample_weight = sample_weight * \
                (np.random.random_sample(sample_weight.shape[0]) \
                    < sample_fraction)
            estimator = self.estimator(**self.params)
            estimator.fit(X, Y, random_sample_weight, **params)
            self.append(estimator)
        return self

    def predict(self, X):

        if len(self) == 0:
            return None
        prediction = np.zeros(X.shape[0], dtype=np.float64)
        for estimator in self:
            prediction += estimator.predict(X)
        prediction /= len(self)
        return prediction
    
    def predict_single(self, x):

        if len(self) == 0:
            return None
        prediction = 0.
        for estimator in self:
            prediction += estimator.predict_single(x)
        prediction /= len(self)
        return prediction
