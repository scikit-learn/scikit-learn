import numpy as np
#import cuv_python as cp
from sklearn.base import BaseEstimator

class FourierSampler(BaseEstimator):
    def fit_transform(self, X):
        return self.fit(X).transform(X)


class RBFSampler(FourierSampler):
    """ Following "Random Features for Large-Scale Kernel Machines" by A, Rahimi and
        Benjamin Recht """
    def __init__(self, gamma=1., D=1.):
        self.gamma = gamma
        self.D = D

    def fit(self, X, y=None):
        """ Samples random projection according to feature lenght of X"""
        n_features = X.shape[1]
        self.omega_ = np.sqrt(self.gamma) * np.random.normal(size=(n_features, self.D))
        self.b = np.random.uniform(0, 2*np.pi, size=self.D)
        return self

    def transform(self, X):
        projection = np.dot(X, self.omega_)
        return np.sqrt(2)/np.sqrt(self.D)*np.cos(projection + self.b)


class SkewedChi2Sampler(FourierSampler):
    """Following "Random Fourier Approximations for Skewed Multiplicative Histogram Kernels"
    by Fuxin Li, Catalin Ionescu and Cristian Sminchisescu
    Parameter "c" is how "skewed" the group is. Fuxin uses cross validation.
    """

    def __init__(self, c, D):
        self.c = c
        self.D = D

    def fit(self, X):
        """Samples random projection according to feature lenght of X"""
        n_features = X.shape[1]
        uniform = np.random.uniform(size=(n_features, self.D))
        # transform by inverse CDF of sech
        self.omega_ = 1./np.pi * np.log(np.tan(np.pi/2. * uniform))
        self.b = np.random.uniform(0, 2*np.pi, size=self.D)
        return self

    def transform(self, X):
        # see paper?
        projection = np.dot(np.log(X + self.c), self.omega_)
        return np.sqrt(2)/np.sqrt(self.D)*np.cos(projection + self.b)


class AdditiveChi2Approx(object):
    def fit_transform(self, X):
        return self.fit(X).transform(X)

class IntersectionApprox(object):
    def fit_transform(self, X):
        return self.fit(X).transform(X)
