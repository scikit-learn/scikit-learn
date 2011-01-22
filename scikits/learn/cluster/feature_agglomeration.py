"""
Feature agglomeration. Classes and functions for performing feature
agglomeration.
"""
# Author: V. Michel, A. Gramfort
# License: BSD 3 clause

import numpy as np
from ..base import BaseEstimator
from .hierarchical import Ward

###############################################################################
# Mixin class for feature agglomeration.

class AgglomerationTransform(BaseEstimator):
    """
    A class for feature agglomeration via the transform interface
    """

    def transform(self, X, pooling_func=np.mean):
        """
        Transform a new matrix using the built clustering

        Parameters
        ---------
        X : array-like, shape = [n_samples, n_features]
            A M by N array of M observations in N dimensions or a length
            M array of M one-dimensional observations.

        pooling_func : a function that takes an array of shape = [M, N] and
                       return an array of value of size M.
                       Defaut is np.mean
        """
        nX = []
        for l in np.unique(self.labels_):
            nX.append(pooling_func(X[:, self.labels_ == l], axis=1))
        return np.array(nX).T

    def inverse_transform(self, Xred):
        """
        Inverse the transformation.
        Return a vector of size nb_features with the values of Xred assigned
        to each group of features

        Parameters
        ----------
        Xred : array of size k
            The values to be assigned to each cluster of samples

        Return
        ------
        X : array of size nb_samples
            A vector of size nb_samples with the values of Xred assigned to
            each of the cluster of samples.
        """
        if np.size((Xred.shape)) == 1:
            X = np.zeros([self.labels_.shape[0]])
        else:
            X = np.zeros([Xred.shape[0], self.labels_.shape[0]])
        unil = np.unique(self.labels_)
        for i in range(len(unil)):
            if np.size((Xred.shape)) == 1:
                X[self.labels_ == unil[i]] = Xred[i]
            else:
                ncol = np.sum(self.labels_ == unil[i])
                X[:, self.labels_ == unil[i]] = np.tile(np.atleast_2d(Xred
                                                        [:, i]).T, ncol)
        return X

###############################################################################
# Cluster based features agglomeration objects

class WardAgglomeration(AgglomerationTransform):
    """Feature agglomeration base on Ward hierarchical clustering

    XXX
    """

    def __init__(self, k=2, adjacency_matrix=None, copy=True):
        self.k = k
        self.adjacency_matrix = adjacency_matrix
        self.copy = copy
        self._ward = Ward(k)

    def fit(self, X, y=None, **params):
        self._set_params(**params)
        self._ward.fit(X.T, k=self.k, adjacency_matrix=self.adjacency_matrix,
                        copy=self.copy, **params)
        self.labels_ = self._ward.labels_
        return self
