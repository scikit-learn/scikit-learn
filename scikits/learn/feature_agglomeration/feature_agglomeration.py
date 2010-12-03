"""
Feature agglomeration. Classes and functions for performing feature
agglomeration.
"""
# Author: V. Michel
# License: BSD 3 clause

import numpy as np
from scikits.learn.base import BaseEstimator
from scikits.learn.cluster import Ward

######################################################################
# General class for feature agglomeration.


class FeatureAgglomeration(BaseEstimator):
    """
    Class for feature agglomeration
    """

    def __init__(self, clustering, *args, **kwargs):
        """
        Initialize the feature agglomeration

        Parameters
        ---------
        clustering: a cluster method, that should have a .fit method, and a
                    label_ attributes.
        """
        self.clustering = clustering(*args, **kwargs)

    def fit(self, X, *args, **kwargs):
        """
        Learn the clustering on the data.

        Parameters
        ---------
        X : array-like, shape = [n_samples, n_features]
            A M by N array of M observations in N dimensions or a length
            M array of M one-dimensional observations.

        """
        cluster = self.clustering.fit(X.T, *args, **kwargs)
        self.label_ = cluster.label_
        return self

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
        for l in np.unique(self.label_):
            nX.append(pooling_func(X[:, self.label_ == l], 1))
        return np.array(nX).T

    def inverse_transform(self, Xred):
        """
        Inverse the transformation.
        Return a vector of size nb_features with the values of Xred assigned
        to each group of features

        Parameters
        ----------
        Xred : array of size k
        The values to be assigned to each cluster of features

        Return
        ------
        X : array of size nb_features
        A vector of size nb_features with the values of Xred assigned to each
        of the cluster of features.
        """
        if np.size((Xred.shape)) == 1:
            X = np.zeros([self.label_.shape[0]])
        else:
            X = np.zeros([Xred.shape[0], self.label_.shape[0]])
        unil = np.unique(self.label_)
        for i in range(len(unil)):
            if np.size((Xred.shape)) == 1:
                X[self.label_ == unil[i]] = Xred[i]
            else:
                ncol = np.sum(self.label_ == unil[i])
                X[:, self.label_ == unil[i]] = np.tile(np.atleast_2d(Xred
                                                        [:, i]).T, ncol)
        return X


######################################################################
# Classes for specific feature agglomeration.


class WardAgglomeration(FeatureAgglomeration):
    """
    Class for ward feature agglomeration
    """
    def __init__(self, *args, **kwargs):
        """
        Initialize the ward feature agglomeration
        """
        self.clustering = Ward(*args, **kwargs)
