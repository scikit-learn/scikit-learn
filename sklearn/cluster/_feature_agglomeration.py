"""
Feature agglomeration. Base classes and functions for performing feature
agglomeration.
"""
# Author: V. Michel, A. Gramfort
# License: BSD 3 clause

import numpy as np

from ..base import TransformerMixin
from ..utils import check_array

import warnings


###############################################################################
# Mixin class for feature agglomeration.

class AgglomerationTransform(TransformerMixin):
    """
    A class for feature agglomeration via the transform interface
    """

    pooling_func = np.mean

    def transform(self, X, pooling_func=None):
        """
        Transform a new matrix using the built clustering

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features] or [n_features]
            A M by N array of M observations in N dimensions or a length
            M array of M one-dimensional observations.

        pooling_func : callable, default=np.mean
            This combines the values of agglomerated features into a single
            value, and should accept an array of shape [M, N] and the keyword
            argument `axis=1`, and reduce it to an array of size [M].

        Returns
        -------
        Y : array, shape = [n_samples, n_clusters] or [n_clusters]
            The pooled values for each feature cluster.
        """
        if pooling_func is not None:
            warnings.warn("The pooling_func parameter is deprecated since 0.15 and will be "
                "removed in 0.18. Pass it to the constructor instead.", DeprecationWarning)
        else:
            pooling_func = self.pooling_func
        X = check_array(X)
        nX = []
        if len(self.labels_) != X.shape[1]:
            raise ValueError("X has a different number of features than "
                             "during fitting.")

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
        Xred : array-like, shape=[n_samples, n_clusters] or [n_clusters,]
            The values to be assigned to each cluster of samples

        Returns
        -------
        X : array, shape=[n_samples, n_features] or [n_features]
            A vector of size n_samples with the values of Xred assigned to
            each of the cluster of samples.
        """
        unil, inverse = np.unique(self.labels_, return_inverse=True)
        return Xred[..., inverse]
