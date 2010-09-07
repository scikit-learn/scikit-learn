# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

from math import sqrt

import numpy as np
from scipy import linalg

from .base import BaseEstimator

class PCA(BaseEstimator):
    """
    Principal component analysis (PCA)

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    Attributes
    ----------
    k : int
        Number of components
        If k is not set all components are kept

    copy : bool
        If False, data passed to fit are overwritten

    components_ : array, [n_features, k]
        Components with maximum variance

    explained_variance_ : array, [k]
        Percentage of variance explained by each of the selected components.
        k is not set then all components are stored and the sum of
        explained variances is equal to 1.0

    Methods
    -------
    fit(X) : self
        Fit the model

    transform(X) : array
        Apply dimension reduction to k components


    Examples
    --------
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = PCA(k=2)
    >>> pca.fit(X)
    PCA(k=2, copy=True)
    >>> print pca.explained_variance_
    [ 0.99244289  0.00755711]

    See also
    --------

    """
    def __init__(self, k=None, copy=True):
        self.k = k
        self.copy = copy

    def fit(self, X, **params):
        self._set_params(**params)
        n_samples = X.shape[0]
        X = np.atleast_2d(X)
        if self.copy:
            X = X.copy()
        # Center data
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_
        U, S, V = linalg.svd(X, full_matrices=False)
        self.explained_variance_ = S**2
        self.explained_variance_ /= self.explained_variance_.sum()
        self.components_ = V
        if self.k is not None:
            self.components_ = self.components_[:,:self.k]
            self.explained_variance_ = self.explained_variance_[:self.k]

        return self

    def transform(self, X):
        Xr = X - self.mean_
        Xr = np.dot(Xr, self.components_)
        return Xr

