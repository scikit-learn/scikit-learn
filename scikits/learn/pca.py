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

    Methods
    -------
    fit(X) : self
        Fit the model

    transform(X) : array
        Apply dimension reduction to k components


    Examples
    --------
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> pca = PCA(k=1)
    >>> pca.fit(X)
    PCA(k=1, copy=True)
    >>> print pca.components_
    [[ 0.83849224]
     [-0.54491354]]

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
        scores = np.dot(U, S)
        S /= sqrt(n_samples-1)
        latent = S**2
        
        self.components_ = V
        if self.k is not None:
            self.components_ = self.components_[:,:self.k]

        return self

    def transform(self, X):
        Xr = X - self.mean_
        Xr = np.dot(Xr, self.components_)
        return Xr

