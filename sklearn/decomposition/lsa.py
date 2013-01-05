"""Latent semantic analysis.
"""

# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: 3-clause BSD.

import numpy as np
from scipy.sparse.linalg import svds

from ..base import BaseEstimator, TransformerMixin
from ..utils import as_float_array, atleast2d_or_csr
from ..utils.extmath import safe_sparse_dot

__all__ = ["LSA"]


class LSA(BaseEstimator, TransformerMixin):
    """Dimensionality reduction using latent semantic analysis (LSA, aka LSI).

    LSA performs linear dimensionality reduction by means of truncated
    singular value decomposition (SVD). It is very similar to PCA, but operates
    on sample vectors directly, instead of on a covariance matrix. This means
    it can work with scipy.sparse matrices efficiently, in particular term
    count/tfidf matrices as returned by the vectorizers in
    sklearn.feature_extraction.text.

    This implementation uses the scipy.sparse.linalg.svds implementation of
    truncated SVD.

    Parameters
    ----------
    n_components : int, optional
        Desired dimensionality of output data.
        Must be strictly less than the number of features.
    tol : float, optional
        Tolerance for underlying eigensolver. 0 means machine precision.

    Attributes
    ----------
    components_ : array, shape (n_components, n_features)

    See also
    --------
    PCA
    RandomizedPCA
    """
    def __init__(self, n_components=100, tol=0.):
        self.n_components = n_components
        self.tol = tol

    def fit(self, X, y=None):
        """Fit LSI model on training data X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        self : object
            Returns the transformer object.
        """
        self._fit(X)
        return self

    def fit_transform(self, X, y=None):
        """Fit LSI model to X and perform dimensionality reduction on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Reduced version of X. This will always be a dense array.
        """
        U, Sigma, VT = self._fit(X)
        Sigma = np.diag(Sigma)

        # or (X * VT.T).T, whichever takes fewer operations...
        return np.dot(U, Sigma.T)

    def _fit(self, X):
        X = as_float_array(X)
        U, Sigma, VT = svds(X, k=self.n_components, tol=self.tol)
        self.components_ = VT
        return U, Sigma, VT

    def transform(self, X):
        """Perform dimensionality reduction on X.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            New data.

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Reduced version of X. This will always be a dense array.
        """
        X = atleast2d_or_csr(X)
        return safe_sparse_dot(X, self.components_.T)
