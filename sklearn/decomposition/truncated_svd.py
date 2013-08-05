# -*- coding: utf-8 -*-

"""Truncated SVD for sparse matrices, aka latent semantic analysis (LSA).
"""

# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
# License: 3-clause BSD.

import numpy as np

try:
    from scipy.sparse.linalg import svds
except ImportError:
    from ..utils.arpack import svds

from ..base import BaseEstimator, TransformerMixin
from ..utils import (array2d, as_float_array, atleast2d_or_csr,
                     check_random_state)
from ..utils.extmath import randomized_svd, safe_sparse_dot, svd_flip

__all__ = ["TruncatedSVD"]


class TruncatedSVD(BaseEstimator, TransformerMixin):
    """Dimensionality reduction using truncated SVD (aka LSA).

    This transformer performs linear dimensionality reduction by means of
    truncated singular value decomposition (SVD). It is very similar to PCA,
    but operates on sample vectors directly, instead of on a covariance matrix.
    This means it can work with scipy.sparse matrices efficiently.

    In particular, truncated SVD works on term count/tf–idf matrices as
    returned by the vectorizers in sklearn.feature_extraction.text. In that
    context, it is known as latent semantic analysis (LSA).

    Parameters
    ----------
    n_components : int, default = 2
        Desired dimensionality of output data.
        Must be strictly less than the number of features.
        The default value is useful for visualisation. For LSA, a value of
        100 is recommended.
    algorithm : string, default = "randomized"
        SVD solver to use. Either "arpack" for the ARPACK wrapper in SciPy
        (scipy.sparse.linalg.svds), or "randomized" for the randomized
        algorithm due to Halko (2009).
    n_iterations : int, optional
        Number of iterations for randomized SVD solver. Not used by ARPACK.
    random_state : int or RandomState, optional
        (Seed for) pseudo-random number generator. If not given, the
        numpy.random singleton is used.
    tol : float, optional
        Tolerance for ARPACK. 0 means machine precision. Ignored by randomized
        SVD solver.

    Attributes
    ----------
    `components_` : array, shape (n_components, n_features)

    See also
    --------
    PCA
    RandomizedPCA

    References
    ----------
    Finding structure with randomness: Stochastic algorithms for constructing
    approximate matrix decompositions
    Halko, et al., 2009 (arXiv:909) http://arxiv.org/pdf/0909.4061

    Notes
    -----
    SVD suffers from a problem called "sign indeterminancy", which means the
    sign of the ``components_`` and the output from transform depend on the
    algorithm and random state. To work around this, fit instances of this
    class to data once, then keep the instance around to do transformations.

    """
    def __init__(self, n_components=2, algorithm="randomized",
                 n_iterations=5, random_state=None, tol=0.):
        self.algorithm = algorithm
        self.n_components = n_components
        self.n_iterations = n_iterations
        self.random_state = random_state
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
        X = as_float_array(X, copy=False)
        random_state = check_random_state(self.random_state)

        if self.algorithm == "arpack":
            U, Sigma, VT = svds(X, k=self.n_components, tol=self.tol)
            # svds doesn't abide by scipy.linalg.svd/randomized_svd
            # conventions, so reverse its outputs.
            Sigma = Sigma[::-1]
            U, VT = svd_flip(U[:, ::-1], VT[::-1])

        elif self.algorithm == "randomized":
            k = self.n_components
            n_features = X.shape[1]
            if k >= n_features:
                raise ValueError("n_components must be < n_features;"
                                 " got %d >= %d" % (k, n_features))
            U, Sigma, VT = randomized_svd(X, self.n_components,
                                          n_iter=self.n_iterations,
                                          random_state=random_state)
        else:
            raise ValueError("unknown algorithm %r" % self.algorithm)

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

    def inverse_transform(self, X):
        """Transform X back to its original space.

        Returns an array X_original whose transform would be X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            New data.

        Returns
        -------
        X_original : array, shape (n_samples, n_features)
            Note that this is always a dense array.
        """
        X = array2d(X)
        return np.dot(X, self.components_)
