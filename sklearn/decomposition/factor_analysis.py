"""Factor Analysis.
A latent variabel model that assumes diagonal, but in contrast to PCA not
isotropic, Gaussian latent variables. 

This implementation is based on Bishop's book chapter 12.2.4.
"""


# Author: Andreas Mueller <amueller@ais.uni-bonn.de>
# Licence: BSD3

import numpy as np
from scipy import linalg

from ..base import BaseEstimator, TransformerMixin
from ..utils import array2d, as_float_array, check_random_state


class FactorAnalysis(BaseEstimator, TransformerMixin):
    """Factor Analysis (FA)

    Linear dimensionality reduction using a latent variabel model. The
    Observations are assumed to be caused by a linear transformation of some
    lower dimensional `latent` factors.
    Factor analysis assumes these factors have a diagonal Gaussian
    distribution.

    If we would restrict the model further, by assuming that the Gaussian is
    even isotropic (all diagonal entries are the same) we would obtain
    :class:`PCA`.

    FactorAnalysis performs a maximum likelihood estimate of the so-called
    `loading` matrix, the transformation of the latent variables to the
    observed ones, using expectation-maximization.

    Parameters
    ----------
    n_components : int
        Dimensionality of latent space, the number of components
        of ``X`` that are obtained after ``transform``.

    copy : bool
        Whether to make a copy of X. If ``False``, the input X gets overwritten
        during fitting.

    tol : float
        Stopping tolerance for EM algorithm.

    random_state : int or RandomState instance or None (default)
        Pseudo Random Number generator seed control. If None, use the
        numpy.random singleton.


    References
    ----------
    .. Christopher M. Bishop: Pattern Recognition and Machine Learning,
        Chapter 12.2.4


    See also
    --------
    PCA: Principal component analysis, a very simliar, slightly simpler model
        that can be computed in closed form.
    ICA: Independent component analysis, a latent variable model with
        non-Gaussian latent variables.
    """
    def __init__(self, n_components, tol=.001, copy=True, random_state=None, max_iter=1000):
        self.n_components = n_components
        self.copy = copy
        self.tol = tol
        self.random_state = random_state
        self.max_iter = max_iter

    def fit(self, X, y=None):
        """Fit the FA model to X using EM.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        self
        """
        self.random_state = check_random_state(self.random_state)

        X = array2d(X)
        n_samples, n_features = X.shape
        X = as_float_array(X, copy=self.copy)
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_
        self.cov_ = np.cov(X.T)

        if self.n_components >= n_features:
            raise ValueError("n_components must be smaller than n_features.")

        # initialize model:
        loadings = self.random_state.uniform(size=(n_features,
            self.n_components))
        uniqueness = np.diag(self.cov_)
        latent = np.zeros((n_samples, self.n_components))

        for i in xrange(self.max_iter):
            # expectation step, find latent representation
            inv = linalg.inv(np.eye(self.n_components) + np.dot(loadings.T /
                uniqueness, loadings))
            latent = np.dot(X, np.dot(loadings, inv.T)
                    / uniqueness[:, np.newaxis])
            latent_cov = n_samples * inv + np.dot(latent.T, latent)
            latent_times_X = np.dot(latent.T, X)
            # maximization step, optimize loadings and cov
            loadings = linalg.solve(latent_cov.T, latent_times_X, sym_pos=True).T
            uniqueness = np.diag(self.cov_
                   - np.dot(loadings, latent_times_X / n_samples))

        self.loadings_ = loadings
        self.uniqueness_ = np.diag(uniqueness)
        self.components_ = np.dot(inv, loadings.T) / uniqueness
        return self

    def transform(self, X):
        """Apply dimensionality reduction to X using the model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
            The latent factors of X.
        """
        X_transformed = X - self.mean_
        X_transformed = np.dot(X_transformed, self.components_.T)
        return X_transformed

    def fit_transform(self, X, y=None):
        """Fit the FA model to X and then apply dimensionality reduction to X
        using the learned model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
            The latent factors of X.
        """
        return self.fit(X).transform(X)
