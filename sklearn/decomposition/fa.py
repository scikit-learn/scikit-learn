"""Factor Analysis.
A latent linear variable model, similar to PPCA.

This implementation is based on David Barber's Book,
Bayesian Reasoning and Machine Learning,
http://www.cs.ucl.ac.uk/staff/d.barber/brml,
Algorithm 21.1
"""


# Author: Christian Osendorfer <osendorf@gmail.com>
# Licence: BSD3

import numpy as np
from scipy import linalg

from ..base import BaseEstimator, TransformerMixin
from ..utils import as_float_array, check_random_state


class FA(BaseEstimator, TransformerMixin):
    """Factor Analysis (FA)

    A simple linear generative model with Guassian latent variables.

    The observations are assumed to be caused by a linear transformation of
    lower dimensional latent factors and added Gaussian noise.
    Wlog the factors are distributed according to a Gaussian with
    zero mean and unit covariance. The noise is also zero mean and
    has an arbitrary diagonal covariance matrix.

    If we would restrict the model further, by assuming that the Gaussian noise is
    even isotropic (all diagonal entries are the same) we would obtain
    :class:`PPCA`.

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
    .. David Barber, Bayesian Reasoning and Machine Learning,
        Algorithm 21.1


    See also
    --------
    PCA: Principal component analysis, a simliar non-probabilistic
        model model that can be computed in closed form.
    PPCA: probabilistic PCA.
    ICA: Independent component analysis, a latent variable model with
        non-Gaussian latent variables.
    """
    def __init__(self, n_components, tol=1e-2, copy=True, random_state=None, max_iter=1000, verbose=0):
        self.n_components = n_components
        self.copy = copy
        self.tol = tol
        self.random_state = random_state
        self.max_iter = max_iter
        self.verbose = verbose

    def fit(self, X, psi=None):
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
        n, d = X.shape

        np.atleast_2d(as_float_array(X, copy=self.copy))
        n_samples, n_features = X.shape
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_

        # some constant terms
        nsqrt = np.sqrt(n)
        llconst = d*np.log(2*np.pi) + self.n_components
        var = np.var(X, 0)

        if psi is None:
            psi = np.ones(d)

        loglike = []
        old_ll = -np.inf
        SMALL = 1e-8
        for i in xrange(self.max_iter):
            # SMALL helps numerics
            sqrt_psi = np.sqrt(psi) + SMALL
            Xtilde = X/(sqrt_psi * nsqrt)
            u, s, v = linalg.svd(Xtilde, full_matrices=False)
            v = v[:self.n_components]
            s *= s
            # Use 'maximum' here to avoid sqrt problems.
            W = np.sqrt(np.maximum(s[:self.n_components] - 1, 0))[:, np.newaxis]*v
            W *= sqrt_psi

            # loglikelihood
            ll = llconst + np.sum(np.log(s[:self.n_components]))
            ll += np.sum(s[self.n_components:]) + np.sum(np.log(psi))
            ll *= -n/2.
            loglike.append(ll)
            if ll - old_ll < self.tol:
                break
            old_ll = ll

            psi = var - np.sum(W**2, axis=0)
        self.W = W
        self.psi = psi
        self.loglike = loglike
        return self

    def transform(self, X):
        """Apply dimensionality reduction to X using the model.
        Compute the expected mean of the latent variables.
        See Barber, 21.2.33 (or Bishop, 12.66).

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
            The latent variables of X.
        """
        Ih = np.eye(self.n_components)

        X_transformed = X - self.mean_

        Wpsi = self.W/self.psi
        cov_z = linalg.inv(Ih + np.dot(Wpsi, self.W.T))
        tmp = np.dot(X_transformed, Wpsi.T)
        X_transformed = np.dot(tmp, cov_z)

        return X_transformed
