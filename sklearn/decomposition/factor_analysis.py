"""Factor Analysis.
A latent linear variable model, similar to ProbabilisticPCA.

This implementation is based on David Barber's Book,
Bayesian Reasoning and Machine Learning,
http://www.cs.ucl.ac.uk/staff/d.barber/brml,
Algorithm 21.1
"""

# Author: Christian Osendorfer <osendorf@gmail.com>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
# Licence: BSD3

from math import sqrt
import numpy as np
from scipy import linalg


from ..base import BaseEstimator, TransformerMixin
from ..externals.six.moves import xrange
from ..utils import array2d, check_arrays
from ..utils.extmath import fast_logdet


class FactorAnalysis(BaseEstimator, TransformerMixin):
    """Factor Analysis (FA)

    A simple linear generative model with Gaussian latent variables.

    The observations are assumed to be caused by a linear transformation of
    lower dimensional latent factors and added Gaussian noise.
    Without loss of generality the factors are distributed according to a
    Gaussian with zero mean and unit covariance. The noise is also zero mean
    and has an arbitrary diagonal covariance matrix.

    If we would restrict the model further, by assuming that the Gaussian
    noise is even isotropic (all diagonal entries are the same) we would obtain
    :class:`PPCA`.

    FactorAnalysis performs a maximum likelihood estimate of the so-called
    `loading` matrix, the transformation of the latent variables to the
    observed ones, using expectation-maximization (EM).

    Parameters
    ----------
    n_components : int | None
        Dimensionality of latent space, the number of components
        of ``X`` that are obtained after ``transform``.
        If None, n_components is set to the number of features.

    tol : float
        Stopping tolerance for EM algorithm.

    copy : bool
        Whether to make a copy of X. If ``False``, the input X gets overwritten
        during fitting.

    max_iter : int
        Maximum number of iterations.

    verbose : int | bool
        Print verbose output.

    noise_variance_init : None | array, shape=(n_features,)
        The initial guess of the noise variance for each feature.
        If None, it defaults to np.ones(n_features)

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Components with maximum variance.

    `loglike_` : list, [n_iterations]
        The log likelihood at each iteration.

    `noise_variance_` : array, shape=(n_features,)
        The estimated noise variance for each feature.

    References
    ----------
    .. David Barber, Bayesian Reasoning and Machine Learning,
        Algorithm 21.1

    .. Christopher M. Bishop: Pattern Recognition and Machine Learning,
        Chapter 12.2.4

    See also
    --------
    PCA: Principal component analysis, a similar non-probabilistic
        model model that can be computed in closed form.
    ProbabilisticPCA: probabilistic PCA.
    FastICA: Independent component analysis, a latent variable model with
        non-Gaussian latent variables.
    """
    def __init__(self, n_components=None, tol=1e-2, copy=True, max_iter=1000,
                 verbose=0, noise_variance_init=None):
        self.n_components = n_components
        self.copy = copy
        self.tol = tol
        self.max_iter = max_iter
        self.verbose = verbose
        self.noise_variance_init = noise_variance_init

    def fit(self, X, y=None):
        """Fit the FactorAnalysis model to X using EM

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        Returns
        -------
        self
        """
        X = array2d(check_arrays(X, copy=self.copy, sparse_format='dense',
                    dtype=np.float)[0])

        n_samples, n_features = X.shape
        n_components = self.n_components
        if n_components is None:
            n_components = n_features
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_

        # some constant terms
        nsqrt = sqrt(n_samples)
        llconst = n_features * np.log(2 * np.pi) + n_components
        var = np.var(X, axis=0)

        if self.noise_variance_init is None:
            psi = np.ones(n_features, dtype=X.dtype)
        else:
            if len(self.noise_variance_init) != n_features:
                raise ValueError("noise_variance_init dimension does not "
                                 "with number of features : %d != %d" %
                                 (len(self.noise_variance_init), n_features))
            psi = np.array(self.noise_variance_init)

        loglike = []
        old_ll = -np.inf
        SMALL = 1e-12
        for i in xrange(self.max_iter):
            # SMALL helps numerics
            sqrt_psi = np.sqrt(psi) + SMALL
            Xtilde = X / (sqrt_psi * nsqrt)
            _, s, V = linalg.svd(Xtilde, full_matrices=False)
            V = V[:n_components]
            s **= 2
            # Use 'maximum' here to avoid sqrt problems.
            W = np.sqrt(np.maximum(s[:n_components] - 1, 0))[:, np.newaxis] * V
            W *= sqrt_psi

            # loglikelihood
            ll = llconst + np.sum(np.log(s[:n_components]))
            ll += np.sum(s[n_components:]) + np.sum(np.log(psi))
            ll *= -n_samples / 2.
            loglike.append(ll)
            if (ll - old_ll) < self.tol:
                break
            old_ll = ll

            psi = np.maximum(var - np.sum(W ** 2, axis=0), SMALL)
        else:
            if self.verbose:
                print("Did not converge")

        self.components_ = W
        self.noise_variance_ = psi
        self.loglike_ = loglike
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
        X = array2d(X)
        Ih = np.eye(len(self.components_))

        X_transformed = X - self.mean_

        Wpsi = self.components_ / self.noise_variance_
        cov_z = linalg.inv(Ih + np.dot(Wpsi, self.components_.T))
        tmp = np.dot(X_transformed, Wpsi.T)
        X_transformed = np.dot(tmp, cov_z)

        return X_transformed

    def get_covariance(self):
        """Compute data covariance with the FactorAnalysis model.

        ``cov = components_.T * components_ + diag(noise_variance)``

        Returns
        -------
        cov : array, shape=(n_features, n_features)
            Estimated covariance of data.
        """
        cov = (np.dot(self.components_.T, self.components_)
               + np.diag(self.noise_variance_))
        return cov

    def score(self, X, y=None):
        """Compute score of X under FactorAnalysis model.

        Parameters
        ----------
        X: array of shape(n_samples, n_features)
            The data to test

        Returns
        -------
        ll: array of shape (n_samples),
            log-likelihood of each row of X under the current model
        """
        Xr = X - self.mean_
        cov = self.get_covariance()
        n_features = X.shape[1]
        log_like = np.zeros(X.shape[0])
        self.precision_ = linalg.inv(cov)
        log_like = -.5 * (Xr * (np.dot(Xr, self.precision_))).sum(axis=1)
        log_like -= .5 * (fast_logdet(cov) + n_features * np.log(2 * np.pi))
        return log_like
