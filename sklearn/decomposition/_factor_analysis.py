"""Factor Analysis.

A latent linear variable model.

FactorAnalysis is similar to probabilistic PCA implemented by PCA.score
While PCA assumes Gaussian noise with the same variance for each
feature, the FactorAnalysis model assumes different variances for
each of them.

This implementation is based on David Barber's Book,
Bayesian Reasoning and Machine Learning,
http://www.cs.ucl.ac.uk/staff/d.barber/brml,
Algorithm 21.1
"""

# Author: Christian Osendorfer <osendorf@gmail.com>
#         Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Denis A. Engemann <denis-alexander.engemann@inria.fr>

# License: BSD3

import warnings
from math import sqrt, log
import numpy as np
from scipy import linalg


from ..base import BaseEstimator, TransformerMixin
from ..utils import check_array, check_random_state
from ..utils.extmath import fast_logdet, randomized_svd, squared_norm
from ..utils.validation import check_is_fitted, _deprecate_positional_args
from ..exceptions import ConvergenceWarning


class FactorAnalysis(TransformerMixin, BaseEstimator):
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
    observed ones, using SVD based approach.

    Read more in the :ref:`User Guide <FA>`.

    .. versionadded:: 0.13

    Parameters
    ----------
    n_components : int | None
        Dimensionality of latent space, the number of components
        of ``X`` that are obtained after ``transform``.
        If None, n_components is set to the number of features.

    tol : float
        Stopping tolerance for log-likelihood increase.

    copy : bool
        Whether to make a copy of X. If ``False``, the input X gets overwritten
        during fitting.

    max_iter : int
        Maximum number of iterations.

    noise_variance_init : None | array, shape=(n_features,)
        The initial guess of the noise variance for each feature.
        If None, it defaults to np.ones(n_features)

    svd_method : {'lapack', 'randomized'}
        Which SVD method to use. If 'lapack' use standard SVD from
        scipy.linalg, if 'randomized' use fast ``randomized_svd`` function.
        Defaults to 'randomized'. For most applications 'randomized' will
        be sufficiently precise while providing significant speed gains.
        Accuracy can also be improved by setting higher values for
        `iterated_power`. If this is not sufficient, for maximum precision
        you should choose 'lapack'.

    iterated_power : int, optional
        Number of iterations for the power method. 3 by default. Only used
        if ``svd_method`` equals 'randomized'

    random_state : int, RandomState instance, default=0
        Only used when ``svd_method`` equals 'randomized'. Pass an int for
        reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    components_ : array, [n_components, n_features]
        Components with maximum variance.

    loglike_ : list, [n_iterations]
        The log likelihood at each iteration.

    noise_variance_ : array, shape=(n_features,)
        The estimated noise variance for each feature.

    n_iter_ : int
        Number of iterations run.

    mean_ : array, shape (n_features,)
        Per-feature empirical mean, estimated from the training set.

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.decomposition import FactorAnalysis
    >>> X, _ = load_digits(return_X_y=True)
    >>> transformer = FactorAnalysis(n_components=7, random_state=0)
    >>> X_transformed = transformer.fit_transform(X)
    >>> X_transformed.shape
    (1797, 7)

    References
    ----------
    .. David Barber, Bayesian Reasoning and Machine Learning,
        Algorithm 21.1

    .. Christopher M. Bishop: Pattern Recognition and Machine Learning,
        Chapter 12.2.4

    See also
    --------
    PCA: Principal component analysis is also a latent linear variable model
        which however assumes equal noise variance for each feature.
        This extra assumption makes probabilistic PCA faster as it can be
        computed in closed form.
    FastICA: Independent component analysis, a latent variable model with
        non-Gaussian latent variables.
    """
    @_deprecate_positional_args
    def __init__(self, n_components=None, *, tol=1e-2, copy=True,
                 max_iter=1000,
                 noise_variance_init=None, svd_method='randomized',
                 iterated_power=3, random_state=0):
        self.n_components = n_components
        self.copy = copy
        self.tol = tol
        self.max_iter = max_iter
        if svd_method not in ['lapack', 'randomized']:
            raise ValueError('SVD method %s is not supported. Please consider'
                             ' the documentation' % svd_method)
        self.svd_method = svd_method

        self.noise_variance_init = noise_variance_init
        self.iterated_power = iterated_power
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the FactorAnalysis model to X using SVD based approach

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data.

        y : Ignored

        Returns
        -------
        self
        """
        X = self._validate_data(X, copy=self.copy, dtype=np.float64)

        n_samples, n_features = X.shape
        n_components = self.n_components
        if n_components is None:
            n_components = n_features
        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_

        # some constant terms
        nsqrt = sqrt(n_samples)
        llconst = n_features * log(2. * np.pi) + n_components
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

        # we'll modify svd outputs to return unexplained variance
        # to allow for unified computation of loglikelihood
        if self.svd_method == 'lapack':
            def my_svd(X):
                _, s, Vt = linalg.svd(X, full_matrices=False)
                return (s[:n_components], Vt[:n_components],
                        squared_norm(s[n_components:]))
        elif self.svd_method == 'randomized':
            random_state = check_random_state(self.random_state)

            def my_svd(X):
                _, s, Vt = randomized_svd(X, n_components,
                                          random_state=random_state,
                                          n_iter=self.iterated_power)
                return s, Vt, squared_norm(X) - squared_norm(s)
        else:
            raise ValueError('SVD method %s is not supported. Please consider'
                             ' the documentation' % self.svd_method)

        for i in range(self.max_iter):
            # SMALL helps numerics
            sqrt_psi = np.sqrt(psi) + SMALL
            s, Vt, unexp_var = my_svd(X / (sqrt_psi * nsqrt))
            s **= 2
            # Use 'maximum' here to avoid sqrt problems.
            W = np.sqrt(np.maximum(s - 1., 0.))[:, np.newaxis] * Vt
            del Vt
            W *= sqrt_psi

            # loglikelihood
            ll = llconst + np.sum(np.log(s))
            ll += unexp_var + np.sum(np.log(psi))
            ll *= -n_samples / 2.
            loglike.append(ll)
            if (ll - old_ll) < self.tol:
                break
            old_ll = ll

            psi = np.maximum(var - np.sum(W ** 2, axis=0), SMALL)
        else:
            warnings.warn('FactorAnalysis did not converge.' +
                          ' You might want' +
                          ' to increase the number of iterations.',
                          ConvergenceWarning)

        self.components_ = W
        self.noise_variance_ = psi
        self.loglike_ = loglike
        self.n_iter_ = i + 1
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
        check_is_fitted(self)

        X = check_array(X)
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
        cov : array, shape (n_features, n_features)
            Estimated covariance of data.
        """
        check_is_fitted(self)

        cov = np.dot(self.components_.T, self.components_)
        cov.flat[::len(cov) + 1] += self.noise_variance_  # modify diag inplace
        return cov

    def get_precision(self):
        """Compute data precision matrix with the FactorAnalysis model.

        Returns
        -------
        precision : array, shape (n_features, n_features)
            Estimated precision of data.
        """
        check_is_fitted(self)

        n_features = self.components_.shape[1]

        # handle corner cases first
        if self.n_components == 0:
            return np.diag(1. / self.noise_variance_)
        if self.n_components == n_features:
            return linalg.inv(self.get_covariance())

        # Get precision using matrix inversion lemma
        components_ = self.components_
        precision = np.dot(components_ / self.noise_variance_, components_.T)
        precision.flat[::len(precision) + 1] += 1.
        precision = np.dot(components_.T,
                           np.dot(linalg.inv(precision), components_))
        precision /= self.noise_variance_[:, np.newaxis]
        precision /= -self.noise_variance_[np.newaxis, :]
        precision.flat[::len(precision) + 1] += 1. / self.noise_variance_
        return precision

    def score_samples(self, X):
        """Compute the log-likelihood of each sample

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The data

        Returns
        -------
        ll : array, shape (n_samples,)
            Log-likelihood of each sample under the current model
        """
        check_is_fitted(self)

        Xr = X - self.mean_
        precision = self.get_precision()
        n_features = X.shape[1]
        log_like = -.5 * (Xr * (np.dot(Xr, precision))).sum(axis=1)
        log_like -= .5 * (n_features * log(2. * np.pi)
                          - fast_logdet(precision))
        return log_like

    def score(self, X, y=None):
        """Compute the average log-likelihood of the samples

        Parameters
        ----------
        X : array, shape (n_samples, n_features)
            The data

        y : Ignored

        Returns
        -------
        ll : float
            Average log-likelihood of the samples under the current model
        """
        return np.mean(self.score_samples(X))
