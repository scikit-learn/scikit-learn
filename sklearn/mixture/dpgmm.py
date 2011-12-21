"""Bayesian Gaussian Mixture Models and 
Dirichlet Process Gaussian Mixture Models"""

# Author: Alexandre Passos (alexandre.tp@gmail.com)
#
# Based on mixture.py by:
#         Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#

import numpy as np
from scipy.special import digamma as _digamma, gammaln as _gammaln
from scipy import linalg

from ..utils import check_random_state
from ..utils.extmath import norm
from .. import cluster
from ..metrics import euclidean_distances
from . gmm import GMM


def sqnorm(v):
    return norm(v) ** 2


def digamma(x):
    return _digamma(x + np.finfo(np.float32).eps)


def gammaln(x):
    return _gammaln(x + np.finfo(np.float32).eps)


def log_normalize(v, axis=0):
    """Normalized probabilities from unnormalized log-probabilites"""
    v = np.rollaxis(v, axis)
    v = v.copy()
    v -= v.max(axis=0)
    out = np.log(np.sum(np.exp(v), axis=0))
    v = np.exp(v - out)
    v += np.finfo(np.float32).eps
    v /= np.sum(v, axis=0)
    return np.swapaxes(v, 0, axis)


def detlog_wishart(a, b, detB, n_features):
    """Expected value of the log of the determinant of a Wishart

    The expected value of the logarithm of the determinant of a
    wishart-distributed random variable with the specified parameters."""
    l = np.sum(digamma(0.5 * (a - np.arange(-1, n_features - 1))))
    l += n_features * np.log(2)
    return l + detB


def wishart_logz(v, s, dets, n_features):
    "The logarithm of the normalization constant for the wishart distribution"
    z = 0.
    z += 0.5 * v * n_features * np.log(2)
    z += (0.25 * (n_features * (n_features - 1))
          * np.log(np.pi))
    z += 0.5 * v * np.log(dets)
    z += np.sum(gammaln(0.5 * (v - np.arange(n_features) + 1)))
    return z


##############################################################################
# Variational bound on the log likelihood of each class

def _bound_state_loglik_spherical(X, initial_bound, bound_prec, precs, means):
    n_components, n_features = means.shape
    n_samples = X.shape[0]
    bound = np.empty((n_samples, n_components))
    bound[:] = bound_prec + initial_bound
    for k in xrange(n_components):
        bound[:, k] -= 0.5 * precs[k] * (((X - means[k]) ** 2).sum(axis=-1)
                                         + n_features)
    return bound


def _bound_state_loglik_diag(X, initial_bound, bound_prec, precs, means):
    n_components, n_features = means.shape
    n_samples = X.shape[0]
    bound = np.empty((n_samples, n_components))
    bound[:] = bound_prec + initial_bound
    for k in xrange(n_components):
        d = X - means[k]
        d **= 2
        bound[:, k] -= 0.5 * np.sum(d * precs[k], axis=1)
    return bound


def _bound_state_loglik_tied(X, initial_bound, bound_prec, precs, means):
    n_components, n_features = means.shape
    n_samples = X.shape[0]
    bound = np.empty((n_samples, n_components))
    bound[:] = bound_prec + initial_bound
    # Transform the data to be able to apply standard Euclidean distance,
    # rather than Mahlanobis distance
    sqrt_cov = linalg.cholesky(precs)
    means = np.dot(means, sqrt_cov.T)
    X = np.dot(X, sqrt_cov.T)
    bound -= 0.5 * euclidean_distances(X, means, squared=True)
    return bound


def _bound_state_loglik_full(X, initial_bound, bound_prec, precs, means):
    n_components, n_features = means.shape
    n_samples = X.shape[0]
    bound = np.empty((n_samples, n_components))
    bound[:] = bound_prec + initial_bound
    for k in xrange(n_components):
        d = X - means[k]
        sqrt_cov = linalg.cholesky(precs[k])
        d = np.dot(d, sqrt_cov.T)
        d **= 2
        bound[:, k] -= 0.5 * d.sum(axis=-1)
    return bound


_BOUND_STATE_LOGLIK_DICT = dict(
    spherical=_bound_state_loglik_spherical,
    diag=_bound_state_loglik_diag,
    tied=_bound_state_loglik_tied,
    full=_bound_state_loglik_full)


class DPGMM(GMM):
    """Variational Inference for the Infinite Gaussian Mixture Model.

    DPGMM stands for Dirichlet Process Gaussian Mixture Model, and it
    is an infinite mixture model with the Dirichlet Process as a prior
    distribution on the number of clusters. In practice the
    approximate inference algorithm uses a truncated distribution with
    a fixed maximum number of components, but almost always the number
    of components actually used depends on the data.

    Stick-breaking Representation of a Gaussian mixture model
    probability distribution. This class allows for easy and efficient
    inference of an approximate posterior distribution over the
    parameters of a gaussian mixture model with a variable number of
    components (smaller than the truncation parameter n_components).

    Initialization is with normally-distributed means and identity
    covariance, for proper convergence.

    Parameters
    ----------
    n_components: int, optional
        Number of mixture components. Defaults to 1.

    covariance_type: string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    alpha: float, optional
        Real number representing the concentration parameter of
        the dirichlet process. Intuitively, the Dirichlet Process
        is as likely to start a new cluster for a point as it is
        to add that point to a cluster with alpha elements. A
        higher alpha means more clusters, as the expected number
        of clusters is alpha*log(N). Defaults to 1.

    thresh : float, optional
        Convergence threshold.

    Attributes
    ----------
    covariance_type : string (read-only)
        String describing the type of covariance parameters used by
        the DP-GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.

    n_components : int (read-only)
        Number of mixture components.

    weights : array, shape (`n_components`,)
        Mixing weights for each mixture component.

    means : array, shape (`n_components`, `n_features`)
        Mean parameters for each mixture component.

    precisions : array
        Precision (inverse covariance) parameters for each mixture
        component.  The shape depends on `covariance_type`::
            (`n_components`,)                             if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_components`, `n_features`)                if 'diag',
            (`n_components`, `n_features`, `n_features`)  if 'full'

    converged_ : bool
        True when convergence was reached in fit(), False
        otherwise.

    Methods
    -------
    eval(X)
        Compute a lower-bound of the log likelihood of X under the model
        and an approximate posterior distribution over mixture components.
    fit(X)
        Estimate the posterior of themodel parameters from X using the
        variational mean-field algorithm.
    predict(X)
        Find most likely mixtures components for each observation in X.
    rvs(n=1)
        Generate `n` samples from the posterior for the model.
    score(X)
        Compute the log likelihood of X under the model.

    See Also
    --------
    GMM : Finite gaussian mixture model fit with EM

    VBGMM : Finite gaussian mixture model fit with a variational
    algorithm, better for situations where there might be too little
    data to get a good estimate of the covariance matrix.

    """

    def __init__(self, n_components=1, covariance_type='diag', alpha=1.0,
                 random_state=None, thresh=1e-2, verbose=False,
                 min_covar=None):
        self.alpha = alpha
        self.verbose = verbose
        super(DPGMM, self).__init__(n_components, covariance_type,
                                    random_state=random_state,
                                    thresh=thresh, min_covar=min_covar)

    def _get_precisions(self):
        """Return precisions as a full matrix."""
        if self.covariance_type == 'full':
            return self.precs_
        elif self.covariance_type == 'diag':
            return [np.diag(cov) for cov in self.precs_]
        elif self.covariance_type == 'tied':
            return [self.precs_] * self.n_components
        elif self.covariance_type == 'spherical':
            return [np.eye(self.means.shape[1]) * f for f in self.precs_]

    def _get_covars(self):
        return [linalg.pinv(c) for c in self._get_precisions()]

    def _set_covars(self, covars):
        raise NotImplementedError("""The variational algorithm does
        not support setting the covariance parameters.""")

    precisions = property(_get_precisions, _set_covars)
    covars = property(_get_covars, _set_covars)

    def eval(self, X):
        """Evaluate the model on data

        Compute the bound on log probability of X under the model
        and return the posterior distribution (responsibilities) of
        each mixture component for each element of X.

        This is done by computing the parameters for the mean-field of
        z for each observation.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in X
        responsibilities: array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        X = np.asarray(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        z = np.zeros((X.shape[0], self.n_components))
        sd = digamma(self._gamma.T[1] + self._gamma.T[2])
        dgamma1 = digamma(self._gamma.T[1]) - sd
        dgamma2 = np.zeros(self.n_components)
        dgamma2[0] = digamma(self._gamma[0, 2]) - digamma(self._gamma[0, 1] +
                                                          self._gamma[0, 2])
        for j in xrange(1, self.n_components):
            dgamma2[j] = dgamma2[j - 1] + digamma(self._gamma[j - 1, 2])
            dgamma2[j] -= sd[j - 1]
        dgamma = dgamma1 + dgamma2
        # Free memory and developers cognitive load:
        del dgamma1, dgamma2, sd

        try:
            _bound_state_loglik = _BOUND_STATE_LOGLIK_DICT[self.covariance_type]
        except KeyError:
            raise NotImplementedError("This ctype is not implemented: %s"
                                      % self.covariance_type)

        p = _bound_state_loglik(X, self._initial_bound,
                        self._bound_prec, self.precs_, self.means_)
        z = p + dgamma
        z = log_normalize(z, axis=-1)
        bound = np.sum(z * p, axis=-1)
        return bound, z

    def _update_concentration(self, z):
        """Update the concentration parameters for each cluster"""
        sz = np.sum(z, axis=0)
        self._gamma.T[1] = 1. + sz
        self._gamma.T[2].fill(0)
        for i in xrange(self.n_components - 2, -1, -1):
            self._gamma[i, 2] = self._gamma[i + 1, 2] + sz[i]
        self._gamma.T[2] += self.alpha

    def _update_means(self, X, z):
        """Update the variational distributions for the means"""
        n_features = X.shape[1]
        for k in xrange(self.n_components):
            if self.covariance_type == 'spherical' or self.covariance_type == 'diag':
                num = np.sum(z.T[k].reshape((-1, 1)) * X, axis=0)
                num *= self.precs_[k]
                den = 1. + self.precs_[k] * np.sum(z.T[k])
                self.means_[k] = num / den
            elif self.covariance_type == 'tied' or self.covariance_type == 'full':
                if self.covariance_type == 'tied':
                    cov = self.precs_
                else:
                    cov = self.precs_[k]
                den = np.identity(n_features) + cov * np.sum(z.T[k])
                num = np.sum(z.T[k].reshape((-1, 1)) * X, axis=0)
                num = np.dot(cov, num)
                self.means_[k] = linalg.lstsq(den, num)[0]

    def _update_precisions(self, X, z):
        """Update the variational distributions for the precisions"""
        n_features = X.shape[1]
        if self.covariance_type == 'spherical':
            self._a = 0.5 * n_features * np.sum(z, axis=0)
            for k in xrange(self.n_components):
                # XXX: how to avoid this huge temporary matrix in memory
                dif = (X - self.means_[k])
                self._b[k] = 1.
                d = np.sum(dif * dif, axis=1)
                self._b[k] += 0.5 * np.sum(z.T[k] * (d + n_features))
                self._bound_prec[k] = (
                    0.5 * n_features * (
                        digamma(self._a[k]) - np.log(self._b[k])))
            self.precs_ = self._a / self._b

        elif self.covariance_type == 'diag':
            for k in xrange(self.n_components):
                self._a[k].fill(1. + 0.5 * np.sum(z.T[k], axis=0))
                ddif = (X - self.means_[k])  # see comment above
                for d in xrange(n_features):
                    self._b[k, d] = 1.
                    dd = ddif.T[d] * ddif.T[d]
                    self._b[k, d] += 0.5 * np.sum(z.T[k] * (dd + 1))
                self.precs_[k] = self._a[k] / self._b[k]
                self._bound_prec[k] = 0.5 * np.sum(digamma(self._a[k])
                                                    - np.log(self._b[k]))
                self._bound_prec[k] -= 0.5 * np.sum(self.precs_[k])

        elif self.covariance_type == 'tied':
            self._a = 2 + X.shape[0] + n_features
            self._B = (X.shape[0] + 1) * np.identity(n_features)
            for i in xrange(X.shape[0]):
                for k in xrange(self.n_components):
                    dif = X[i] - self.means_[k]
                    self._B += z[i, k] * np.dot(dif.reshape((-1, 1)),
                                                dif.reshape((1, -1)))
            self._B = linalg.pinv(self._B)
            self.precs_ = self._a * self._B
            self._detB = linalg.det(self._B)
            self._bound_prec = 0.5 * detlog_wishart(
                self._a, self._B, self._detB, n_features)
            self._bound_prec -= 0.5 * self._a * np.trace(self._B)

        elif self.covariance_type == 'full':
            for k in xrange(self.n_components):
                T = np.sum(z.T[k])
                self._a[k] = 2 + T + n_features
                self._B[k] = (T + 1) * np.identity(n_features)
                for i in xrange(X.shape[0]):
                    dif = X[i] - self.means_[k]
                    self._B[k] += z[i, k] * np.dot(dif.reshape((-1, 1)),
                                                   dif.reshape((1, -1)))
                self._B[k] = linalg.pinv(self._B[k])
                self.precs_[k] = self._a[k] * self._B[k]
                self._detB[k] = linalg.det(self._B[k])
                self._bound_prec[k] = 0.5 * detlog_wishart(self._a[k],
                                                           self._B[k],
                                                           self._detB[k],
                                                           n_features)
                self._bound_prec[k] -= 0.5 * self._a[k] * np.trace(self._B[k])

    def _monitor(self, X, z, n, end=False):
        """Monitor the lower bound during iteration

        Debug method to help see exactly when it is failing to converge as
        expected.

        Note: this is very expensive and should not be used by default."""
        if self.verbose:
            print "Bound after updating %8s: %f" % (n, self.lower_bound(X, z))
            if end == True:
                print "Cluster proportions:", self._gamma.T[1]
                print "covariance_type:", self._covariance_type

    def _do_mstep(self, X, z, params):
        """Maximize the variational lower bound

        Update each of the parameters to maximize the lower bound."""
        self._monitor(X, z, "z")
        self._update_concentration(z)
        self._monitor(X, z, "gamma")
        if 'm' in params:
            self._update_means(X, z)
        self._monitor(X, z, "mu")
        if 'c' in params:
            self._update_precisions(X, z)
        self._monitor(X, z, "a and b", end=True)

    def _initialize_gamma(self):
        "Initializes the concentration parameters"
        self._gamma = self.alpha * np.ones((self.n_components, 3))

    def _bound_concentration(self):
        "The variational lower bound for the concentration parameter."
        logprior = 0.
        for k in xrange(self.n_components):
            logprior = gammaln(self.alpha)
            logprior += (self.alpha - 1) * (digamma(self._gamma[k, 2]) -
                                            digamma(self._gamma[k, 1] +
                                                    self._gamma[k, 2]))
            logprior += -gammaln(self._gamma[k, 1] + self._gamma[k, 2])
            logprior += gammaln(self._gamma[k, 1]) + gammaln(self._gamma[k, 2])
            logprior -= (self._gamma[k, 1] - 1) * (digamma(self._gamma[k, 1]) -
                                                   digamma(self._gamma[k, 1] +
                                                           self._gamma[k, 2]))
            logprior -= (self._gamma[k, 2] - 1) * (digamma(self._gamma[k, 2]) -
                                                   digamma(self._gamma[k, 1] +
                                                           self._gamma[k, 2]))
        return logprior

    def _bound_means(self):
        "The variational lower bound for the mean parameters"
        logprior = 0.
        logprior -= 0.5 * sqnorm(self.means_)
        logprior -= 0.5 * self.means.shape[1] * self.n_components
        return logprior

    def _bound_wishart(self, a, B, detB):
        n_features = self.means.shape[1]
        logprior = wishart_logz(a, B, detB, n_features)
        logprior -= wishart_logz(n_features,
                                 np.identity(n_features),
                                 1, n_features)
        logprior += 0.5 * (a - 1) * detlog_wishart(a, B, detB, n_features)
        logprior += 0.5 * a * np.trace(B)
        return logprior

    def _bound_precisions(self):
        logprior = 0.
        if self.covariance_type == 'spherical':
            for k in xrange(self.n_components):
                logprior += gammaln(self._a[k])
                logprior -= (self._a[k] - 1) * digamma(max(0.5, self._a[k]))
                logprior += - np.log(self._b[k]) + self._a[k] - self.precs_[k]
        elif self.covariance_type == 'diag':
            for k in xrange(self.n_components):
                for d in xrange(self.means.shape[1]):
                    logprior += gammaln(self._a[k, d])
                    logprior -= (self._a[k, d] - 1) * digamma(self._a[k, d])
                    logprior -= np.log(self._b[k, d])
                    logprior += self._a[k, d] - self.precs_[k, d]
        elif self.covariance_type == 'tied':
            logprior += self._bound_wishart(self._a, self._B, self._detB)
        elif self.covariance_type == 'full':
            for k in xrange(self.n_components):
                logprior += self._bound_wishart(self._a[k],
                                                self._B[k],
                                                self._detB[k])
        return logprior

    def _bound_proportions(self, z):
        dg12 = digamma(self._gamma.T[1] + self._gamma.T[2])
        dg1 = digamma(self._gamma.T[1]) - dg12
        dg2 = digamma(self._gamma.T[2]) - dg12

        cz = np.cumsum(z[:, ::-1], axis=-1)[:, -2::-1]
        logprior = np.sum(cz * dg2[:-1]) + np.sum(z * dg1)
        del cz  # Save memory
        z_non_zeros = z[z > np.finfo(np.float32).eps]
        logprior -= np.sum(z_non_zeros * np.log(z_non_zeros))
        return logprior

    def _logprior(self, z):
        logprior = self._bound_concentration()
        logprior += self._bound_means()
        logprior += self._bound_precisions()
        logprior += self._bound_proportions(z)
        return logprior

    def lower_bound(self, X, z):
        try:
            _bound_state_loglik = _BOUND_STATE_LOGLIK_DICT[self.covariance_type]
        except KeyError:
            raise NotImplementedError("This ctype is not implemented: %s"
                                      % self.covariance_type)
        X = np.asarray(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        c = np.sum(z * _bound_state_loglik(
                X, self._initial_bound, self._bound_prec, self.precs_, 
                self.means_))

        return c + self._logprior(z)

    def fit(self, X, n_iter=10, params='wmc', init_params='wmc'):
        """Estimate model parameters with the variational
        algorithm.

        For a full derivation and description of the algorithm see
        doc/dp-derivation/dp-derivation.tex

        A initialization step is performed before entering the em
        algorithm. If you want to avoid this step, set the keyword
        argument init_params to the empty string ''. Likewise, if you
        would like just to do an initialization, call this method with
        n_iter=0.

        Parameters
        ----------
        X : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.
        n_iter : int, optional
             Maximum number of iterations to perform before convergence.
       params : string, optional
            Controls which parameters are updated in the training
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
       init_params : string, optional
            Controls which parameters are updated in the initialization
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
        """
        self.random_state = check_random_state(self.random_state)

        ## initialization step
        X = np.asarray(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]

        n_features = X.shape[1]
        z = np.ones((X.shape[0], self.n_components))
        z /= self.n_components

        self._initial_bound = -0.5 * n_features * np.log(2 * np.pi)
        self._initial_bound -= np.log(2 * np.pi * np.e)

        if init_params != '':
            self._initialize_gamma()

        if 'm' in init_params or not hasattr(self, 'means'):
            self.means_ = cluster.KMeans(
                k=self.n_components, random_state=self.random_state
            ).fit(X).cluster_centers_[::-1]

        if 'w' in init_params or not hasattr(self, 'weights'):
            self.weights = np.tile(1.0 / self.n_components, self.n_components)

        if 'c' in init_params or not hasattr(self, 'covars'):
            if self.covariance_type == 'spherical':
                self._a = np.ones(self.n_components)
                self._b = np.ones(self.n_components)
                self.precs_ = np.ones(self.n_components)
                self._bound_prec = (0.5 * n_features *
                                     (digamma(self._a) -
                                      np.log(self._b)))
            elif self.covariance_type == 'diag':
                self._a = 1 + 0.5 * n_features
                self._a *= np.ones((self.n_components, n_features))
                self._b = np.ones((self.n_components, n_features))
                self.precs_ = np.ones((self.n_components, n_features))
                self._bound_prec = np.zeros(self.n_components)
                for k in xrange(self.n_components):
                    self._bound_prec[k] = 0.5 * np.sum(digamma(self._a[k])
                                                        - np.log(self._b[k]))
                    self._bound_prec[k] -= 0.5 * np.sum(self.precs_[k])
            elif self.covariance_type == 'tied':
                self._a = 1.
                self._B = np.identity(n_features)
                self.precs_ = np.identity(n_features)
                self._detB = 1.
                self._bound_prec = 0.5 * detlog_wishart(
                    self._a, self._B, self._detB, n_features)
                self._bound_prec -= 0.5 * self._a * np.trace(self._B)
            elif self.covariance_type == 'full':
                self._a = (1 + self.n_components + X.shape[0])
                self._a *= np.ones(self.n_components)
                self._B = [2 * np.identity(n_features)
                           for i in xrange(self.n_components)]
                self.precs_ = [np.identity(n_features)
                                for i in xrange(self.n_components)]
                self._detB = np.ones(self.n_components)
                self._bound_prec = np.zeros(self.n_components)
                for k in xrange(self.n_components):
                    self._bound_prec[k] = detlog_wishart(
                        self._a[k], self._B[k], self._detB[k], n_features)
                    self._bound_prec[k] -= self._a[k] * np.trace(self._B[k])
                    self._bound_prec[k] *= 0.5

        logprob = []
        # reset self.converged_ to False
        self.converged_ = False
        for i in xrange(n_iter):
            # Expectation step
            curr_logprob, _ = self.eval(X)
            logprob.append(curr_logprob.sum() + self._logprior(z))

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < self.thresh:
                self.converged_ = True
                break

            # Maximization step
            self._do_mstep(X, z, params)

        return self


class VBGMM(DPGMM):
    """Variational Inference for the Gaussian Mixture Model

    Variational inference for a Gaussian mixture model probability
    distribution. This class allows for easy and efficient inference
    of an approximate posterior distribution over the parameters of a
    gaussian mixture model with a fixed number of components.

    Initialization is with normally-distributed means and identity
    covariance, for proper convergence.

    Parameters
    ----------
    n_components: int, optional
        Number of mixture components. Defaults to 1.

    covariance_type: string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    alpha: float, optional
        Real number representing the concentration parameter of
        the dirichlet distribution. Intuitively, the higher the
        value of alpha the more likely the variational mixture of
        gaussians model will use all components it can. Defaults
        to 1.


    Attributes
    ----------
    covariance_type : string (read-only)
        String describing the type of covariance parameters used by
        the DP-GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    n_features : int
        Dimensionality of the Gaussians.
    n_components : int (read-only)
        Number of mixture components.
    weights : array, shape (`n_components`,)
        Mixing weights for each mixture component.
    means : array, shape (`n_components`, `n_features`)
        Mean parameters for each mixture component.
    precisions : array
        Precision (inverse covariance) parameters for each mixture
        component.  The shape depends on `covariance_type`:
            (`n_components`,)                             if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_components`, `n_features`)                if 'diag',
            (`n_components`, `n_features`, `n_features`)  if 'full'
    converged_ : bool
        True when convergence was reached in fit(), False
        otherwise.

    Methods
    -------
    eval(X)
        Compute a lower-bound of the log likelihood of X under the model
        and an approximate posterior distribution over mixture components.
    fit(X)
        Estimate the posterior of themodel parameters from X using the
        variational mean-field algorithm.
    predict(X)
        Find most likely mixtures components for each observation in X.
    rvs(n=1)
        Generate `n` samples from the posterior for the model.
    score(X)
        Compute the log likelihood of X under the model.


    See Also
    --------
    GMM : Finite gaussian mixture model fit with EM

    DPGMM : Ininite gaussian mixture model, using the dirichlet
    process, fit with a variational algorithm
    """

    def __init__(self, n_components=1, covariance_type='diag', alpha=1.0,
                 random_state=None, thresh=1e-2, verbose=False,
                 min_covar=None):
        super(VBGMM, self).__init__(
            n_components, covariance_type, random_state=random_state, thresh=thresh,
            verbose=verbose, min_covar=min_covar)
        self.alpha = float(alpha) / n_components

    def eval(self, X):
        """Evaluate the model on data

        Compute the bound on log probability of X under the model
        and return the posterior distribution (responsibilities) of
        each mixture component for each element of X.

        This is done by computing the parameters for the mean-field of
        z for each observation.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in X
        responsibilities: array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        X = np.asarray(X)
        if X.ndim == 1:
            X = X[:, np.newaxis]
        z = np.zeros((X.shape[0], self.n_components))
        p = np.zeros(self.n_components)
        bound = np.zeros(X.shape[0])
        dg = digamma(self._gamma) - digamma(np.sum(self._gamma))
        try:
            _bound_state_loglik = _BOUND_STATE_LOGLIK_DICT[self.covariance_type]
        except KeyError:
            raise NotImplementedError("This ctype is not implemented: %s"
                                      % self.covariance_type)

        p = _bound_state_loglik(X, self._initial_bound,
                                self._bound_prec, self.precs_, self.means_)
        z = p + dg
        z = log_normalize(z, axis=-1)
        bound = np.sum(z * p, axis=-1)
        return bound, z

    def _update_concentration(self, z):
        for i in xrange(self.n_components):
            self._gamma[i] = self.alpha + np.sum(z.T[i])

    def _initialize_gamma(self):
        self._gamma = self.alpha * np.ones(self.n_components)

    def _bound_proportions(self, z):
        logprior = 0.
        dg = digamma(self._gamma)
        dg -= digamma(np.sum(self._gamma))
        logprior += np.sum(dg.reshape((-1, 1)) * z.T)
        z_non_zeros = z[z > np.finfo(np.float32).eps]
        logprior -= np.sum(z_non_zeros * np.log(z_non_zeros))
        return logprior

    def _bound_concentration(self):
        logprior = 0.
        logprior = gammaln(np.sum(self._gamma)) - gammaln(self.n_components
                                                          * self.alpha)
        logprior -= np.sum(gammaln(self._gamma) - gammaln(self.alpha))
        sg = digamma(np.sum(self._gamma))
        logprior += np.sum((self._gamma - self.alpha)
                           * (digamma(self._gamma) - sg))
        return logprior

    def _monitor(self, X, z, n, end=False):
        """Monitor the lower bound during iteration

        Debug method to help see exactly when it is failing to converge as
        expected.

        Note: this is very expensive and should not be used by default."""
        if self.verbose:
            print "Bound after updating %8s: %f" % (n, self.lower_bound(X, z))
            if end == True:
                print "Cluster proportions:", self._gamma
                print "covariance_type:", self._covariance_type
