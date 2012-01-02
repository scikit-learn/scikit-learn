"""Dirichlet Process Gaussian Mixture Models"""

# Author: Alexandre Passos (alexandre.tp@gmail.com)
#
# Based on mixture.py by:
#         Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#

import numpy as np
from scipy.special import digamma as _digamma, gammaln as _gammaln
from scipy import linalg
from scipy.spatial.distance import cdist

from ..utils import check_random_state
from ..utils.extmath import norm
from .. import cluster
from ..metrics import euclidean_distances
from .gmm import GMM


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


def _sym_quad_form(x, mu, A):
    """helper function to calculate symmetric quadratic form x.T * A * x"""
    q = (cdist(x, mu[np.newaxis], "mahalanobis", VI=A) ** 2).reshape(-1)
    return q


def _bound_state_loglik_full(X, initial_bound, bound_prec, precs, means):
    n_components, n_features = means.shape
    n_samples = X.shape[0]
    bound = np.empty((n_samples, n_components))
    bound[:] = bound_prec + initial_bound
    for k in xrange(n_components):
        bound[:, k] -= 0.5 * _sym_quad_form(X, means[k], precs[k])
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
    parameters of a Gaussian mixture model with a variable number of
    components (smaller than the truncation parameter n_components).

    Initialization is with normally-distributed means and identity
    covariance, for proper convergence.

    Parameters
    ----------
    n_components: int, optional
        Number of mixture components. Defaults to 1.

    cvtype: string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    alpha: float, optional
        Real number representing the concentration parameter of
        the dirichlet process. Intuitively, the Dirichler Process
        is as likely to start a new cluster for a point as it is
        to add that point to a cluster with alpha elements. A
        higher alpha means more clusters, as the expected number
        of clusters is ``alpha*log(N)``. Defaults to 1.

    thresh : float, optional
        Convergence threshold.

    Attributes
    ----------
    cvtype : string (read-only)
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
        component.  The shape depends on `cvtype`::

            (`n_components`,)                             if 'spherical',
            (`n_features`, `n_features`)                  if 'tied',
            (`n_components`, `n_features`)                if 'diag',
            (`n_components`, `n_features`, `n_features`)  if 'full'

    `converged_` : bool
        True when convergence was reached in fit(), False otherwise.

    See Also
    --------
    GMM : Finite Gaussian mixture model fit with EM

    VBGMM : Finite Gaussian mixture model fit with a variational
    algorithm, better for situations where there might be too little
    data to get a good estimate of the covariance matrix.

    """

    def __init__(self, n_components=1, cvtype='diag', alpha=1.0,
                 random_state=None, thresh=1e-2, verbose=False,
                 min_covar=None):
        self.alpha = alpha
        self.verbose = verbose
        super(DPGMM, self).__init__(n_components, cvtype,
                                    random_state=random_state,
                                    thresh=thresh, min_covar=min_covar)

    def _get_precisions(self):
        """Return precisions as a full matrix."""
        if self.cvtype == 'full':
            return self._precs
        elif self.cvtype == 'diag':
            return [np.diag(cov) for cov in self._precs]
        elif self.cvtype == 'tied':
            return [self._precs] * self.n_components
        elif self.cvtype == 'spherical':
            return [np.eye(self.n_features) * f for f in self._precs]

    def _get_covars(self):
        return [linalg.pinv(c) for c in self._get_precisions()]

    def _set_covars(self, covars):
        raise NotImplementedError("""The variational algorithm does
        not support setting the covariance parameters.""")

    precisions = property(_get_precisions, _set_covars)
    covars = property(_get_covars, _set_covars)

    def eval(self, obs=None):
        """Evaluate the model on data

        Compute the bound on log probability of `obs` under the model
        and return the posterior distribution (responsibilities) of
        each mixture component for each element of `obs`.

        This is done by computing the parameters for the mean-field of
        z for each observation.

        Parameters
        ----------
        obs : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in `obs`
        posteriors: array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        if obs is None:
            z = self._z
            obs = self._X
        else:
            z = np.zeros((obs.shape[0], self.n_components))
        obs = np.asarray(obs)
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
            _bound_state_loglik = _BOUND_STATE_LOGLIK_DICT[self.cvtype]
        except KeyError:
            raise NotImplementedError("This ctype is not implemented: %s"
                                      % self.cvtype)

        p = _bound_state_loglik(obs, self._initial_bound,
                        self._bound_prec, self._precs, self._means)
        z = p + dgamma
        self._z = z = log_normalize(z, axis=-1)
        bound = np.sum(z * p, axis=-1)
        return bound, z

    def _update_concentration(self):
        """Update the concentration parameters for each cluster"""
        sz = np.sum(self._z, axis=0)
        self._gamma.T[1] = 1. + sz
        self._gamma.T[2].fill(0)
        for i in xrange(self.n_components - 2, -1, -1):
            self._gamma[i, 2] = self._gamma[i + 1, 2] + sz[i]
        self._gamma.T[2] += self.alpha

    def _update_means(self):
        """Update the variational distributions for the means"""
        for k in xrange(self.n_components):
            if self.cvtype == 'spherical' or self.cvtype == 'diag':
                num = np.sum(self._z.T[k].reshape((-1, 1)) * self._X, axis=0)
                num *= self._precs[k]
                den = 1. + self._precs[k] * np.sum(self._z.T[k])
                self._means[k] = num / den
            elif self.cvtype == 'tied' or self.cvtype == 'full':
                if self.cvtype == 'tied':
                    cov = self._precs
                else:
                    cov = self._precs[k]
                den = np.identity(self.n_features) + cov * np.sum(self._z.T[k])
                num = np.sum(self._z.T[k].reshape((-1, 1)) * self._X, axis=0)
                num = np.dot(cov, num)
                self._means[k] = linalg.lstsq(den, num)[0]

    def _update_precisions(self):
        """Update the variational distributions for the precisions"""
        if self.cvtype == 'spherical':
            self._a = 0.5 * self.n_features * np.sum(self._z, axis=0)
            for k in xrange(self.n_components):
                # XXX: how to avoid this huge temporary matrix in memory
                dif = (self._X - self._means[k])
                self._b[k] = 1.
                d = np.sum(dif * dif, axis=1)
                self._b[k] += 0.5 * np.sum(
                    self._z.T[k] * (d + self.n_features))
                self._bound_prec[k] = (
                    0.5 * self.n_features * (
                        digamma(self._a[k]) - np.log(self._b[k])))
            self._precs = self._a / self._b

        elif self.cvtype == 'diag':
            for k in xrange(self.n_components):
                self._a[k].fill(1. + 0.5 * np.sum(self._z.T[k], axis=0))
                ddif = (self._X - self._means[k])  # see comment above
                for d in xrange(self.n_features):
                    self._b[k, d] = 1.
                    dd = ddif.T[d] * ddif.T[d]
                    self._b[k, d] += 0.5 * np.sum(self._z.T[k] * (dd + 1))
                self._precs[k] = self._a[k] / self._b[k]
                self._bound_prec[k] = 0.5 * np.sum(digamma(self._a[k])
                                                    - np.log(self._b[k]))
                self._bound_prec[k] -= 0.5 * np.sum(self._precs[k])

        elif self.cvtype == 'tied':
            self._a = 2 + self._X.shape[0] + self.n_features
            self._B = (self._X.shape[0] + 1) * np.identity(self.n_features)
            for i in xrange(self._X.shape[0]):
                for k in xrange(self.n_components):
                    dif = self._X[i] - self._means[k]
                    self._B += self._z[i, k] * np.dot(dif.reshape((-1, 1)),
                                                      dif.reshape((1, -1)))
            self._B = linalg.pinv(self._B)
            self._precs = self._a * self._B
            self._detB = linalg.det(self._B)
            self._bound_prec = 0.5 * detlog_wishart(
                self._a, self._B, self._detB, self.n_features)
            self._bound_prec -= 0.5 * self._a * np.trace(self._B)

        elif self.cvtype == 'full':
            for k in xrange(self.n_components):
                T = np.sum(self._z.T[k])
                self._a[k] = 2 + T + self.n_features
                self._B[k] = (T + 1) * np.identity(self.n_features)
                dx = self._X - self._means[k]
                self._B[k] += np.dot((self._z[:, k] * dx.T), dx)
                self._B[k] = linalg.inv(self._B[k])
                self._precs[k] = self._a[k] * self._B[k]
                self._detB[k] = linalg.det(self._B[k])
                self._bound_prec[k] = 0.5 * detlog_wishart(self._a[k],
                                                           self._B[k],
                                                           self._detB[k],
                                                           self.n_features)
                self._bound_prec[k] -= 0.5 * self._a[k] * np.trace(self._B[k])

    def _monitor(self, n, end=False):
        """Monitor the lower bound during iteration

        Debug method to help see exactly when it is failing to converge as
        expected.

        Note: this is very expensive and should not be used by default."""
        if self.verbose:
            print "Bound after updating %8s: %f" % (n, self.lower_bound())
            if end == True:
                print "Cluster proportions:", self._gamma.T[1]
                print "cvtype:", self._cvtype

    def _do_mstep(self, params):
        """Maximize the variational lower bound

        Update each of the parameters to maximize the lower bound."""
        self._monitor("z")
        self._update_concentration()
        self._monitor("gamma")
        if 'm' in params:
            self._update_means()
        self._monitor("mu")
        if 'c' in params:
            self._update_precisions()
        self._monitor("a and b", end=True)

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
        logprior -= 0.5 * sqnorm(self._means)
        logprior -= 0.5 * self.n_features * self.n_components
        return logprior

    def _bound_wishart(self, a, B, detB):
        logprior = wishart_logz(a, B, detB, self.n_features)
        logprior -= wishart_logz(self.n_features,
                                 np.identity(self.n_features),
                                 1, self.n_features)
        logprior += 0.5 * (a - 1) * detlog_wishart(a, B, detB, self.n_features)
        logprior += 0.5 * a * np.trace(B)
        return logprior

    def _bound_precisions(self):
        logprior = 0.
        if self.cvtype == 'spherical':
            for k in xrange(self.n_components):
                logprior += gammaln(self._a[k])
                logprior -= (self._a[k] - 1) * digamma(max(0.5, self._a[k]))
                logprior += - np.log(self._b[k]) + self._a[k] - self._precs[k]
        elif self.cvtype == 'diag':
            for k in xrange(self.n_components):
                for d in xrange(self.n_features):
                    logprior += gammaln(self._a[k, d])
                    logprior -= (self._a[k, d] - 1) * digamma(self._a[k, d])
                    logprior -= np.log(self._b[k, d])
                    logprior += self._a[k, d] - self._precs[k, d]
        elif self.cvtype == 'tied':
            logprior += self._bound_wishart(self._a, self._B, self._detB)
        elif self.cvtype == 'full':
            for k in xrange(self.n_components):
                logprior += self._bound_wishart(self._a[k],
                                                self._B[k],
                                                self._detB[k])
        return logprior

    def _bound_proportions(self):
        dg12 = digamma(self._gamma.T[1] + self._gamma.T[2])
        dg1 = digamma(self._gamma.T[1]) - dg12
        dg2 = digamma(self._gamma.T[2]) - dg12

        cz = np.cumsum(self._z[:, ::-1], axis=-1)[:, -2::-1]
        logprior = np.sum(cz * dg2[:-1]) + np.sum(self._z * dg1)
        del cz  # Save memory
        z_non_zeros = self._z[self._z > np.finfo(np.float32).eps]
        logprior -= np.sum(z_non_zeros * np.log(z_non_zeros))
        return logprior

    def _logprior(self):
        logprior = self._bound_concentration()
        logprior += self._bound_means()
        logprior += self._bound_precisions()
        logprior += self._bound_proportions()
        return logprior

    def lower_bound(self):
        try:
            _bound_state_loglik = _BOUND_STATE_LOGLIK_DICT[self.cvtype]
        except KeyError:
            raise NotImplementedError("This ctype is not implemented: %s"
                                      % self.cvtype)

        c = np.sum(self._z * _bound_state_loglik(self._X, self._initial_bound,
                        self._bound_prec, self._precs, self._means))

        return c + self._logprior()

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

        self._X = np.asarray(X)
        if hasattr(self, 'n_features') and self.n_features != self._X.shape[1]:
            raise ValueError('Unexpected number of dimensions, got %s but '
                             'expected %s' % (self._X.shape[1],
                                              self.n_features))

        self.n_features = self._X.shape[1]
        self._z = np.ones((self._X.shape[0], self.n_components))
        self._z /= self.n_components

        self._initial_bound = -0.5 * self.n_features * np.log(2 * np.pi)
        self._initial_bound -= np.log(2 * np.pi * np.e)

        if init_params != '':
            self._initialize_gamma()

        if 'm' in init_params or not hasattr(self, 'means'):
            self._means = cluster.KMeans(
                k=self.n_components, random_state=self.random_state
            ).fit(X).cluster_centers_[::-1]

        if 'w' in init_params or not hasattr(self, 'weights'):
            self.weights = np.tile(1.0 / self.n_components, self.n_components)

        if 'c' in init_params or not hasattr(self, 'covars'):
            if self.cvtype == 'spherical':
                self._a = np.ones(self.n_components)
                self._b = np.ones(self.n_components)
                self._precs = np.ones(self.n_components)
                self._bound_prec = (0.5 * self.n_features *
                                     (digamma(self._a) -
                                      np.log(self._b)))
            elif self.cvtype == 'diag':
                self._a = 1 + 0.5 * self.n_features
                self._a *= np.ones((self.n_components, self.n_features))
                self._b = np.ones((self.n_components, self.n_features))
                self._precs = np.ones((self.n_components, self.n_features))
                self._bound_prec = np.zeros(self.n_components)
                for k in xrange(self.n_components):
                    self._bound_prec[k] = 0.5 * np.sum(digamma(self._a[k])
                                                        - np.log(self._b[k]))
                    self._bound_prec[k] -= 0.5 * np.sum(self._precs[k])
            elif self.cvtype == 'tied':
                self._a = 1.
                self._B = np.identity(self.n_features)
                self._precs = np.identity(self.n_features)
                self._detB = 1.
                self._bound_prec = 0.5 * detlog_wishart(
                    self._a, self._B, self._detB, self.n_features)
                self._bound_prec -= 0.5 * self._a * np.trace(self._B)
            elif self.cvtype == 'full':
                self._a = (1 + self.n_components + self._X.shape[0])
                self._a *= np.ones(self.n_components)
                self._B = [2 * np.identity(self.n_features)
                           for i in xrange(self.n_components)]
                self._precs = [np.identity(self.n_features)
                                for i in xrange(self.n_components)]
                self._detB = np.ones(self.n_components)
                self._bound_prec = np.zeros(self.n_components)
                for k in xrange(self.n_components):
                    self._bound_prec[k] = detlog_wishart(
                        self._a[k], self._B[k], self._detB[k], self.n_features)
                    self._bound_prec[k] -= self._a[k] * np.trace(self._B[k])
                    self._bound_prec[k] *= 0.5

        logprob = []
        # reset self.converged_ to False
        self.converged_ = False
        for i in xrange(n_iter):
            # Expectation step
            curr_logprob, _ = self.eval()
            logprob.append(curr_logprob.sum() + self._logprior())

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < self.thresh:
                self.converged_ = True
                break

            # Maximization step
            self._do_mstep(params)

        return self


class VBGMM(DPGMM):
    """Variational Inference for the Gaussian Mixture Model

    Variational inference for a Gaussian mixture model probability
    distribution. This class allows for easy and efficient inference
    of an approximate posterior distribution over the parameters of a
    Gaussian mixture model with a fixed number of components.

    Initialization is with normally-distributed means and identity
    covariance, for proper convergence.

    Parameters
    ----------
    n_components: int, optional
        Number of mixture components. Defaults to 1.

    cvtype: string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    alpha: float, optional
        Real number representing the concentration parameter of
        the dirichlet distribution. Intuitively, the higher the
        value of alpha the more likely the variational mixture of
        Gaussians model will use all components it can. Defaults
        to 1.


    Attributes
    ----------
    cvtype : string (read-only)
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
        component.  The shape depends on `cvtype`::
            (`n_components`,)                             if 'spherical',
            (`n_features`, `n_features`)                  if 'tied',
            (`n_components`, `n_features`)                if 'diag',
            (`n_components`, `n_features`, `n_features`)  if 'full'

    `converged_` : bool
        True when convergence was reached in fit(), False
        otherwise.

    See Also
    --------
    GMM : Finite Gaussian mixture model fit with EM

    DPGMM : Ininite Gaussian mixture model, using the dirichlet
    process, fit with a variational algorithm
    """

    def __init__(self, n_components=1, cvtype='diag', alpha=1.0,
                 random_state=None, thresh=1e-2, verbose=False,
                 min_covar=None):
        super(VBGMM, self).__init__(
            n_components, cvtype, random_state=random_state, thresh=thresh,
            verbose=verbose, min_covar=min_covar)
        self.alpha = float(alpha) / n_components

    def eval(self, obs=None):
        """Evaluate the model on data

        Compute the bound on log probability of `obs` under the model
        and return the posterior distribution (responsibilities) of
        each mixture component for each element of `obs`.

        This is done by computing the parameters for the mean-field of
        z for each observation.

        Parameters
        ----------
        obs : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in `obs`

        posteriors: array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        if obs is None:
            z = self._z
            obs = self._X
        else:
            z = np.zeros((obs.shape[0], self.n_components))
        obs = np.asarray(obs)
        p = np.zeros(self.n_components)
        bound = np.zeros(obs.shape[0])
        dg = digamma(self._gamma) - digamma(np.sum(self._gamma))
        try:
            _bound_state_loglik = _BOUND_STATE_LOGLIK_DICT[self.cvtype]
        except KeyError:
            raise NotImplementedError("This ctype is not implemented: %s"
                                      % self.cvtype)

        p = _bound_state_loglik(obs, self._initial_bound,
                                self._bound_prec, self._precs, self._means)
        z = p + dg
        self._z = z = log_normalize(z, axis=-1)
        bound = np.sum(z * p, axis=-1)
        return bound, z

    def _update_concentration(self):
        for i in xrange(self.n_components):
            self._gamma[i] = self.alpha + np.sum(self._z.T[i])

    def _initialize_gamma(self):
        self._gamma = self.alpha * np.ones(self.n_components)

    def _bound_proportions(self):
        logprior = 0.
        dg = digamma(self._gamma)
        dg -= digamma(np.sum(self._gamma))
        logprior += np.sum(dg.reshape((-1, 1)) * self._z.T)
        z_non_zeros = self._z[self._z > np.finfo(np.float32).eps]
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

    def _monitor(self, n, end=False):
        """Monitor the lower bound during iteration

        Debug method to help see exactly when it is failing to converge as
        expected.

        Note: this is very expensive and should not be used by default."""
        if self.verbose:
            print "Bound after updating %8s: %f" % (n, self.lower_bound())
            if end == True:
                print "Cluster proportions:", self._gamma
                print "cvtype:", self._cvtype
