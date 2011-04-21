import numpy as np

from . import gmm as mixture

from scikits.learn import cluster

from scipy.special import digamma, gammaln

from scipy import linalg

# Author: Alexandre Passos (alexandre.tp@gmail.com)
#
# Based on mixture.py by:
#         Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#


def norm(v):
    """The squared norm of vector v, as a scalar. As dot is not always
    scalar, using np.sum."""
    return np.sum(v * v)


def diagnorm(x, Sigma):
    """The x^T Sigma x norm when x is a vector representing a diagonal
    matrix."""
    return np.sum(x * Sigma * x)


def squarenorm(x, Sigma):
    """The x^T Sigma x norm when Sigma is a matrix."""
    return np.sum(np.dot(np.dot(x, Sigma), x))


def lognormalize(v):
    """Given a vector of unnormalized log-probabilites v returns a
 vector of normalized probabilities"""
    v = np.exp(v - np.logaddexp.reduce(v))
    return v / np.sum(v)


class DPGMM(mixture.GMM):
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
    components (smaller than the truncation parameter n_states).

    Initialization is with normally-distributed means and identity
    covariance, for proper convergence.

    Parameters
    ----------
    n_states : int, optional
        Number of mixture components. Defaults to 1.

    cvtype : string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    alpha : float, optional
        Real number representing the concentration parameter of
        the dirichlet process. Intuitively, the Dirichler Process
        is as likely to start a new cluster for a point as it is 
        to add that point to a cluster with alpha elements. A
        higher alpha means more clusters, as the expected number 
        of clusters is alpha*log(N). Defaults to 1.


    Attributes
    ----------
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the DP-GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    n_features : int
        Dimensionality of the Gaussians.
    n_states : int (read-only)
        Number of mixture components.
    weights : array, shape (`n_states`,)
        Mixing weights for each mixture component.
    means : array, shape (`n_states`, `n_features`)
        Mean parameters for each mixture component.
    precisions : array
        Precision (inverse covariance) parameters for each mixture
        component.  The shape depends on `cvtype`:
            (`n_states`,)                             if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_states`, `n_features`)                if 'diag',
            (`n_states`, `n_features`, `n_features`)  if 'full'
    converged_ : bool
        True when convergence was reached in fit(), False
        otherwise.

    Methods
    -------
    decode(X)
        Find most likely mixture components for each point in `X`.
    eval(X)
        Compute a lower-bound of the log likelihood of `X` under the model
        and an approximate posterior distribution over mixture components.
    fit(X)
        Estimate the posterior of themodel parameters from `X` using the
        variational mean-field algorithm.
    predict(X)
        Like decode, find most likely mixtures components for each
        observation in `X`.
    rvs(n=1)
        Generate `n` samples from the posterior for the model.
    score(X)
        Compute the log likelihood of `X` under the model.

    See Also
    --------
    GMM : Finite gaussian mixture model fit with EM

    VBGMM : Finite gaussian mixture model fit with a variational
    algorithm, better for situations where there might be too little
    data to get a good estimate of the covariance matrix.

    """

    def __init__(self, n_states=1, cvtype='diag', alpha=1.0, rng=np.random):
        self.alpha = alpha
        super(DPGMM, self).__init__(n_states, cvtype, rng=rng)

    def _get_precisions(self):
        """Return precisions as a full matrix."""
        if self.cvtype == 'full':
            return self._covars
        elif self.cvtype == 'diag':
            return [np.diag(cov) for cov in self._covars]
        elif self.cvtype == 'tied':
            return [self._covars] * self._n_states
        elif self.cvtype == 'spherical':
            return [np.eye(self.n_features) * f for f in self._covars]

    def _set_precisions(self, covars):
        covars = np.asanyarray(covars)
        mixture._validate_covars(covars, self._cvtype, 
                                 self._n_states, self.n_features)
        self._covars = covars

    def _get_covars(self):
        return [linalg.pinv(c) for c in self._get_precisions()]

    def _set_covars(self, covars):
        self._covars = [linalg.pinv(c) for c in covars]

    precisions = property(_get_precisions, _set_covars)
    covars = property(_get_covars, _set_covars)

    def _wishart_detlogw(self, a, b, detB):
        l = 0.
        for i in xrange(self.n_features):
            l += digamma(0.5 * (a - i + 1))
        l += self.n_features * np.log(2)
        return l + detB

    def _bound_pxgivenz(self, x, k):
        bound = -0.5 * self.n_features * np.log(2 * np.pi)
        bound -= np.log(2 * np.pi * np.e)
        if self.cvtype == 'spherical':
            bound += 0.5 * self.n_features * (digamma(self._a[k])
                                              - np.log(self._b[k]))
            bound -= 0.5 * (self._covars[k]) * (norm(x - self._means[k])
                                                + self.n_features)
        elif self.cvtype == 'diag':
            bound += 0.5 * np.sum(digamma(self._a[k]) - np.log(self._b[k]))
            bound -= 0.5 * diagnorm(x - self._means[k], self._covars[k])
            bound -= 0.5 * np.sum(self._covars[k])
        elif self.cvtype == 'tied' or self.cvtype == 'full':
            if self.cvtype == 'tied':
                a, B, detB, c = self._a, self._B, self._detB, self._covars
            else:
                a, B, detB, c = (self._a[k], self._B[k],
                                 self._detB[k], self._covars[k])
            bound += 0.5 * self._wishart_detlogw(a, B, detB)
            bound -= 0.5 * squarenorm(x - self._means[k], c)
            bound -= 0.5 * a * np.trace(B)
        else:
            raise NotImplementedError("This ctype is not implemented: "
                                      + self.cvtype)
        return bound

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
        posteriors: array_like, shape (n_samples, n_states)
            Posterior probabilities of each mixture component for each
            observation
        """
        if obs is None:
            z = self._z
            obs = self._X
        else:
            z = np.zeros((obs.shape[0], self.n_states))
        obs = np.asanyarray(obs)
        p = np.zeros(self.n_states)
        bound = np.zeros(obs.shape[0])
        for i in xrange(obs.shape[0]):
            for k in xrange(self.n_states):
                p[k] = z[i, k] = self._bound_pxgivenz(obs[i], k)
                z[i, k] += digamma(self._gamma[k, 1])
                z[i, k] -= digamma(self._gamma[k, 1] + self._gamma[k, 2])
                for j in xrange(k):
                    z[i, k] += digamma(self._gamma[j, 2])
                    z[i, k] -= digamma(self._gamma[j, 1] + self._gamma[j, 2])
            z[i] = lognormalize(z[i])
            bound[i] = np.sum(z[i] * p)
        return bound, z

    def _update_gamma(self):
        for i in xrange(self.n_states):
            self._gamma[i, 1] = 1. + np.sum(self._z.T[i])
            self._gamma[i, 2] = self.alpha
            for k in xrange(i + 1, self.n_states):
                self._gamma[i, 2] += np.sum(self._z.T[k])

    def _update_mu(self):
        for k in xrange(self.n_states):
            if self.cvtype == 'spherical' or self.cvtype == 'diag':
                num = self._X[0] * self._z[0, k]
                for i in xrange(1, self._X.shape[0]):
                    num += self._X[i] * self._z[i, k]
                num *= self._covars[k]
                den = 1. + self._covars[k] * np.sum(self._z.T[k])
                self._means[k] = num / den
            elif self.cvtype == 'tied' or self.cvtype == 'full':
                if self.cvtype == 'tied':
                    cov = self._covars
                else:
                    cov = self._covars[k]
                den = np.identity(self.n_features) + cov * np.sum(self._z.T[k])
                num = self._X[0] * self._z[0, k]
                for i in xrange(1, self._X.shape[0]):
                    num += self._X[i] * self._z[i, k]
                num = np.dot(cov, num)
                self._means[k] = linalg.lstsq(den, num)[0]

    def _update_ab(self):
        if self.cvtype == 'spherical':
            self._a = 0.5 * self.n_features * np.sum(self._z, axis=0)
            for k in xrange(self.n_states):
                self._b[k] = 1.
                for i in xrange(self._X.shape[0]):
                    dif = norm((self._X[i] - self._means[k]))
                    dif += self.n_features
                    self._b[k] += 0.5 * self._z[i, k] * dif
            self._covars = self._a / self._b
        elif self.cvtype == 'diag':
            for k in xrange(self.n_states):
                self._a[k].fill(1. + 0.5 * np.sum(self._z.T[k], axis=0))
                for d in xrange(self.n_features):
                    self._b[k, d] = 1.
                    for i in xrange(self._X.shape[0]):
                        dif = self._X[i, d] - self._means[k, d]
                        self._b[k, d] += 0.5 * self._z[i, k] * (dif * dif + 1)
                self._b[k, d] = min(10 * self._a[k, d], self._b[k, d])
                self._covars[k] = self._a[k] / self._b[k]
        elif self.cvtype == 'tied':
            self._a = 2 + self._X.shape[0] + self.n_features
            self._B = (self._X.shape[0] + 1) * np.identity(self.n_features)
            for i in xrange(self._X.shape[0]):
                for k in xrange(self.n_states):
                    dif = self._X[i] - self._means[k]
                    self._B += self._z[i, k] * np.dot(dif.reshape((-1, 1)),
                                                      dif.reshape((1, -1)))
            self._B = linalg.pinv(self._B)
            self._covars = self._a * self._B
            self._detB = linalg.det(self._B)
        elif self.cvtype == 'full':
            for k in xrange(self.n_states):
                T = np.sum(self._z.T[k])
                self._a[k] = 2 + T + self.n_features
                self._B[k] = (T + 1) * np.identity(self.n_features)
                for i in xrange(self._X.shape[0]):
                    dif = self._X[i] - self._means[k]
                    self._B[k] += self._z[i, k] * np.dot(dif.reshape((-1, 1)),
                                                         dif.reshape((1, -1)))
                self._B[k] = linalg.pinv(self._B[k])
                self._covars[k] = self._a[k] * self._B[k]
                self._detB[k] = linalg.det(self._B[k])

    def _monitor(self, monitor, n, end=False):
        if monitor:
            print n, self.lower_bound()
            if end == True:
                print self._gamma.T[1]

    def _do_mstep(self, params, monitor=False):
        self._monitor(monitor, 1)
        self._update_gamma()
        self._monitor(monitor, 2)
        if 'm' in params:
            self._update_mu()
        self._monitor(monitor, 3)
        if 'c' in params:
            self._update_ab()
        self._monitor(monitor, 4, end=True)

    def _initialize_gamma(self):
        self._gamma = self.alpha * np.ones((self.n_states, 3))

    def _bound_gamma(self):
        logprior = 0.
        for k in xrange(self.n_states):
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

    def _bound_mu(self):
        logprior = 0.
        for k in xrange(self.n_states):
            logprior -= 0.5 * norm(self._means[k]) + 0.5 * self.n_features
        return logprior

    def _wishart_logz(self, v, s, dets):
            z = 0.
            z += 0.5 * v * self.n_features * np.log(2)
            z += (0.25 * (self.n_features * (self.n_features - 1))
                  * np.log(np.pi))
            z += 0.5 * v * np.log(dets)
            for i in xrange(self.n_features):
                z += gammaln(0.5 * (v - i + 1))
            return z

    def _bound_wishart(self, a, B, detB):
        logprior = self._wishart_logz(a, B, detB)
        logprior -= self._wishart_logz(self.n_features,
                                       np.identity(self.n_features),
                                       1)
        logprior += 0.5 * (a - 1) * self._wishart_detlogw(a, B, detB)
        logprior += 0.5 * a * np.trace(B)
        return logprior

    def _bound_ab(self):
        logprior = 0.
        if self.cvtype == 'spherical':
            for k in xrange(self.n_states):
                logprior += gammaln(self._a[k])
                logprior -= (self._a[k] - 1) * digamma(self._a[k])
                logprior += - np.log(self._b[k]) + self._a[k] - self._covars[k]
        elif self.cvtype == 'diag':
            for k in xrange(self.n_states):
                for d in xrange(self.n_features):
                    logprior += gammaln(self._a[k, d])
                    logprior -= (self._a[k, d] - 1) * digamma(self._a[k, d])
                    logprior -= np.log(self._b[k, d])
                    logprior += self._a[k, d] - self._covars[k, d]
        elif self.cvtype == 'tied':
            logprior += self._bound_wishart(self._a, self._B, self._detB)
        elif self.cvtype == 'full':
            for k in xrange(self.n_states):
                logprior += self._bound_wishart(self._a[k],
                                                self._B[k],
                                                self._detB[k])
        return logprior

    def _bound_z(self):
        logprior = 0.
        for i in xrange(self._z.shape[0]):
            for k in xrange(self.n_states):
                a = 0.
                for j in xrange(k + 1, self.n_states):
                    a += self._z[i, j]
                logprior += a * (digamma(self._gamma[k, 2])
                                 - digamma(self._gamma[k, 1]
                                           + self._gamma[k, 2]))
                logprior += self._z[i, k] * (digamma(self._gamma[k, 1])
                                             - digamma(self._gamma[k, 1] +
                                                       self._gamma[k, 2]))
                if self._z[i, k] != 0:
                    logprior -= self._z[i, k] * np.log(self._z[i, k])
        return logprior

    def _logprior(self):
        logprior = self._bound_gamma()
        logprior += self._bound_mu()
        logprior += self._bound_ab()
        logprior += self._bound_z()
        return logprior

    def lower_bound(self):
        c = 0.
        for i in xrange(self._X.shape[0]):
            for k in xrange(self.n_states):
                c += self._z[i, k] * self._bound_pxgivenz(self._X[i], k)
        return c + self._logprior()

    def fit(self, X, n_iter=30, thresh=1e-2, params='wmc',
            init_params='wmc', monitor=False, min_covar=None):
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
            Number of EM iterations to perform.
        thresh : float, optional
            Convergence threshold.
        params : string, optional
            Controls which parameters are updated in the training
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
        init_params : string, optional
            Controls which parameters are updated in the initialization
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
        monitor: boolean
            Prints the lower bound at every step, to help monitor
            convergence.
        """

        ## initialization step

        self._X = np.asanyarray(X)
        if hasattr(self, 'n_features') and self.n_features != self._X.shape[1]:
            raise ValueError('Unexpected number of dimensions, got %s but '
                             'expected %s' % (self._X.shape[1],
                                              self.n_features))

        self.n_features = self._X.shape[1]
        self._z = np.ones((self._X.shape[0], self.n_states))
        self._z /= self.n_states

        if init_params != '':
            self._initialize_gamma()

        if 'm' in init_params or not hasattr(self, 'means'):
            self._means = cluster.KMeans(
                k=self._n_states,rng=self.rng).fit(X).cluster_centers_[::-1]

        if 'w' in init_params or not hasattr(self, 'weights'):
            self.weights = np.tile(1.0 / self._n_states, self._n_states)

        if 'c' in init_params or not hasattr(self, 'covars'):
            if self.cvtype == 'spherical':
                self._a = np.ones(self.n_states)
                self._b = np.ones(self.n_states)
                self._covars = np.ones(self.n_states)
            elif self.cvtype == 'diag':
                self._a = 1 + 0.5 * self.n_features
                self._a *= np.ones((self.n_states, self.n_features))
                self._b = np.ones((self.n_states, self.n_features))
                self._covars = np.ones((self.n_states, self.n_features))
            elif self.cvtype == 'tied':
                self._a = 1.
                self._B = np.identity(self.n_features)
                self._covars = np.identity(self.n_features)
                self._detB = 1.
            elif self.cvtype == 'full':
                self._a = (1 + self.n_states + self._X.shape[0])
                self._a *= np.ones(self.n_states)
                self._B = [2 * np.identity(self.n_features)
                           for i in xrange(self.n_states)]
                self._covars = [np.identity(self.n_features)
                                for i in xrange(self.n_states)]
                self._detB = np.ones(self.n_states)

        logprob = []
        # reset self.converged_ to False
        self.converged_ = False
        for i in xrange(n_iter):
            # Expectation step
            curr_logprob, z = self.eval()
            logprob.append(curr_logprob.sum() + self._logprior())

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < thresh:
                self.converged_ = True
                break

            # Maximization step
            self._do_mstep(params, monitor)

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
    n_states : int, optional
        Number of mixture components. Defaults to 1.

    cvtype : string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    alpha : float, optional
        Real number representing the concentration parameter of
        the dirichlet distribution. Intuitively, the higher the
        value of alpha the more likely the variational mixture of
        gaussians model will use all components it can. Defaults 
        to 1.


    Attributes
    ----------
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the DP-GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    n_features : int
        Dimensionality of the Gaussians.
    n_states : int (read-only)
        Number of mixture components.
    weights : array, shape (`n_states`,)
        Mixing weights for each mixture component.
    means : array, shape (`n_states`, `n_features`)
        Mean parameters for each mixture component.
    precisions : array
        Precision (inverse covariance) parameters for each mixture
        component.  The shape depends on `cvtype`:
            (`n_states`,)                             if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_states`, `n_features`)                if 'diag',
            (`n_states`, `n_features`, `n_features`)  if 'full'
    converged_ : bool
        True when convergence was reached in fit(), False
        otherwise.

    Methods
    -------
    decode(X)
        Find most likely mixture components for each point in `X`.
    eval(X)
        Compute a lower-bound of the log likelihood of `X` under the model
        and an approximate posterior distribution over mixture components.
    fit(X)
        Estimate the posterior of themodel parameters from `X` using the
        variational mean-field algorithm.
    predict(X)
        Like decode, find most likely mixtures components for each
        observation in `X`.
    rvs(n=1)
        Generate `n` samples from the posterior for the model.
    score(X)
        Compute the log likelihood of `X` under the model.


    See Also
    --------
    GMM : Finite gaussian mixture model fit with EM

    DPGMM : Ininite gaussian mixture model, using the dirichlet
    process, fit with a variational algorithm
    """

    def __init__(self, n_states=1, cvtype='diag', alpha=1.0):
        super(VBGMM, self).__init__(n_states, cvtype)
        self.alpha = float(alpha) / n_states

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
        posteriors: array_like, shape (n_samples, n_states)
            Posterior probabilities of each mixture component for each
            observation
        """
        if obs is None:
            z = self._z
            obs = self._X
        else:
            z = np.zeros((obs.shape[0], self.n_states))
        obs = np.asanyarray(obs)
        p = np.zeros(self.n_states)
        bound = np.zeros(obs.shape[0])
        for i in xrange(obs.shape[0]):
            for k in xrange(self.n_states):
                p[k] = z[i, k] = self._bound_pxgivenz(obs[i], k)
                z[i, k] += digamma(self._gamma[k])
                z[i, k] -= digamma(np.sum(self._gamma))
            z[i] = lognormalize(z[i])
            bound[i] = np.sum(z[i] * p)
        return bound, z

    def _update_gamma(self):
        for i in xrange(self.n_states):
            self._gamma[i] = self.alpha + np.sum(self._z.T[i])

    def _initialize_gamma(self):
        self._gamma = self.alpha * np.ones(self.n_states)

    def _bound_z(self):
        logprior = 0.
        for i in xrange(self._z.shape[0]):
            for k in xrange(self.n_states):
                logprior += self._z[i, k] * (digamma(self._gamma[k])
                                             - digamma(np.sum(self._gamma)))
                logprior -= self._z[i, k] * np.log(self._z[i, k])
        return logprior

    def _bound_gamma(self):
        logprior = 0.
        logprior = gammaln(np.sum(self._gamma)) - gammaln(self.n_states
                                                          * self.alpha)
        for k in xrange(self.n_states):
            logprior -= gammaln(self._gamma[k]) - gammaln(self.alpha)
            a = (self._gamma[k] - self.alpha) * (digamma(self._gamma[k]) -
                                                 digamma(np.sum(self._gamma)))
            logprior += a
        return logprior

    def _monitor(self, monitor, n, end=False):
        if monitor:
            print n, self.lower_bound()
            if end == True:
                print self._gamma
