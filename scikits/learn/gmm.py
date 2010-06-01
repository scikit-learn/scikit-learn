import itertools
import time
import logging

import numpy as np
from scipy import cluster

ZEROLOGPROB = -1e200

log = logging.getLogger('gmm.gmm')

def almost_equal(actual, desired, decimal=7):
    """Check that two floats are approximately equal."""
    return abs(desired - actual) < 0.5 * 10**(-decimal)

def logsum(A, axis=None):
    """Computes the sum of A assuming A is in the log domain.

    Returns log(sum(exp(A), axis)) while minimizing the possibility of
    over/underflow.
    """
    Amax = A.max(axis)
    if axis and A.ndim > 1:
        shape = list(A.shape)
        shape[axis] = 1
        Amax.shape = shape
    Asum = np.log(np.sum(np.exp(A - Amax), axis))
    Asum += Amax.reshape(Asum.shape)
    if axis:
        # Look out for underflow.
        Asum[np.isnan(Asum)] = -np.Inf
    return Asum

def normalize(A, axis=None):
    Asum = A.sum(axis)
    if axis and A.ndim > 1:
        # Make sure we don't divide by zero.
        Asum[Asum == 0] = 1
        shape = list(A.shape)
        shape[axis] = 1
        Asum.shape = shape
    return A / Asum

def lmvnpdf(obs, means, covars, cvtype='diag'):
    """Compute the log probability under a multivariate Gaussian distribution.

    Parameters
    ----------
    obs : array_like, shape (O, D)
        List of D-dimensional data points.  Each row corresponds to a
        single data point.
    means : array_like, shape (C, D)
        List of D-dimensional mean vectors for C Gaussians.  Each row
        corresponds to a single mean vector.
    covars : array_like
        List of C covariance parameters for each Gaussian.  The shape
        depends on `cvtype`:
            (C,)      if 'spherical',
            (D, D)    if 'tied',
            (C, D)    if 'diag',
            (C, D, D) if 'full'
    cvtype : string
        Type of the covariance parameters.  Must be one of
        'spherical', 'tied', 'diag', 'full'.  Defaults to 'diag'.

    Returns
    -------
    lpr : array_like, shape (O, C)
        Array containing the log probabilities of each data point in
        `obs` under each of the C multivariate Gaussian distributions.
    """
    lmvnpdf_dict = {'spherical': _lmvnpdfspherical,
                    'tied': _lmvnpdftied,
                    'diag': _lmvnpdfdiag,
                    'full': _lmvnpdffull}
    return lmvnpdf_dict[cvtype](obs, means, covars)


def sample_gaussian(mean, covar, cvtype='diag', n=1):
    """Generate random samples from a Gaussian distribution.

    Parameters
    ----------
    mean : array_like, shape (ndim,)
        Mean of the distribution.
    covars : array_like
        Covariance of the distribution.  The shape depends on `cvtype`:
            scalar  if 'spherical',
            (D)     if 'diag',
            (D, D)  if 'tied', or 'full'
    cvtype : string
        Type of the covariance parameters.  Must be one of
        'spherical', 'tied', 'diag', 'full'.  Defaults to 'diag'.
    n : int
        Number of samples to generate.

    Returns
    -------
    obs : array, shape (n, ndim)
        Randomly generated sample
    """
    ndim = len(mean)
    rand = np.random.randn(ndim, n)
    if n == 1:
        rand.shape = (ndim,)

    if cvtype == 'spherical':
        rand *= np.sqrt(covar)
    elif cvtype == 'diag':
        rand = np.dot(np.diag(np.sqrt(covar)), rand)
    else:
        U, s, V = np.linalg.svd(covar)
        sqrtS = np.diag(np.sqrt(s))
        sqrt_covar = np.dot(U, np.dot(sqrtS, V))
        rand = np.dot(sqrt_covar, rand)

    return (rand.T + mean).T


class GMM(object):
    """Gaussian Mixture Model

    Representation of a Gaussian mixture model probability distribution.
    This class allows for easy evaluation of, sampling from, and
    maximum-likelihood estimation of the parameters of a GMM distribution.

    Attributes
    ----------
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    ndim : int (read-only)
        Dimensionality of the Gaussians.
    nstates : int (read-only)
        Number of states (mixture components).
    weights : array, shape (`nstates`,)
        Mixing weights for each mixture component.
    means : array, shape (`nstates`, `ndim`)
        Mean parameters for each mixture component.
    covars : array
        Covariance parameters for each mixture component.  The shape
        depends on `cvtype`:
            (`nstates`,)                if 'spherical',
            (`ndim`, `ndim`)            if 'tied',
            (`nstates`, `ndim`)         if 'diag',
            (`nstates`, `ndim`, `ndim`) if 'full'
    labels : list, len `nstates`
        Optional labels for each mixture component.

    Methods
    -------
    eval(obs)
        Compute the log likelihood of `obs` under the model.
    decode(obs)
        Find most likely mixture components for each point in `obs`.
    rvs(n=1)
        Generate `n` samples from the model.
    init(obs)
        Initialize model parameters from `obs`.
    fit(obs)
        Estimate model parameters from `obs` using the EM algorithm.

    Examples
    --------
    >>> gmm = GMM(2, ndim=1)
    >>> obs = numpy.concatenate((numpy.random.randn(100, 1),
    ...                          10 + numpy.random.randn(300, 1)))
    >>> # Roughly initialize the model parameters.
    >>> gmm.init(obs)
    >>> gmm.fit(obs)
    >>> gmm.weights, gmm.means, gmm.covars
    (array([ 0.25,  0.75]),
     array([[ -0.22744484],
           [ 10.07096441]]),
     array([[ 1.02857617],
           [ 1.11389491]]))
    >>> gmm.decode([0, 2, 9, 10])
    array([0, 0, 1, 1])
    >>> # Refit the model on new data (initial parameters remain the same).
    >>> gmm.fit(numpy.concatenate((20 * [0], 20 * [10])))
    """

    def __init__(self, nstates=1, ndim=1, cvtype='diag'):
        """Create a Gaussian mixture model

        Initializes parameters such that every mixture component has
        zero mean and identity covariance.

        Parameters
        ----------
        ndim : int
            Dimensionality of the mixture components.
        nstates : int
            Number of mixture components.
        cvtype : string (read-only)
            String describing the type of covariance parameters to
            use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
            Defaults to 'diag'.
        """

        self._nstates = nstates
        self.ndim = ndim
        self._cvtype = cvtype

        self.weights = np.tile(1.0 / nstates, nstates)
        self.means = np.zeros((nstates, ndim))
        self.covars = _distribute_covar_matrix_to_match_cvtype(
            np.eye(ndim), cvtype, nstates)
        
        self.labels = [None] * nstates

    # Read-only properties.
    @property
    def cvtype(self):
        """Covariance type of the model.

        Must be one of 'spherical', 'tied', 'diag', 'full'.
        """
        return self._cvtype

    @property
    def nstates(self):
        """Number of mixture components in the model."""
        return self._nstates

    def _get_weights(self):
        """Mixing weights for each mixture component."""
        return np.exp(self._log_weights)

    def _set_weights(self, weights):
        if len(weights) != self._nstates:
            raise ValueError, 'weights must have length nstates'
        if not almost_equal(np.sum(weights), 1.0):
            raise ValueError, 'weights must sum to 1.0'
        
        self._log_weights = np.log(np.asarray(weights).copy())

    weights = property(_get_weights, _set_weights)
    
    def eval(self, obs):
        """Evaluate the model on data

        Compute the log probability of `obs` under the model and
        return the posterior distribution (responsibilities) of each
        mixture component for each element of `obs`

        Parameters
        ----------
        obs : array_like, shape (n, ndim)
            List of ndim-dimensional data points.  Each row corresponds to a
            single data point.

        Returns
        -------
        logprob : array_like, shape (n,)
            Log probabilities of each data point in `obs`
        posteriors: array_like, shape (n, nstates)
            Posterior probabilities of each mixture component for each
            observation
        """
        lpr = (lmvnpdf(obs, self.means, self.covars, self._cvtype)
               + self._log_weights)
        logprob = logsum(lpr, axis=1)
        posteriors = np.exp(lpr - logprob[:,np.newaxis])
        return logprob, posteriors

    def lpdf(self, obs):
        """Compute the log probability under the model.

        Parameters
        ----------
        obs : array_like, shape (n, ndim)
            List of ndim-dimensional data points.  Each row corresponds to a
            single data point.

        Returns
        -------
        logprob : array_like, shape (n,)
            Log probabilities of each data point in `obs`
        """
        logprob, posteriors = self.eval(obs)
        return logprob

    def decode(self, obs):
        """Find most likely mixture components for each point in `obs`.

        Parameters
        ----------
        obs : array_like, shape (n, ndim)
            List of ndim-dimensional data points.  Each row corresponds to a
            single data point.

        Returns
        -------
        components : array_like, shape (n,)
            Index of the most likelihod mixture components for each observation
        """
        logprob, posteriors = self.eval(obs)
        return logprob, posteriors.argmax(axis=1)
        
    def rvs(self, n=1):
        """Generate random samples from the model.

        Parameters
        ----------
        n : int
            Number of samples to generate.

        Returns
        -------
        obs : array_like, shape (n, ndim)
            List of samples
        """
        weight_pdf = self.weights
        weight_cdf = np.cumsum(weight_pdf)

        obs = np.empty((n, self.ndim))
        for x in xrange(n):
            rand = np.random.rand()
            c = (weight_cdf > rand).argmax()
            if self._cvtype == 'tied':
                cv = self.covars
            else:
                cv = self.covars[c]
            obs[x] = sample_gaussian(self.means[c], cv, self._cvtype)
        return obs

    def init(self, obs, params='wmc', **kwargs):
        """Initialize model parameters from data using the k-means algorithm

        Parameters
        ----------
        obs : array_like, shape (n, ndim)
            List of ndim-dimensional data points.  Each row corresponds to a
            single data point.
        params : string
            Controls which parameters are updated in the training
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.
        **kwargs :
            Keyword arguments to pass through to the k-means function 
            (scipy.cluster.vq.kmeans2)

        See Also
        --------
        scipy.cluster.vq.kmeans2
        """
        
        if 'm' in params:
            self.means, tmp = cluster.vq.kmeans2(obs, self._nstates,
                                                     **kwargs)
        if 'w' in params:
            self.weights = np.tile(1.0 / self._nstates, self._nstates)
        if 'c' in params:
            cv = np.cov(obs.T)
            if not cv.shape:
                cv.shape = (1, 1)
            self.covars = _distribute_covar_matrix_to_match_cvtype(
                cv, self._cvtype, self._nstates)

    def fit(self, obs, iter=10, min_covar=1.0, thresh=1e-2, params='wmc'):
        """Estimate model parameters with the expectation-maximization
        algorithm.

        Parameters
        ----------
        obs : array_like, shape (n, ndim)
            List of ndim-dimensional data points.  Each row corresponds to a
            single data point.
        iter : int
            Number of EM iterations to perform.
        min_covar : float
            Floor on the diagonal of the covariance matrix to prevent
            overfitting.  Defaults to 1.0.
        thresh : float
            Convergence threshold.
        params : string
            Controls which parameters are updated in the training
            process.  Can contain any combination of 'w' for weights,
            'm' for means, and 'c' for covars.  Defaults to 'wmc'.

        Returns
        -------
        logprob : list
            Log probabilities of each data point in `obs` for each iteration
        """
        covar_mstep_fun = {'spherical': _covar_mstep_spherical,
                           'diag': _covar_mstep_diag,
                           #'tied': _covar_mstep_tied,
                           #'full': _covar_mstep_full,
                           'tied': _covar_mstep_slow,
                           'full': _covar_mstep_slow,
                           }[self._cvtype]

        T = time.time()
        logprob = []
        for i in xrange(iter):
            # Expectation step
            curr_logprob, posteriors = self.eval(obs)
            logprob.append(curr_logprob.sum())

            currT = time.time()
            log.info('Iteration %d: log likelihood = %f (took %f seconds).'
                      % (i, logprob[-1], currT - T))
            T = currT


            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < thresh:
                log.info('Converged at iteration %d.' % i)
                break

            # Maximization step
            w = posteriors.sum(axis=0)
            avg_obs = np.dot(posteriors.T, obs)
            norm = 1.0 / w[:,np.newaxis]
            
            if 'w' in params:
                self.weights = w / w.sum()
            if 'm' in params:
                self.means = avg_obs * norm
            if 'c' in params:
                self.covars = covar_mstep_fun(self, obs, posteriors,
                                               avg_obs, norm, min_covar)

        return logprob


def _lmvnpdfdiag(obs, means=0.0, covars=1.0):
    nobs, ndim = obs.shape
    # (x-y).T A (x-y) = x.T A x - 2x.T A y + y.T A y
    #lpr = -0.5 * (np.tile((np.sum((means**2) / covars, 1)
    #                      + np.sum(np.log(covars), 1))[np.newaxis,:], (nobs,1))
    lpr = -0.5 * (ndim * np.log(2 * np.pi) + np.sum(np.log(covars), 1)
                  + np.sum((means**2) / covars, 1)
                  - 2 * np.dot(obs, (means / covars).T)
                  + np.dot(obs**2, (1.0 / covars).T))
    return lpr

def _lmvnpdfspherical(obs, means=0.0, covars=1.0):
    cv = covars.copy()
    if covars.ndim == 1:
        cv = cv[:,np.newaxis]
    return _lmvnpdfdiag(obs, means, np.tile(cv, (1, obs.shape[-1])))

def _lmvnpdftied(obs, means, covars):
    nobs, ndim = obs.shape
    nmix = len(means)
    # (x-y).T A (x-y) = x.T A x - 2x.T A y + y.T A y
    icv = np.linalg.inv(covars)
    lpr = -0.5 * (ndim * np.log(2 * np.pi) + np.log(np.linalg.det(covars))
                  + np.sum(obs * np.dot(obs, icv), 1)[:,np.newaxis]
                  - 2 * np.dot(np.dot(obs, icv), means.T)
                  + np.sum(means * np.dot(means, icv), 1))
    return lpr

def _lmvnpdffull(obs, means, covars):
    # FIXME: this representation of covars is going to lose for caching
    nobs, ndim = obs.shape
    nmix = len(means)
    lpr = np.empty((nobs,nmix))
    for c, (mu, cv) in enumerate(itertools.izip(means, covars)):
        icv = np.linalg.inv(cv)
        lpr[:,c] = -0.5 * (ndim * np.log(2 * np.pi) + np.log(np.linalg.det(cv)))
        for o, currobs in enumerate(obs):
            dzm = (currobs - mu)
            lpr[o,c] += -0.5 * np.dot(np.dot(dzm, icv), dzm.T)
        #dzm = (obs - mu)
        #lpr[:,c] = -0.5 * (np.dot(np.dot(dzm, np.linalg.inv(cv)), dzm.T)
        #                   + np.log(2 * np.pi) + np.linalg.det(cv)).diagonal()
    return lpr

def _validate_covars(covars, cvtype, nmix, ndim):
    if cvtype == 'spherical':
        if len(covars) != nmix:
            raise ValueError, "'spherical' covars must have length nmix"
        elif np.any(covars <= 0):
            raise ValueError, "'spherical' covars must be non-negative"
    elif cvtype == 'tied':
        if covars.shape != (ndim, ndim):
            raise ValueError, "'tied' covars must have shape (ndim, ndim)"
        elif (not np.all(almost_equal(covars, covars.T))
              or np.any(np.linalg.eigvalsh(covars) <= 0)):
            raise (ValueError,
                   "'tied' covars must be symmetric, positive-definite")
    elif cvtype == 'diag':
        if covars.shape != (nmix, ndim):
            raise ValueError, "'diag' covars must have shape (nmix, ndim)"
        elif np.any(covars <= 0):
            raise ValueError, "'diag' covars must be non-negative"
    elif cvtype == 'full':
        if covars.shape != (nmix, ndim, ndim):
            raise (ValueError,
                   "'full' covars must have shape (nmix, ndim, ndim)")
        for n,cv in enumerate(covars):
            if (not np.all(almost_equal(cv, cv.T))
                or np.any(np.linalg.eigvalsh(cv) <= 0)):
                raise (ValueError,
                       "component %d of 'full' covars must be symmetric,"
                       "positive-definite" % n)

def _distribute_covar_matrix_to_match_cvtype(tiedcv, cvtype, nstates):
    if cvtype == 'spherical':
        cv = np.tile(np.diag(tiedcv).mean(), nstates)
    elif cvtype == 'tied':
        cv = tiedcv
    elif cvtype == 'diag':
        cv = np.tile(np.diag(tiedcv), (nstates, 1))
    elif cvtype == 'full':
        cv = np.tile(tiedcv, (nstates, 1, 1))
    else:
        raise (ValueError,
               "cvtype must be one of 'spherical', 'tied', 'diag', 'full'")
    return cv

def _covar_mstep_diag(gmm, obs, posteriors, avg_obs, norm, min_covar):
    # For column vectors:
    # covars_c = average((obs(t) - means_c) (obs(t) - means_c).T,
    #                    weights_c)
    # (obs(t) - means_c) (obs(t) - means_c).T
    #     = obs(t) obs(t).T - 2 obs(t) means_c.T + means_c means_c.T
    #
    # But everything here is a row vector, so all of the
    # above needs to be transposed.
    avg_obs2 = np.dot(posteriors.T, obs * obs) * norm
    avg_means2 = gmm.means**2 
    avg_obs_means = gmm.means * avg_obs * norm
    return avg_obs2 - 2 * avg_obs_means + avg_means2 + min_covar

def _covar_mstep_spherical(*args):
    return _covar_mstep_diag(*args).mean(axis=1)

def _covar_mstep_full(gmm, obs, posteriors, avg_obs, norm, min_covar):
    print "THIS IS BROKEN"
    # Eq. 12 from K. Murphy, "Fitting a Conditional Linear Gaussian
    # Distribution"
    avg_obs2 = np.dot(obs.T, obs)
    #avg_obs2 = np.dot(obs.T, avg_obs)
    cv = np.empty((gmm._nstates, gmm.ndim, gmm.ndim))
    for c in xrange(gmm._nstates):
        wobs = obs.T * posteriors[:,c]
        avg_obs2 = np.dot(wobs, obs) / posteriors[:,c].sum()
        mu = gmm.means[c][np.newaxis]
        cv[c] = (avg_obs2 - np.dot(mu, mu.T)
                 + min_covar * np.eye(gmm.ndim))
    return cv

def _covar_mstep_tied2(*args):
    return _covar_mstep_full(*args).mean(axis=0)

def _covar_mstep_tied(gmm, obs, posteriors, avg_obs, norm, min_covar):
    print "THIS IS BROKEN"
    # Eq. 15 from K. Murphy, "Fitting a Conditional Linear Gaussian
    # Distribution"
    avg_obs2 = np.dot(obs.T, obs)
    avg_means2 = np.dot(gmm.means.T, gmm.means)
    return (avg_obs2 - avg_means2 + min_covar * np.eye(gmm.ndim))

def _covar_mstep_slow(gmm, obs, posteriors, avg_obs, norm, min_covar):
    w = posteriors.sum(axis=0)
    covars = np.zeros(gmm.covars.shape)
    for c in xrange(gmm._nstates):
        mu = gmm.means[c]
        #cv = np.dot(mu.T, mu)
        avg_obs2 = np.zeros((gmm.ndim, gmm.ndim))
        for t,o in enumerate(obs):
            avg_obs2 += posteriors[t,c] * np.outer(o, o)
        cv = (avg_obs2 / w[c]
              - 2 * np.outer(avg_obs[c] / w[c], mu)
              + np.outer(mu, mu)
              + min_covar * np.eye(gmm.ndim))
        if gmm.cvtype == 'spherical':
            covars[c] = np.diag(cv).mean()
        elif gmm.cvtype == 'diag':
            covars[c] = np.diag(cv)
        elif gmm.cvtype == 'full':
            covars[c] = cv
        elif gmm.cvtype == 'tied':
            covars += cv / gmm._nstates
    return covars
