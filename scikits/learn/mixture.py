"""
Gaussian Mixture Models
"""

# Author: Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#

import numpy as np

from .base import BaseEstimator
from . import cluster


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
        Asum[np.isnan(Asum)] = - np.Inf
    return Asum


# TODO: this lacks a docstring
def normalize(A, axis=None):
    A += np.finfo(float).eps
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


def sample_gaussian(mean, covar, cvtype='diag', n_samples=1):
    """Generate random samples from a Gaussian distribution.

    Parameters
    ----------
    mean : array_like, shape (n_features,)
        Mean of the distribution.

    covars : array_like, optional
        Covariance of the distribution.  The shape depends on `cvtype`:
            scalar  if 'spherical',
            (D)     if 'diag',
            (D, D)  if 'tied', or 'full'

    cvtype : string, optional
        Type of the covariance parameters.  Must be one of
        'spherical', 'tied', 'diag', 'full'.  Defaults to 'diag'.

    n_samples : int, optional
        Number of samples to generate. Defaults to 1.

    Returns
    -------
    obs : array, shape (n_features, n_samples)
        Randomly generated sample
    """
    n_dim = len(mean)
    rand = np.random.randn(n_dim, n_samples)
    if n_samples == 1:
        rand.shape = (n_dim,)

    if cvtype == 'spherical':
        rand *= np.sqrt(covar)
    elif cvtype == 'diag':
        rand = np.dot(np.diag(np.sqrt(covar)), rand)
    else:
        from scipy import linalg
        U, s, V = linalg.svd(covar)
        sqrtS = np.diag(np.sqrt(s))
        sqrt_covar = np.dot(U, np.dot(sqrtS, V))
        rand = np.dot(sqrt_covar, rand)

    return (rand.T + mean).T


class GMM(BaseEstimator):
    """Gaussian Mixture Model

    Representation of a Gaussian mixture model probability distribution.
    This class allows for easy evaluation of, sampling from, and
    maximum-likelihood estimation of the parameters of a GMM distribution.

    Initializes parameters such that every mixture component has zero
    mean and identity covariance.


    Parameters
    ----------
    n_states : int, optional
        Number of mixture components. Defaults to 1.

    cvtype : string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.


    Attributes
    ----------
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    n_features : int
        Dimensionality of the Gaussians.
    n_states : int (read-only)
        Number of mixture components.
    weights : array, shape (`n_states`,)
        Mixing weights for each mixture component.
    means : array, shape (`n_states`, `n_features`)
        Mean parameters for each mixture component.
    covars : array
        Covariance parameters for each mixture component.  The shape
        depends on `cvtype`:
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
        Compute the log likelihood of `X` under the model and the
        posterior distribution over mixture components.
    fit(X)
        Estimate model parameters from `X` using the EM algorithm.
    predict(X)
        Like decode, find most likely mixtures components for each
        observation in `X`.
    rvs(n=1)
        Generate `n` samples from the model.
    score(X)
        Compute the log likelihood of `X` under the model.

    Examples
    --------
    >>> import numpy as np
    >>> from scikits.learn import mixture
    >>> g = mixture.GMM(n_states=2)

    >>> # Generate random observations with two modes centered on 0
    >>> # and 10 to use for training.
    >>> np.random.seed(0)
    >>> obs = np.concatenate((np.random.randn(100, 1),
    ...                       10 + np.random.randn(300, 1)))
    >>> g.fit(obs)
    GMM(cvtype='diag', n_states=2)
    >>> g.weights
    array([ 0.25,  0.75])
    >>> g.means
    array([[ 0.05980802],
           [ 9.94199467]])
    >>> g.covars
    [array([[ 1.01682662]]), array([[ 0.96080513]])]
    >>> np.round(g.weights, 2)
    array([ 0.25,  0.75])
    >>> np.round(g.means, 2)
    array([[ 0.06],
           [ 9.94]])
    >>> np.round(g.covars, 2)
    ... #doctest: +NORMALIZE_WHITESPACE
    array([[[ 1.02]],
           [[ 0.96]]])
    >>> g.predict([[0], [2], [9], [10]])
    array([0, 0, 1, 1])
    >>> np.round(g.score([[0], [2], [9], [10]]), 2)
    array([-2.32, -4.16, -1.65, -1.19])

    >>> # Refit the model on new data (initial parameters remain the
    >>> # same), this time with an even split between the two modes.
    >>> g.fit(20 * [[0]] +  20 * [[10]])
    GMM(cvtype='diag', n_states=2)
    >>> np.round(g.weights, 2)
    array([ 0.5,  0.5])
    """

    def __init__(self, n_states=1, cvtype='diag'):
        self._n_states = n_states
        self._cvtype = cvtype

        if not cvtype in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError('bad cvtype')

        self.weights = np.ones(self._n_states) / self._n_states

        # flag to indicate exit status of fit() method: converged (True) or
        # n_iter reached (False)
        self.converged_ = False

    # Read-only properties.
    @property
    def cvtype(self):
        """Covariance type of the model.

        Must be one of 'spherical', 'tied', 'diag', 'full'.
        """
        return self._cvtype

    @property
    def n_states(self):
        """Number of mixture components in the model."""
        return self._n_states

    def _get_covars(self):
        """Return covars as a full matrix."""
        if self.cvtype == 'full':
            return self._covars
        elif self.cvtype == 'diag':
            return [np.diag(cov) for cov in self._covars]
        elif self.cvtype == 'tied':
            return [self._covars] * self._n_states
        elif self.cvtype == 'spherical':
            return [np.eye(self.n_features) * f for f in self._covars]

    def _set_covars(self, covars):
        covars = np.asanyarray(covars)
        _validate_covars(covars, self._cvtype, self._n_states, self.n_features)
        self._covars = covars

    covars = property(_get_covars, _set_covars)

    def _get_means(self):
        """Mean parameters for each mixture component."""
        return self._means

    def _set_means(self, means):
        means = np.asarray(means)
        if hasattr(self, 'n_features') and \
               means.shape != (self._n_states, self.n_features):
            raise ValueError('means must have shape (n_states, n_features)')
        self._means = means.copy()
        self.n_features = self._means.shape[1]

    means = property(_get_means, _set_means)

    def _get_weights(self):
        """Mixing weights for each mixture component."""
        return np.exp(self._log_weights)

    def _set_weights(self, weights):
        if len(weights) != self._n_states:
            raise ValueError('weights must have length n_states')
        if not np.allclose(np.sum(weights), 1.0):
            raise ValueError('weights must sum to 1.0')

        self._log_weights = np.log(np.asarray(weights).copy())

    weights = property(_get_weights, _set_weights)

    def eval(self, obs):
        """Evaluate the model on data

        Compute the log probability of `obs` under the model and
        return the posterior distribution (responsibilities) of each
        mixture component for each element of `obs`.

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
        obs = np.asanyarray(obs)
        lpr = (lmvnpdf(obs, self._means, self._covars, self._cvtype)
               + self._log_weights)
        logprob = logsum(lpr, axis=1)
        posteriors = np.exp(lpr - logprob[:, np.newaxis])
        return logprob, posteriors

    def score(self, obs):
        """Compute the log probability under the model.

        Parameters
        ----------
        obs : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in `obs`
        """
        logprob, posteriors = self.eval(obs)
        return logprob

    def decode(self, obs):
        """Find most likely mixture components for each point in `obs`.

        Parameters
        ----------
        obs : array_like, shape (n, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprobs : array_like, shape (n_samples,)
            Log probability of each point in `obs` under the model.
        components : array_like, shape (n_samples,)
            Index of the most likelihod mixture components for each observation
        """
        logprob, posteriors = self.eval(obs)
        return logprob, posteriors.argmax(axis=1)

    def predict(self, X):
        """Predict label for data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = (n_samples,)
        """
        logprob, components = self.decode(X)
        return components

    def predict_proba(self, X):
        """Predict posterior probability of data under each Gaussian
        in the model.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        T : array-like, shape = (n_samples, n_states)
            Returns the probability of the sample for each Gaussian
            (state) in the model.
        """
        logprob, posteriors = self.eval(X)
        return posteriors

    def rvs(self, n_samples=1):
        """Generate random samples from the model.

        Parameters
        ----------
        n_samples : int, optional
            Number of samples to generate. Defaults to 1.

        Returns
        -------
        obs : array_like, shape (n_samples, n_features)
            List of samples
        """
        weight_pdf = self.weights
        weight_cdf = np.cumsum(weight_pdf)

        obs = np.empty((n_samples, self.n_features))
        rand = np.random.rand(n_samples)
        # decide which component to use for each sample
        comps = weight_cdf.searchsorted(rand)
        # for each component, generate all needed samples
        for comp in xrange(self._n_states):
            # occurrences of current component in obs
            comp_in_obs = (comp == comps)
            # number of those occurrences
            num_comp_in_obs = comp_in_obs.sum()
            if num_comp_in_obs > 0:
                if self._cvtype == 'tied':
                    cv = self._covars
                else:
                    cv = self._covars[comp]
                obs[comp_in_obs] = sample_gaussian(
                    self._means[comp], cv, self._cvtype, num_comp_in_obs).T
        return obs

    def fit(self, X, n_iter=10, min_covar=1e-3, thresh=1e-2, params='wmc',
            init_params='wmc'):
        """Estimate model parameters with the expectation-maximization
        algorithm.

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
        min_covar : float, optional
            Floor on the diagonal of the covariance matrix to prevent
            overfitting.  Defaults to 1e-3.
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
        """

        ## initialization step

        X = np.asanyarray(X)

        if hasattr(self, 'n_features') and self.n_features != X.shape[1]:
            raise ValueError('Unexpected number of dimensions, got %s but '
                             'expected %s' % (X.shape[1], self.n_features))

        self.n_features = X.shape[1]

        if 'm' in init_params:
            self._means = cluster.KMeans(
                k=self._n_states).fit(X).cluster_centers_
        elif not hasattr(self, 'means'):
                self._means = np.zeros((self.n_states, self.n_features))

        if 'w' in init_params or not hasattr(self, 'weights'):
            self.weights = np.tile(1.0 / self._n_states, self._n_states)

        if 'c' in init_params:
            cv = np.cov(X.T)
            if not cv.shape:
                cv.shape = (1, 1)
            self._covars = _distribute_covar_matrix_to_match_cvtype(
                cv, self._cvtype, self._n_states)
        elif not hasattr(self, 'covars'):
                self.covars = _distribute_covar_matrix_to_match_cvtype(
                    np.eye(self.n_features), self.cvtype, self.n_states)

        # EM algorithm
        logprob = []
        # reset self.converged_ to False
        self.converged_ = False
        for i in xrange(n_iter):
            # Expectation step
            curr_logprob, posteriors = self.eval(X)
            logprob.append(curr_logprob.sum())

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < thresh:
                self.converged_ = True
                break

            # Maximization step
            self._do_mstep(X, posteriors, params, min_covar)

        return self

    def _do_mstep(self, X, posteriors, params, min_covar=0):
            w = posteriors.sum(axis=0)
            avg_obs = np.dot(posteriors.T, X)
            norm = 1.0 / (w[:, np.newaxis] + 1e-200)

            if 'w' in params:
                self._log_weights = np.log(w / w.sum())
            if 'm' in params:
                self._means = avg_obs * norm
            if 'c' in params:
                covar_mstep_func = _covar_mstep_funcs[self._cvtype]
                self._covars = covar_mstep_func(self, X, posteriors,
                                                avg_obs, norm, min_covar)

            return w


##
## some helper routines
##


def _lmvnpdfdiag(obs, means=0.0, covars=1.0):
    n_obs, n_dim = obs.shape
    # (x-y).T A (x-y) = x.T A x - 2x.T A y + y.T A y
    #lpr = -0.5 * (np.tile((np.sum((means**2) / covars, 1)
    #                  + np.sum(np.log(covars), 1))[np.newaxis,:], (n_obs,1))
    lpr = -0.5 * (n_dim * np.log(2 * np.pi) + np.sum(np.log(covars), 1)
                  + np.sum((means ** 2) / covars, 1)
                  - 2 * np.dot(obs, (means / covars).T)
                  + np.dot(obs ** 2, (1.0 / covars).T))
    return lpr


def _lmvnpdfspherical(obs, means=0.0, covars=1.0):
    cv = covars.copy()
    if covars.ndim == 1:
        cv = cv[:, np.newaxis]
    return _lmvnpdfdiag(obs, means, np.tile(cv, (1, obs.shape[-1])))


def _lmvnpdftied(obs, means, covars):
    from scipy import linalg
    n_obs, n_dim = obs.shape
    # (x-y).T A (x-y) = x.T A x - 2x.T A y + y.T A y
    icv = linalg.pinv(covars)
    lpr = -0.5 * (n_dim * np.log(2 * np.pi) + np.log(linalg.det(covars))
                  + np.sum(obs * np.dot(obs, icv), 1)[:, np.newaxis]
                  - 2 * np.dot(np.dot(obs, icv), means.T)
                  + np.sum(means * np.dot(means, icv), 1))
    return lpr


def _lmvnpdffull(obs, means, covars):
    """
    Log probability for full covariance matrices.
    """
    from scipy import linalg
    import itertools
    if hasattr(linalg, 'solve_triangular'):
        # only in scipy since 0.9
        solve_triangular = linalg.solve_triangular
    else:
        # slower, but works
        solve_triangular = linalg.solve
    n_obs, n_dim = obs.shape
    nmix = len(means)
    log_prob = np.empty((n_obs, nmix))
    for c, (mu, cv) in enumerate(itertools.izip(means, covars)):
        cv_chol = linalg.cholesky(cv, lower=True)
        cv_log_det = 2 * np.sum(np.log(np.diagonal(cv_chol)))
        cv_sol = solve_triangular(cv_chol, (obs - mu).T, lower=True).T
        log_prob[:, c] = -.5 * (np.sum(cv_sol ** 2, axis=1) + \
                           n_dim * np.log(2 * np.pi) + cv_log_det)

    return log_prob


def _validate_covars(covars, cvtype, nmix, n_dim):
    from scipy import linalg
    if cvtype == 'spherical':
        if len(covars) != nmix:
            raise ValueError("'spherical' covars must have length nmix")
        elif np.any(covars <= 0):
            raise ValueError("'spherical' covars must be non-negative")
    elif cvtype == 'tied':
        if covars.shape != (n_dim, n_dim):
            raise ValueError("'tied' covars must have shape (n_dim, n_dim)")
        elif (not np.allclose(covars, covars.T)
              or np.any(linalg.eigvalsh(covars) <= 0)):
            raise ValueError("'tied' covars must be symmetric, "
                             "positive-definite")
    elif cvtype == 'diag':
        if covars.shape != (nmix, n_dim):
            raise ValueError("'diag' covars must have shape (nmix, n_dim)")
        elif np.any(covars <= 0):
            raise ValueError("'diag' covars must be non-negative")
    elif cvtype == 'full':
        if covars.shape != (nmix, n_dim, n_dim):
            raise ValueError("'full' covars must have shape "
                             "(nmix, n_dim, n_dim)")
        for n, cv in enumerate(covars):
            if (not np.allclose(cv, cv.T)
                or np.any(linalg.eigvalsh(cv) <= 0)):
                raise ValueError("component %d of 'full' covars must be "
                                 "symmetric, positive-definite" % n)


def _distribute_covar_matrix_to_match_cvtype(tiedcv, cvtype, n_states):
    if cvtype == 'spherical':
        cv = np.tile(np.diag(tiedcv).mean(), n_states)
    elif cvtype == 'tied':
        cv = tiedcv
    elif cvtype == 'diag':
        cv = np.tile(np.diag(tiedcv), (n_states, 1))
    elif cvtype == 'full':
        cv = np.tile(tiedcv, (n_states, 1, 1))
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
    avg_means2 = gmm._means ** 2
    avg_obs_means = gmm._means * avg_obs * norm
    return avg_obs2 - 2 * avg_obs_means + avg_means2 + min_covar


def _covar_mstep_spherical(*args):
    return _covar_mstep_diag(*args).mean(axis=1)


def _covar_mstep_full(gmm, obs, posteriors, avg_obs, norm, min_covar):
    # Eq. 12 from K. Murphy, "Fitting a Conditional Linear Gaussian
    # Distribution"
    cv = np.empty((gmm._n_states, gmm.n_features, gmm.n_features))
    for c in xrange(gmm._n_states):
        post = posteriors[:, c]
        avg_cv = np.dot(post * obs.T, obs) / post.sum()
        mu = gmm._means[c][np.newaxis]
        cv[c] = (avg_cv - np.dot(mu.T, mu)
                 + min_covar * np.eye(gmm.n_features))
    return cv


def _covar_mstep_tied2(*args):
    return _covar_mstep_full(*args).mean(axis=0)


def _covar_mstep_tied(gmm, obs, posteriors, avg_obs, norm, min_covar):
    print "THIS IS BROKEN"
    # Eq. 15 from K. Murphy, "Fitting a Conditional Linear Gaussian
    avg_obs2 = np.dot(obs.T, obs)
    avg_means2 = np.dot(gmm._means.T, gmm._means)
    return (avg_obs2 - avg_means2 + min_covar * np.eye(gmm.n_features))


def _covar_mstep_slow(gmm, obs, posteriors, avg_obs, norm, min_covar):
    w = posteriors.sum(axis=0)
    covars = np.zeros(gmm._covars.shape)
    for c in xrange(gmm._n_states):
        mu = gmm._means[c]
        #cv = np.dot(mu.T, mu)
        avg_obs2 = np.zeros((gmm.n_features, gmm.n_features))
        for t, o in enumerate(obs):
            avg_obs2 += posteriors[t, c] * np.outer(o, o)
        cv = (avg_obs2 / w[c]
              - 2 * np.outer(avg_obs[c] / w[c], mu)
              + np.outer(mu, mu)
              + min_covar * np.eye(gmm.n_features))
        if gmm.cvtype == 'spherical':
            covars[c] = np.diag(cv).mean()
        elif gmm.cvtype == 'diag':
            covars[c] = np.diag(cv)
        elif gmm.cvtype == 'full':
            covars[c] = cv
        elif gmm.cvtype == 'tied':
            covars += cv / gmm._n_states
    return covars


_covar_mstep_funcs = {'spherical': _covar_mstep_spherical,
                      'diag': _covar_mstep_diag,
                      #'tied': _covar_mstep_tied,
                      'full': _covar_mstep_full,
                      'tied': _covar_mstep_slow,
                      }
