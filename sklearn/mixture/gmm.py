"""
Gaussian Mixture Models
"""

# Author: Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#

import numpy as np

from ..base import BaseEstimator
from ..utils import check_random_state
from ..utils.extmath import logsumexp
from .. import cluster


# FIXME this lacks a proper docstring
def normalize(A, axis=None):
    """ Normalize the input array so that it sums to 1.

        WARNING: Modifies inplace the array
    """
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


def sample_gaussian(mean, covar, cvtype='diag', n_samples=1,
                    random_state=None):
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
    rng = check_random_state(random_state)
    n_dim = len(mean)
    rand = rng.randn(n_dim, n_samples)
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
    n_components : int, optional
        Number of mixture components. Defaults to 1.

    cvtype : string (read-only), optional
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        Defaults to 'diag'.

    rng : numpy.random object, optional
        Must support the full numpy random number generator API.

    min_covar : float, optional
        Floor on the diagonal of the covariance matrix to prevent
        overfitting.  Defaults to 1e-3.

    thresh : float, optional
        Convergence threshold.

    Attributes
    ----------

    n_features : int
        Dimensionality of the Gaussians.

    n_states : int (read-only)
        Number of mixture components.

    `converged_` : bool
        True when convergence was reached in fit(), False
        otherwise.

    weights : property - this string will be replaced

    means : property - this string will be replaced

    cvtype : property - this string will be replaced

    covars : property - this string will be replaced


    See Also
    --------

    DPGMM : Ininite gaussian mixture model, using the dirichlet
        process, fit with a variational algorithm


    VBGMM : Finite gaussian mixture model fit with a variational
        algorithm, better for situations where there might be too little
        data to get a good estimate of the covariance matrix.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn import mixture
    >>> np.random.seed(1)
    >>> g = mixture.GMM(n_components=2)
    >>> # Generate random observations with two modes centered on 0
    >>> # and 10 to use for training.
    >>> obs = np.concatenate((np.random.randn(100, 1),
    ...                       10 + np.random.randn(300, 1)))
    >>> g.fit(obs)
    GMM(cvtype='diag', n_components=2)
    >>> np.round(g.weights, 2)
    array([ 0.75,  0.25])
    >>> np.round(g.means, 2)
    array([[ 10.05],
           [  0.06]])
    >>> np.round(g.covars, 2) #doctest: +SKIP
    array([[[ 1.02]],
           [[ 0.96]]])
    >>> g.predict([[0], [2], [9], [10]])
    array([1, 1, 0, 0])
    >>> np.round(g.score([[0], [2], [9], [10]]), 2)
    array([-2.19, -4.58, -1.75, -1.21])
    >>> # Refit the model on new data (initial parameters remain the
    >>> # same), this time with an even split between the two modes.
    >>> g.fit(20 * [[0]] +  20 * [[10]])
    GMM(cvtype='diag', n_components=2)
    >>> np.round(g.weights, 2)
    array([ 0.5,  0.5])

    """

    def __init__(self, n_components=1, cvtype='diag', random_state=None,
                 thresh=1e-2, min_covar=1e-3):
        self.n_components = n_components
        self._cvtype = cvtype
        self.thresh = thresh
        self.min_covar = min_covar
        self.random_state = random_state

        if not cvtype in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError('bad cvtype: ' + str(cvtype))

        self.weights = np.ones(self.n_components) / self.n_components

        # flag to indicate exit status of fit() method: converged (True) or
        # n_iter reached (False)
        self.converged_ = False

    # Read-only properties.
    @property
    def cvtype(self):
        """Covariance type of the model.
        String describing the type of covariance parameters used by
        the GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
        """
        return self._cvtype

    def _get_covars(self):
        """Covariance parameters for each mixture component.
        The shape depends on `cvtype`::

            (`n_states`,)                             if 'spherical',
            (`n_features`, `n_features`)              if 'tied',
            (`n_states`, `n_features`)                if 'diag',
            (`n_states`, `n_features`, `n_features`)  if 'full'

        """
        if self.cvtype == 'full':
            return self._covars
        elif self.cvtype == 'diag':
            return [np.diag(cov) for cov in self._covars]
        elif self.cvtype == 'tied':
            return [self._covars] * self.n_components
        elif self.cvtype == 'spherical':
            return [np.eye(self.n_features) * f for f in self._covars]

    def _set_covars(self, covars):
        covars = np.asarray(covars)
        _validate_covars(covars, self._cvtype, self.n_components,
                self.n_features)
        self._covars = covars

    covars = property(_get_covars, _set_covars)

    def _get_means(self):
        """Mean parameters for each mixture component.
        array, shape ``(n_states, n_features)``.
        """
        return self._means

    def _set_means(self, means):
        means = np.asarray(means)
        if hasattr(self, 'n_features') and \
               means.shape != (self.n_components, self.n_features):
            raise ValueError('means must have shape ' +
                    '(n_components, n_features)')
        self._means = means.copy()
        self.n_features = self._means.shape[1]

    means = property(_get_means, _set_means)

    def __repr__(self):
        return "GMM(cvtype='%s', n_components=%s)" % (self._cvtype,
                self.n_components)

    def _get_weights(self):
        """Mixing weights for each mixture component.
        array, shape ``(n_states,)``
        """
        return np.exp(self._log_weights)

    def _set_weights(self, weights):
        if len(weights) != self.n_components:
            raise ValueError('weights must have length n_components')
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
        obs: array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob: array_like, shape (n_samples,)
            Log probabilities of each data point in `obs`

        posteriors: array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        obs = np.asarray(obs)
        lpr = (lmvnpdf(obs, self._means, self._covars, self._cvtype)
               + self._log_weights)
        logprob = logsumexp(lpr, axis=1)
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
        logprob, _ = self.eval(obs)
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
        T : array-like, shape = (n_samples, n_components)
            Returns the probability of the sample for each Gaussian
            (state) in the model.
        """
        logprob, posteriors = self.eval(X)
        return posteriors

    def rvs(self, n_samples=1, random_state=None):
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
        if random_state is None:
            random_state = self.random_state
        random_state = check_random_state(random_state)
        weight_pdf = self.weights
        weight_cdf = np.cumsum(weight_pdf)

        obs = np.empty((n_samples, self.n_features))
        rand = random_state.rand(n_samples)
        # decide which component to use for each sample
        comps = weight_cdf.searchsorted(rand)
        # for each component, generate all needed samples
        for comp in xrange(self.n_components):
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
                    self._means[comp], cv, self._cvtype, num_comp_in_obs,
                    random_state=random_state
                ).T
        return obs

    def fit(self, X, n_iter=10, thresh=1e-2, params='wmc',
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

        X = np.asarray(X)

        if hasattr(self, 'n_features') and self.n_features != X.shape[1]:
            raise ValueError('Unexpected number of dimensions, got %s but '
                             'expected %s' % (X.shape[1], self.n_features))

        self.n_features = X.shape[1]

        if 'm' in init_params:
            self._means = cluster.KMeans(
                k=self.n_components).fit(X).cluster_centers_
        elif not hasattr(self, 'means'):
                self._means = np.zeros((self.n_components, self.n_features))

        if 'w' in init_params or not hasattr(self, 'weights'):
            self.weights = np.tile(1.0 / self.n_components, self.n_components)

        if 'c' in init_params:
            cv = np.cov(X.T)
            if not cv.shape:
                cv.shape = (1, 1)
            self._covars = _distribute_covar_matrix_to_match_cvtype(
                cv, self._cvtype, self.n_components)
        elif not hasattr(self, 'covars'):
                self.covars = _distribute_covar_matrix_to_match_cvtype(
                    np.eye(self.n_features), self.cvtype, self.n_components)

        # EM algorithm
        logprob = []
        # reset self.converged_ to False
        self.converged_ = False
        for i in xrange(n_iter):
            # Expectation step
            curr_logprob, posteriors = self.eval(X)
            logprob.append(curr_logprob.sum())

            # Check for convergence.
            if i > 0 and abs(logprob[-1] - logprob[-2]) < self.thresh:
                self.converged_ = True
                break

            # Maximization step
            self._do_mstep(X, posteriors, params, self.min_covar)

        return self

    def _do_mstep(self, X, posteriors, params, min_covar=0):
            w = posteriors.sum(axis=0)
            avg_obs = np.dot(posteriors.T, X)
            norm = 1.0 / (w[:, np.newaxis] + 10 * np.finfo(np.float).eps)

            if 'w' in params:
                self._log_weights = np.log(w /
                        (w.sum() + 10 * np.finfo(np.float).eps)
                        + np.finfo(np.float).eps)
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
    lpr = -0.5 * (n_dim * np.log(2 * np.pi) + np.log(linalg.det(covars) + 0.1)
                  + np.sum(obs * np.dot(obs, icv), 1)[:, np.newaxis]
                  - 2 * np.dot(np.dot(obs, icv), means.T)
                  + np.sum(means * np.dot(means, icv), 1))
    return lpr


def _lmvnpdffull(obs, means, covars):
    """
    Log probability for full covariance matrices.

    WARNING: In certain cases, this function will modify in-place
    some of the covariance matrices

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
        try:
            cv_chol = linalg.cholesky(cv, lower=True)
        except linalg.LinAlgError:
            # The model is most probabily stuck in a component with too
            # few observations, we need to reinitialize this components
            cv[:] = 10 * np.eye(cv.shape[0])
            cv_chol = cv
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


def _distribute_covar_matrix_to_match_cvtype(tiedcv, cvtype, n_components):
    if cvtype == 'spherical':
        cv = np.tile(np.diag(tiedcv).mean(), n_components)
    elif cvtype == 'tied':
        cv = tiedcv
    elif cvtype == 'diag':
        cv = np.tile(np.diag(tiedcv), (n_components, 1))
    elif cvtype == 'full':
        cv = np.tile(tiedcv, (n_components, 1, 1))
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
    cv = np.empty((gmm.n_components, gmm.n_features, gmm.n_features))
    for c in xrange(gmm.n_components):
        post = posteriors[:, c]
        avg_cv = np.dot(post * obs.T, obs) / (post.sum() +
                                10 * np.finfo(np.float).eps)
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
    for c in xrange(gmm.n_components):
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
            covars += cv / gmm.n_components
    return covars


_covar_mstep_funcs = {'spherical': _covar_mstep_spherical,
                      'diag': _covar_mstep_diag,
                      #'tied': _covar_mstep_tied,
                      'full': _covar_mstep_full,
                      'tied': _covar_mstep_slow,
                      }
