"""
Gaussian Mixture Models.

This implementation corresponds to frequentist (non-Bayesian) formulation
of Gaussian Mixture Models.
"""

# Author: Ron Weiss <ronweiss@gmail.com>
#         Fabian Pedregosa <fabian.pedregosa@inria.fr>
#

import numpy as np

from ..base import BaseEstimator
from ..utils import check_random_state
from ..utils.extmath import logsumexp
from .. import cluster

INF_EPS = np.finfo(float).eps


def normalize(A, axis=None):
    """ Normalize the input array so that it sums to 1.
    
    Parameters
    ----------
    A: array, shape
    axis: int, dimension along which normalization is performed
    
    Returns
    -------
    normalized_A: A with values normalized along the prescribed axis

    WARNING: Modifies inplace the array
    """
    A += INF_EPS
    Asum = A.sum(axis)
    if axis and A.ndim > 1:
        # Make sure we don't divide by zero.
        Asum[Asum == 0] = 1
        shape = list(A.shape)
        shape[axis] = 1
        Asum.shape = shape
    return A / Asum


def log_multivariate_normal_density(X, means, covars, cvtype='diag'):
    """Compute the log probability under a multivariate Gaussian distribution.

    Parameters
    ----------
    X : array_like, shape (n_samples, n_features)
        List of n_features-dimensional data points.  Each row corresponds to a
        single data point.
    means : array_like, shape (n_components, n_features)
        List of n_features-dimensional mean vectors for n_components Gaussians.
        Each row corresponds to a single mean vector.
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
    lpr : array_like, shape (n_samples, n_components)
        Array containing the log probabilities of each data point in
        X under each of the n_components multivariate Gaussian distributions.
    """
    log_multivariate_normal_density_dict = {
        'spherical': _log_multivariate_normal_density_spherical,
        'tied': _log_multivariate_normal_density_tied,
        'diag': _log_multivariate_normal_density_diag,
        'full': _log_multivariate_normal_density_full}
    return log_multivariate_normal_density_dict[cvtype](X, means, covars)


def sample_gaussian(mean, covar, cvtype='diag', n_samples=1,
                    random_state=None):
    """Generate random samples from a Gaussian distribution.

    Parameters
    ----------
    mean : array_like, shape (n_features,)
        Mean of the distribution.

    covars : array_like, optional
        Covariance of the distribution.  The shape depends on `cvtype`:
            scalar if 'spherical',
            (n_features) if 'diag',
            (n_features, n_features)  if 'tied', or 'full'

    cvtype : string, optional
        Type of the covariance parameters.  Must be one of
        'spherical', 'tied', 'diag', 'full'.  Defaults to 'diag'.

    n_samples : int, optional
        Number of samples to generate. Defaults to 1.

    Returns
    -------
    X : array, shape (n_features, n_samples)
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
    cvtype : string (read-only)
        String describing the type of covariance parameters used by
        the GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.
    weights : array, shape (`n_components`,)
        Mixing weights for each mixture component.
    means : array, shape (`n_components`, `n_features`)
        Mean parameters for each mixture component.
    covars : array
        Covariance parameters for each mixture component.  The shape
        depends on `cvtype`:
            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'
    converged_ : bool
        True when convergence was reached in fit(), False
        otherwise.

    Methods
    -------
    eval(X) -> predict_proba
        Compute the log likelihood of X under the model and the
        posterior distribution over mixture components.
    fit(X)
        Estimate model parameters from X using the EM algorithm.
    predict(X)
        Find most likely mixtures components for each observation in X.
    rvs(n=1, random_state=None)
        Generate `n` samples from the model.
    score(X)
        Compute the log likelihood of X under the model.

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

        Must be one of 'spherical', 'tied', 'diag', 'full'.
        """
        return self._cvtype

    def _get_covars(self):
        """Return covars as a full matrix."""
        if self.cvtype == 'full':
            return self.covars_
        elif self.cvtype == 'diag':
            return [np.diag(cov) for cov in self.covars_]
        elif self.cvtype == 'tied':
            return [self.covars_] * self.n_components
        elif self.cvtype == 'spherical':
            return [np.eye(self.means.shape[1]) * f for f in self.covars_]

    def _set_covars(self, covars):
        covars = np.asarray(covars)
        _validate_covars(covars, self._cvtype, self.n_components)
        self.covars_ = covars

    covars = property(_get_covars, _set_covars)

    def _get_means(self):
        """Mean parameters for each mixture component."""
        return self.means_

    def _set_means(self, means):
        """Provide values for means"""
        means = np.asarray(means)
        if means.shape[0] != self.n_components:
            raise ValueError('means must have shape ' +
                    '(n_components, n_features)')
        self.means_ = means.copy()

    means = property(_get_means, _set_means)

    def __repr__(self):
        return "GMM(cvtype='%s', n_components=%s)" % (self._cvtype,
                self.n_components)

    def _get_weights(self):
        """Mixing weights for each mixture component."""
        return np.exp(self.log_weights_)

    def _set_weights(self, weights):
        """Provide value for micture weights"""
        if len(weights) != self.n_components:
            raise ValueError('weights must have length n_components')
        if not np.allclose(np.sum(weights), 1.0):
            raise ValueError('weights must sum to 1.0')

        self.log_weights_ = np.log(np.asarray(weights).copy())

    weights = property(_get_weights, _set_weights)

    def eval(self, X):
        """Evaluate the model on data

        Compute the log probability of X under the model and
        return the posterior distribution (responsibilities) of each
        mixture component for each element of X.

        Parameters
        ----------
        X: array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob: array_like, shape (n_samples,)
            Log probabilities of each data point in X
        posteriors: array_like, shape (n_samples, n_components)
            Posterior probabilities of each mixture component for each
            observation
        """
        X = np.asarray(X)
        lpr = (log_multivariate_normal_density(
                X, self.means_, self.covars_, self._cvtype) + self.log_weights_)
        logprob = logsumexp(lpr, axis=1)
        posteriors = np.exp(lpr - logprob[:, np.newaxis])
        return logprob, posteriors

    def score(self, X):
        """Compute the log probability under the model.

        Parameters
        ----------
        X : array_like, shape (n_samples, n_features)
            List of n_features-dimensional data points.  Each row
            corresponds to a single data point.

        Returns
        -------
        logprob : array_like, shape (n_samples,)
            Log probabilities of each data point in X
        """
        logprob, _ = self.eval(X)
        return logprob

    def predict(self, X):
        """Predict label for data.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = (n_samples,)
        """
        logprob, posteriors = self.eval(X)
        return posteriors.argmax(axis=1)

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
        X : array_like, shape (n_samples, n_features)
            List of samples
        """
        if random_state is None:
            random_state = self.random_state
        random_state = check_random_state(random_state)
        weight_pdf = self.weights
        weight_cdf = np.cumsum(weight_pdf)

        X = np.empty((n_samples, self.means.shape[1]))
        rand = random_state.rand(n_samples)
        # decide which component to use for each sample
        comps = weight_cdf.searchsorted(rand)
        # for each component, generate all needed samples
        for comp in xrange(self.n_components):
            # occurrences of current component in X
            comp_in_X = (comp == comps)
            # number of those occurrences
            num_comp_in_X = comp_in_X.sum()
            if num_comp_in_X > 0:
                if self._cvtype == 'tied':
                    cv = self.covars_
                else:
                    cv = self.covars_[comp]
                X[comp_in_X] = sample_gaussian(
                    self.means_[comp], cv, self._cvtype, num_comp_in_X,
                    random_state=random_state).T
        return X

    def fit(self, X, n_iter=100, n_init=1, thresh=1e-2, params='wmc',
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
        n_init : int, optional
            number of initializations to perform. the best results is kept
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
        if X.ndim == 1:
            X = X[:, np.newaxis]
        if X.shape[0] < self.n_components:
            raise ValueError('GMM estimation with %s components, but '
                             'got only %s' % (self.n_components, X.shape[0]))
        n_features = X.shape[1]
        
        max_log_prob = - np.infty
        if n_init < 1:
            raise ValueError('GMM estimation requires at least one run')
        for _ in range(n_init):
            if 'm' in init_params or not hasattr(self, 'means'):
                self.means_ = cluster.KMeans(
                    k=self.n_components).fit(X).cluster_centers_

            if 'w' in init_params or not hasattr(self, 'weights'):
                self.weights = np.tile(1.0 / self.n_components, 
                                       self.n_components)

            if 'c' in init_params or not hasattr(self, 'covars'):
                cv = np.cov(X.T)
                if not cv.shape:
                    cv.shape = (1, 1)
                self.covars_ = _distribute_covar_matrix_to_match_cvtype(
                    cv, self._cvtype, self.n_components)

            # EM algorithm
            log_likelihood = []
            # reset self.converged_ to False
            self.converged_ = False
            for i in xrange(n_iter):
                # Expectation step
                curr_log_likelihood, posteriors = self.eval(X)
                log_likelihood.append(curr_log_likelihood.sum())

                # Check for convergence.
                if i > 0 and abs(log_likelihood[-1] - log_likelihood[-2]) < \
                        self.thresh:
                    self.converged_ = True
                    break

                # Maximization step
                self._do_mstep(X, posteriors, params, self.min_covar)

            # if the results is better, keep it   
            if log_likelihood[-1] > max_log_prob:
                max_log_prob = log_likelihood[-1]
                best_params = {'weights': self.weights, 
                               'means': self.means_, 
                               'covars': self.covars_}
        self.covars_ = best_params['covars']
        self.means_ = best_params['means']
        self.weights = best_params['weights']
        return self

    def _do_mstep(self, X, posteriors, params, min_covar=0):
        """ Perform the Mstep of the EM algorithm and return the class weihgts.
        """
        weights = posteriors.sum(axis=0)
        weighted_X_sum = np.dot(posteriors.T, X)
        
        inverse_weights = 1.0 / (
            weights[:, np.newaxis] + 10 * INF_EPS)
        
        if 'w' in params:
            self.log_weights_ = np.log(
                weights / (weights.sum() + 10 * INF_EPS) + INF_EPS)
        if 'm' in params:
            self.means_ = weighted_X_sum * inverse_weights
        if 'c' in params:
            covar_mstep_func = _covar_mstep_funcs[self._cvtype]
            self.covars_ = covar_mstep_func(
                self, X, posteriors, weighted_X_sum, inverse_weights, min_covar)

        return weights


#########################################################################
## some helper routines
#########################################################################


def _log_multivariate_normal_density_diag(X, means=0.0, covars=1.0):
    """Compute Gaussian density at X for a diagonal model"""
    n_samples, n_dim = X.shape
    lpr = -0.5 * (n_dim * np.log(2 * np.pi) + np.sum(np.log(covars), 1)
                  + np.sum((means ** 2) / covars, 1)
                  - 2 * np.dot(X, (means / covars).T)
                  + np.dot(X ** 2, (1.0 / covars).T))
    return lpr


def _log_multivariate_normal_density_spherical(X, means=0.0, covars=1.0):
    """Compute Gaussian density at X for a spherical model"""
    cv = covars.copy()
    if covars.ndim == 1:
        cv = cv[:, np.newaxis]
    return _log_multivariate_normal_density_diag(X, means, 
                                                 np.tile(cv, (1, X.shape[-1])))


def _log_multivariate_normal_density_tied(X, means, covars):
    """Compute Gaussian density at X for a tied model"""
    from scipy import linalg
    n_samples, n_dim = X.shape
    icv = linalg.pinv(covars)
    lpr = -0.5 * (n_dim * np.log(2 * np.pi) + np.log(linalg.det(covars) + 0.1)
                  + np.sum(X * np.dot(X, icv), 1)[:, np.newaxis]
                  - 2 * np.dot(np.dot(X, icv), means.T)
                  + np.sum(means * np.dot(means, icv), 1))
    return lpr


def _log_multivariate_normal_density_full(X, means, covars):
    """Log probability for full covariance matrices.

    WARNING
    -------
    In certain cases, this function will modify in-place
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
    n_samples, n_dim = X.shape
    nmix = len(means)
    log_prob = np.empty((n_samples, nmix))
    for c, (mu, cv) in enumerate(itertools.izip(means, covars)):
        try:
            cv_chol = linalg.cholesky(cv, lower=True)
        except linalg.LinAlgError:
            # The model is most probabily stuck in a component with too
            # few observations, we need to reinitialize this components
            cv[:] = np.diag(np.var(X, 0))
            cv_chol = np.diag(np.std(X, 0))
        cv_log_det = 2 * np.sum(np.log(np.diagonal(cv_chol)))
        cv_sol = solve_triangular(cv_chol, (X - mu).T, lower=True).T
        log_prob[:, c] = -.5 * (np.sum(cv_sol ** 2, axis=1) + \
                           n_dim * np.log(2 * np.pi) + cv_log_det)

    return log_prob


def _validate_covars(covars, cvtype, n_components):
    """Do basic checks on matrix covariance sizes and values
    """
    from scipy import linalg
    if cvtype == 'spherical':
        if len(covars) != n_components:
            raise ValueError("'spherical' covars must have length n_components")
        elif np.any(covars <= 0):
            raise ValueError("'spherical' covars must be non-negative")
    elif cvtype == 'tied':
        if covars.shape[0] != covars.shape[1]: 
            raise ValueError("'tied' covars must have shape (n_dim, n_dim)")
        elif (not np.allclose(covars, covars.T)
              or np.any(linalg.eigvalsh(covars) <= 0)):
            raise ValueError("'tied' covars must be symmetric, "
                             "positive-definite")
    elif cvtype == 'diag':
        if len(covars.shape) != 2:
            raise ValueError("'diag' covars must have shape"
                             "(n_components, n_dim)")
        elif np.any(covars <= 0):
            raise ValueError("'diag' covars must be non-negative")
    elif cvtype == 'full':
        if len(covars.shape) != 3:
            raise ValueError("'full' covars must have shape "
                             "(n_components, n_dim, n_dim)")
        elif covars.shape[2] != covars.shape[3]:
            raise ValueError("'full' covars must have shape "
                             "(n_components, n_dim, n_dim)")
        for n, cv in enumerate(covars):
            if (not np.allclose(cv, cv.T)
                or np.any(linalg.eigvalsh(cv) <= 0)):
                raise ValueError("component %d of 'full' covars must be "
                                 "symmetric, positive-definite" % n)


def _distribute_covar_matrix_to_match_cvtype(tiedcv, cvtype, n_components):
    """Create all the covariance matrices from a given template 
    """
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


def _covar_mstep_diag(gmm, X, posteriors, weighted_X_sum, norm, min_covar):
    """Performing the covariance M step for diagonal cases"""
    avg_X2 = np.dot(posteriors.T, X * X) * norm
    avg_means2 = gmm.means_ ** 2
    avg_X_means = gmm.means_ * weighted_X_sum * norm
    return avg_X2 - 2 * avg_X_means + avg_means2 + min_covar


def _covar_mstep_spherical(*args):
    """Performing the covariance M step for spherical cases"""
    return _covar_mstep_diag(*args).mean(axis=1)


def _covar_mstep_full(gmm, X, posteriors, weighted_X_sum, norm, min_covar):
    """Performing the covariance M step for full cases"""
    # Eq. 12 from K. Murphy, "Fitting a Conditional Linear Gaussian
    # Distribution"
    n_features = X.shape[1]
    cv = np.empty((gmm.n_components, n_features, n_features))
    for c in xrange(gmm.n_components):
        post = posteriors[:, c]
        # Underflow Errors in doing  post * X.T are not important
        err_mgt = np.seterr(under='ignore')
        avg_cv = np.dot(post * X.T, X) / (post.sum() + 10 * INF_EPS)
        mu = gmm.means_[c][np.newaxis]
        cv[c] = (avg_cv - np.dot(mu.T, mu) + min_covar * np.eye(n_features))
    return cv


def _covar_mstep_tied2(*args):
    return _covar_mstep_full(*args).mean(axis=0)


def _covar_mstep_tied(gmm, X, posteriors, weighted_X_sum, norm, min_covar):
    # Eq. 15 from K. Murphy, "Fitting a Conditional Linear Gaussian
    n_features = X.shape[1]
    avg_X2 = np.dot(X.T, X)
    avg_means2 = np.dot(gmm.means_.T, weighted_X_sum)
    return (avg_X2 - avg_means2 + min_covar * np.eye(n_features)) / X.shape[0]


_covar_mstep_funcs = {'spherical': _covar_mstep_spherical,
                      'diag': _covar_mstep_diag,
                      'tied': _covar_mstep_tied,
                      'full': _covar_mstep_full,
                      }
