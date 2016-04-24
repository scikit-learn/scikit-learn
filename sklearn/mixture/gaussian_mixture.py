"""Gaussian Mixture Model."""

# Author: Wei Xue <xuewei4d@gmail.com>
# Modified by Thierry Guillemot <thierry.guillemot.work@gmail.com>

import numpy as np

from scipy import linalg

from .base import BaseMixture, _check_shape
from ..externals.six.moves import zip
from ..utils import check_array
from ..utils.validation import check_is_fitted


###############################################################################
# Gaussian mixture shape checkers used by the GaussianMixture class

def _check_weights(weights, n_components):
    """Check the user provided 'weights'.

    Parameters
    ----------
    weights : array-like, shape (n_components,)
        The proportions of components of each mixture.

    n_components : int
        Number of components.

    Returns
    -------
    weights : array, shape (n_components,)
    """
    weights = check_array(weights, dtype=[np.float64, np.float32],
                          ensure_2d=False)
    _check_shape(weights, (n_components,), 'weights')

    # check range
    if (any(np.less(weights, 0)) or
            any(np.greater(weights, 1))):
        raise ValueError("The parameter 'weights' should be in the range "
                         "[0, 1], but got max value %.5f, min value %.5f"
                         % (np.min(weights), np.max(weights)))

    # check normalization
    if not np.allclose(np.abs(1 - np.sum(weights)), 0.0):
        raise ValueError("The parameter 'weights' should be normalized, "
                         "but got sum(weights) = %.5f" % np.sum(weights))
    return weights


def _check_means(means, n_components, n_features):
    """Validate the provided 'means'.

    Parameters
    ----------
    means : array-like, shape (n_components, n_features)
        The centers of the current components.

    n_components : int
        Number of components.

    n_features : int
        Number of features.

    Returns
    -------
    means : array, (n_components, n_features)
    """
    means = check_array(means, dtype=[np.float64, np.float32], ensure_2d=False)
    _check_shape(means, (n_components, n_features), 'means')
    return means


def _check_covariance_matrix(covariance, covariance_type):
    """Check a covariance matrix is symmetric and positive-definite."""
    if (not np.allclose(covariance, covariance.T) or
            np.any(np.less_equal(linalg.eigvalsh(covariance), .0))):
        raise ValueError("'%s covariance' should be symmetric, "
                         "positive-definite" % covariance_type)


def _check_covariance_positivity(covariance, covariance_type):
    """Check a covariance vector is positive-definite."""
    if np.any(np.less_equal(covariance, 0.0)):
        raise ValueError("'%s covariance' should be "
                         "positive" % covariance_type)


def _check_covariances_full(covariances, covariance_type):
    """Check the covariance matrices are symmetric and positive-definite."""
    for k, cov in enumerate(covariances):
        _check_covariance_matrix(cov, covariance_type)


def _check_covariances(covariances, covariance_type, n_components, n_features):
    """Validate user provided covariances.

    Parameters
    ----------
    covariances : array-like,
        'full' : shape of (n_components, n_features, n_features)
        'tied' : shape of (n_features, n_features)
        'diag' : shape of (n_components, n_features)
        'spherical' : shape of (n_components,)

    covariance_type : string

    n_components : int
        Number of components.

    n_features : int
        Number of features.

    Returns
    -------
    covariances : array
    """
    covariances = check_array(covariances, dtype=[np.float64, np.float32],
                              ensure_2d=False,
                              allow_nd=covariance_type is 'full')

    covariances_shape = {'full': (n_components, n_features, n_features),
                         'tied': (n_features, n_features),
                         'diag': (n_components, n_features),
                         'spherical': (n_components,)}
    _check_shape(covariances, covariances_shape[covariance_type],
                 '%s covariance' % covariance_type)

    check_functions = {'full': _check_covariances_full,
                       'tied': _check_covariance_matrix,
                       'diag': _check_covariance_positivity,
                       'spherical': _check_covariance_positivity}
    check_functions[covariance_type](covariances, covariance_type)

    return covariances


###############################################################################
# Gaussian mixture parameters estimators (used by the M-Step)

def _estimate_gaussian_covariance_full(resp, X, nk, means, reg_covar):
    """Estimate the full covariance matrices.

    Parameters
    ----------
    resp : array-like, shape (n_samples, n_components)

    X : array-like, shape (n_samples, n_features)

    nk : array-like, shape (n_components,)

    means : array-like, shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_components, n_features, n_features)
    """
    n_features = X.shape[1]
    n_components = means.shape[0]
    covariances = np.empty((n_components, n_features, n_features))
    for k in range(n_components):
        diff = X - means[k]
        covariances[k] = np.dot(resp[:, k] * diff.T, diff) / nk[k]
        covariances[k].flat[::n_features + 1] += reg_covar
    return covariances


def _estimate_gaussian_covariance_tied(resp, X, nk, means, reg_covar):
    """Estimate the tied covariance matrix.

    Parameters
    ----------
    resp : array-like, shape (n_samples, n_components)

    X : array-like, shape (n_samples, n_features)

    nk : array-like, shape (n_components,)

    means : array-like, shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_features, n_features)
    """
    avg_X2 = np.dot(X.T, X)
    avg_means2 = np.dot(nk * means.T, means)
    covariances = avg_X2 - avg_means2
    covariances /= X.shape[0]
    covariances.flat[::len(covariances) + 1] += reg_covar
    return covariances


def _estimate_gaussian_covariance_diag(resp, X, nk, means, reg_covar):
    """Estimate the diagonal covariance matrices.

    Parameters
    ----------
    responsibilities : array-like, shape (n_samples, n_components)

    X : array-like, shape (n_samples, n_features)

    nk : array-like, shape (n_components,)

    means : array-like, shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_components, n_features)
    """
    avg_X2 = np.dot(resp.T, X * X) / nk[:, np.newaxis]
    avg_means2 = means ** 2
    avg_X_means = means * np.dot(resp.T, X) / nk[:, np.newaxis]
    return avg_X2 - 2 * avg_X_means + avg_means2 + reg_covar


def _estimate_gaussian_covariance_spherical(resp, X, nk, means, reg_covar):
    """Estimate the spherical covariance matrices.

    Parameters
    ----------
    responsibilities : array-like, shape (n_samples, n_components)

    X : array-like, shape (n_samples, n_features)

    nk : array-like, shape (n_components,)

    means : array-like, shape (n_components, n_features)

    reg_covar : float

    Returns
    -------
    covariances : array, shape (n_components,)
    """
    covariances = _estimate_gaussian_covariance_diag(resp, X, nk, means,
                                                     reg_covar)
    return covariances.mean(axis=1)


def _estimate_gaussian_parameters(X, resp, reg_covar, covariance_type):
    """Estimate the Gaussian distribution parameters.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The input data array.

    resp : array-like, shape (n_samples, n_features)
        The responsibilities for each data sample in X.

    reg_covar : float
        The regularization added to each covariance matrices.

    covariance_type : {'full', 'tied', 'diag', 'spherical'}
        The type of covariance matrices.

    Returns
    -------
    nk : array, shape (n_components,)
        The numbers of data samples in the current components.

    means : array, shape (n_components, n_features)
        The centers of the current components.

    covariances : array
        The sample covariances of the current components.
        The shape depends of the covariance_type.
    """
    compute_covariance = {
        "full": _estimate_gaussian_covariance_full,
        "tied": _estimate_gaussian_covariance_tied,
        "diag": _estimate_gaussian_covariance_diag,
        "spherical": _estimate_gaussian_covariance_spherical}

    nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
    means = np.dot(resp.T, X) / nk[:, np.newaxis]
    covariances = compute_covariance[covariance_type](
        resp, X, nk, means, reg_covar)

    return nk, means, covariances


###############################################################################
# Gaussian mixture probability estimators

def _estimate_log_gaussian_prob_full(X, means, covariances):
    """Estimate the log Gaussian probability for 'full' covariance.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    means : array-like, shape (n_components, n_features)

    covariances : array-like, shape (n_components, n_features, n_features)

    Returns
    -------
    log_prob : array, shape (n_samples, n_components)
    """
    n_samples, n_features = X.shape
    n_components = means.shape[0]
    log_prob = np.empty((n_samples, n_components))
    for k, (mu, cov) in enumerate(zip(means, covariances)):
        try:
            cov_chol = linalg.cholesky(cov, lower=True)
        except linalg.LinAlgError:
            raise ValueError("The algorithm has diverged because of too "
                             "few samples per components. "
                             "Try to decrease the number of components, or "
                             "increase reg_covar.")
        cv_log_det = 2. * np.sum(np.log(np.diagonal(cov_chol)))
        cv_sol = linalg.solve_triangular(cov_chol, (X - mu).T, lower=True).T
        log_prob[:, k] = - .5 * (n_features * np.log(2. * np.pi) +
                                 cv_log_det +
                                 np.sum(np.square(cv_sol), axis=1))
    return log_prob


def _estimate_log_gaussian_prob_tied(X, means, covariances):
    """Estimate the log Gaussian probability for 'tied' covariance.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    means : array-like, shape (n_components, n_features)

    covariances : array-like, shape (n_features, n_features)

    Returns
    -------
    log_prob : array-like, shape (n_samples, n_components)
    """
    n_samples, n_features = X.shape
    n_components = means.shape[0]
    log_prob = np.empty((n_samples, n_components))
    try:
        cov_chol = linalg.cholesky(covariances, lower=True)
    except linalg.LinAlgError:
        raise ValueError("The algorithm has diverged because of too "
                         "few samples per components. "
                         "Try to decrease the number of components, or "
                         "increase reg_covar.")
    cv_log_det = 2. * np.sum(np.log(np.diagonal(cov_chol)))
    for k, mu in enumerate(means):
        cv_sol = linalg.solve_triangular(cov_chol, (X - mu).T,
                                         lower=True).T
        log_prob[:, k] = np.sum(np.square(cv_sol), axis=1)
    log_prob = - .5 * (n_features * np.log(2. * np.pi) + cv_log_det + log_prob)
    return log_prob


def _estimate_log_gaussian_prob_diag(X, means, covariances):
    """Estimate the log Gaussian probability for 'diag' covariance.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    means : array-like, shape (n_components, n_features)

    covariances : array-like, shape (n_components, n_features)

    Returns
    -------
    log_prob : array-like, shape (n_samples, n_components)
    """
    if np.any(np.less_equal(covariances, 0.0)):
        raise ValueError("The algorithm has diverged because of too "
                         "few samples per components. "
                         "Try to decrease the number of components, or "
                         "increase reg_covar.")
    n_samples, n_features = X.shape
    log_prob = - .5 * (n_features * np.log(2. * np.pi) +
                       np.sum(np.log(covariances), 1) +
                       np.sum((means ** 2 / covariances), 1) -
                       2. * np.dot(X, (means / covariances).T) +
                       np.dot(X ** 2, (1. / covariances).T))
    return log_prob


def _estimate_log_gaussian_prob_spherical(X, means, covariances):
    """Estimate the log Gaussian probability for 'spherical' covariance.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)

    means : array-like, shape (n_components, n_features)

    covariances : array-like, shape (n_components, )

    Returns
    -------
    log_prob : array-like, shape (n_samples, n_components)
    """
    if np.any(np.less_equal(covariances, 0.0)):
        raise ValueError("The algorithm has diverged because of too "
                         "few samples per components. "
                         "Try to decrease the number of components, or "
                         "increase reg_covar.")
    n_samples, n_features = X.shape
    log_prob = - .5 * (n_features * np.log(2 * np.pi) +
                       n_features * np.log(covariances) +
                       np.sum(means ** 2, 1) / covariances -
                       2 * np.dot(X, means.T / covariances) +
                       np.outer(np.sum(X ** 2, axis=1), 1. / covariances))
    return log_prob


class GaussianMixture(BaseMixture):
    """Gaussian Mixture.

    Representation of a Gaussian mixture model probability distribution.
    This class allows to estimate the parameters of a Gaussian mixture
    distribution.

    Parameters
    ----------
    n_components : int, defaults to 1.
        The number of mixture components.

    covariance_type : {'full', 'tied', 'diag', 'spherical'},
        defaults to 'full'.
        String describing the type of covariance parameters to use.
        Must be one of::
        'full' (each component has its own general covariance matrix).
        'tied' (all components share the same general covariance matrix),
        'diag' (each component has its own diagonal covariance matrix),
        'spherical' (each component has its own single variance),

    tol : float, defaults to 1e-3.
        The convergence threshold. EM iterations will stop when the
        log_likelihood average gain is below this threshold.

    reg_covar : float, defaults to 0.
        Non-negative regularization added to the diagonal of covariance.
        Allows to assure that the covariance matrices are all positive.

    max_iter : int, defaults to 100.
        The number of EM iterations to perform.

    n_init : int, defaults to 1.
        The number of initializations to perform. The best results is kept.

    init_params : {'kmeans', 'random'}, defaults to 'kmeans'.
        The method used to initialize the weights, the means and the
        covariances.
        Must be one of::
        'kmeans' : responsibilities are initialized using kmeans.
        'random' : responsibilities are initialized randomly.

    weights_init : array-like, shape (n_components, ), optional
        The user-provided initial weights, defaults to None.
        If it None, weights are initialized using the `init_params` method.

    means_init: array-like, shape (n_components, n_features), optional
        The user-provided initial means, defaults to None,
        If it None, means are initialized using the `init_params` method.

    covariances_init: array-like, optional.
        The user-provided initial covariances, defaults to None.
        If it None, covariances are initialized using the 'init_params' method.
        The shape depends on 'covariance_type'::
            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'

    random_state: RandomState or an int seed, defaults to None.
        A random number generator instance.

    warm_start : bool, default to False.
        If 'warm_start' is True, the solution of the last fitting is used as
        initialization for the next call of fit(). This can speed up
        convergence when fit is called several time on similar problems.

    verbose : int, default to 0.
        Enable verbose output. If 1 then it prints the current
        initialization and each iteration step. If greater than 1 then
        it prints also the log probability and the time needed
        for each step.

    Attributes
    ----------
    weights_ : array, shape (n_components,)
        The weights of each mixture components.
        `weights_` will not exist before a call to fit.

    means_ : array, shape (n_components, n_features)
        The mean of each mixture component.
        `means_` will not exist before a call to fit.

    covariances_ : array
        The covariance of each mixture component.
        The shape depends on `covariance_type`::
            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'
        `covariances_` will not exist before a call to fit.

    converged_ : bool
        True when convergence was reached in fit(), False otherwise.
        `converged_` will not exist before a call to fit.

    n_iter_ : int
        Number of step used by the best fit of EM to reach the convergence.
        `n_iter_`  will not exist before a call to fit.
    """

    def __init__(self, n_components=1, covariance_type='full', tol=1e-3,
                 reg_covar=1e-6, max_iter=100, n_init=1, init_params='kmeans',
                 weights_init=None, means_init=None, covariances_init=None,
                 random_state=None, warm_start=False,
                 verbose=0, verbose_interval=10):
        super(GaussianMixture, self).__init__(
            n_components=n_components, tol=tol, reg_covar=reg_covar,
            max_iter=max_iter, n_init=n_init, init_params=init_params,
            random_state=random_state, warm_start=warm_start,
            verbose=verbose, verbose_interval=verbose_interval)

        self.covariance_type = covariance_type
        self.weights_init = weights_init
        self.means_init = means_init
        self.covariances_init = covariances_init

    def _check_parameters(self, X):
        """Check the Gaussian mixture parameters are well defined."""
        if self.covariance_type not in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError("Invalid value for 'covariance_type': %s "
                             "'covariance_type' should be in "
                             "['spherical', 'tied', 'diag', 'full']"
                             % self.covariance_type)

        if self.weights_init is not None:
            self.weights_init = _check_weights(self.weights_init,
                                               self.n_components)

        if self.means_init is not None:
            self.means_init = _check_means(self.means_init,
                                           self.n_components, X.shape[1])

        if self.covariances_init is not None:
            self.covariances_init = _check_covariances(self.covariances_init,
                                                       self.covariance_type,
                                                       self.n_components,
                                                       X.shape[1])

    def _initialize(self, X, resp):
        """Initialization of the Gaussian mixture parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        resp : array-like, shape (n_samples, n_components)
        """
        weights, means, covariances = _estimate_gaussian_parameters(
            X, resp, self.reg_covar, self.covariance_type)
        weights /= X.shape[0]

        self.weights_ = (weights if self.weights_init is None
                         else self.weights_init)
        self.means_ = means if self.means_init is None else self.means_init
        self.covariances_ = (covariances if self.covariances_init is None
                             else self.covariances_init)

    def _e_step(self, X):
        log_prob_norm, _, log_resp = self._estimate_log_prob_resp(X)
        return np.mean(log_prob_norm), np.exp(log_resp)

    def _m_step(self, X, resp):
        self.weights_, self.means_, self.covariances_ = (
            _estimate_gaussian_parameters(X, resp, self.reg_covar,
                                          self.covariance_type))
        self.weights_ /= X.shape[0]

    def _estimate_log_prob(self, X):
        estimate_log_prob_functions = {
            "full": _estimate_log_gaussian_prob_full,
            "tied": _estimate_log_gaussian_prob_tied,
            "diag": _estimate_log_gaussian_prob_diag,
            "spherical": _estimate_log_gaussian_prob_spherical
        }
        return estimate_log_prob_functions[self.covariance_type](
            X, self.means_, self.covariances_)

    def _estimate_log_weights(self):
        return np.log(self.weights_)

    def _check_is_fitted(self):
        check_is_fitted(self, ['weights_', 'means_', 'covariances_'])

    def _get_parameters(self):
        return self.weights_, self.means_, self.covariances_

    def _set_parameters(self, params):
        self.weights_, self.means_, self.covariances_ = params

    def _n_parameters(self):
        """Return the number of free parameters in the model."""
        ndim = self.means_.shape[1]
        if self.covariance_type == 'full':
            cov_params = self.n_components * ndim * (ndim + 1) / 2.
        elif self.covariance_type == 'diag':
            cov_params = self.n_components * ndim
        elif self.covariance_type == 'tied':
            cov_params = ndim * (ndim + 1) / 2.
        elif self.covariance_type == 'spherical':
            cov_params = self.n_components
        mean_params = ndim * self.n_components
        return int(cov_params + mean_params + self.n_components - 1)

    def bic(self, X):
        """Bayesian information criterion for the current model on the input X.

        Parameters
        ----------
        X : array of shape (n_samples, n_dimensions)

        Returns
        -------
        bic: float
            The greater the better.
        """
        return (-2 * self.score(X) * X.shape[0] +
                self._n_parameters() * np.log(X.shape[0]))

    def aic(self, X):
        """Akaike information criterion for the current model on the input X.

        Parameters
        ----------
        X : array of shape(n_samples, n_dimensions)

        Returns
        -------
        aic: float
            The greater the better.
        """
        return -2 * self.score(X) * X.shape[0] + 2 * self._n_parameters()
