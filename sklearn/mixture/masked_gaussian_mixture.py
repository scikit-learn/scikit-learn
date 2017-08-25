from __future__ import division
import warnings
import numpy as np

from .gaussian_mixture import GaussianMixture, _compute_precision_cholesky, \
    _estimate_gaussian_parameters
from .base import _check_X
from ..exceptions import ConvergenceWarning
from ..utils import check_array, check_random_state


def _check_empty_axis(X_mask):
    """Use the input data mask X_mask to check for empty rows and/or columns.

    Parameters
    ----------
    X_mask : array-like, shape (n_samples, n_features)
    """
    if X_mask.sum(axis=0) == X_mask.shape[0]:
        raise ValueError("Please remove empty column(s) present in the "
                         "dataset.")
    if X_mask.sum(axis=1) == X_mask.shape[1]:
        raise ValueError("Please remove empty row(s) present in the "
                         "dataset.")
    return X_mask


def _check_inf(X):
    """Check for +/- inf in the input data X.

    Parameters
    ----------
    X_mask : array-like, shape (n_samples, n_features)
    """
    if np.any(np.isinf(X)):
        raise ValueError("+/- inf values are not allowed.")
    return X


def _get_mask(X, value_to_mask):
    """Compute the boolean mask X == missing_values."""
    if value_to_mask == "NaN" or np.isnan(value_to_mask):
        return np.isnan(X)
    else:
        return X == value_to_mask


# See Rubin (2002) p. 148
# The _sweep() function below is an adaptation of the code by Jaime at
# https://stackoverflow.com/questions/15767435
def _sweep(g, k_index, reverse=False):
    r = -1 if reverse else 1
    k_index = np.array(
        [k_index]) if type(k_index) == int else np.array(k_index)
    for k in k_index:
        h = g - np.outer(g[:, k], g[k, :]) / g[k, k]
        h[:, k] = r * g[:, k] / g[k, k]
        h[k, :] = h[:, k]
        h[k, k] = -1.0 / g[k, k]
        g = h
    return g


def _I(s, I_all):
    return np.where(I_all == s)[0]

# Returns observed columns in pattern s
def _O(s, R):
    return np.where(R[s, ])[0]

# Returns missing columns in pattern s
def _M(s, R):
    return np.where(~R[s, ])[0]


def _exp_sufficient_stats(X, theta, T, S):
    # Based on the algorithm by Schafer (1997)
    _, p = X.shape  # Feature count
    c = np.zeros(p)
    for s in range(S):
        for j in range(0, p):
            if j in _O(s) and np.diag(theta)[j + 1] > 0:  # j to j+1
                theta = _sweep(theta, j + 1)
            if j in _M(s) and np.diag(theta)[j + 1] < 0:
                theta = _sweep(theta, j + 1, reverse=True)
        for i in _I(s):
            for j in _M(s):
                c[j] = theta[0, j + 1]
                for k in _O(s):
                    c[j] = c[j] + theta[k + 1, j + 1] * X[i, k]  # CHECK
            for j in _M(s):
                T[0, j + 1] = T[0, j + 1] + c[j]
                T[j + 1, 0] = T[0, j + 1]  # Added to maintain symmetry
                for k in _O(s):
                    T[k + 1, j + 1] = T[k + 1, j + 1] + c[j] * X[i, k]
                    T[j + 1, k + 1] = T[k + 1, j + 1]  # Added
                for k in _M(s):
                    if k >= j:
                        T[k + 1, j + 1] = T[k + 1, j + 1] + theta[
                            k + 1, j + 1] + c[k] * c[j]
                        T[j + 1, k + 1] = T[k + 1, j + 1]  # Added
    return T


def _estimate_masked_gaussian_parameters(X, resp, reg_covar,
                                         covariance_type, sufficient_stats):
    """Estimate the Gaussian distribution parameters.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        The input data array.

    resp : array-like, shape (n_samples, n_components)
        The responsibilities for each data sample in X.

    reg_covar : float
        The regularization added to the diagonal of the covariance matrices.

    covariance_type : {'full', 'tied', 'diag', 'spherical'}
        The type of precision matrices.

    Returns
    -------
    nk : array-like, shape (n_components,)
        The numbers of data samples in the current components.

    means : array-like, shape (n_components, n_features)
        The centers of the current components.

    covariances : array-like
        The covariance matrix of the current components.
        The shape depends of the covariance_type.
    """
    nk = resp.sum(axis=0) + 10 * np.finfo(resp.dtype).eps
    _, n_components = resp.shape
    n_features = sufficient_stats.shape[1] - 1
    means = np.empty((n_components, n_features))
    covariances = np.empty((n_components, n_features, n_features))
    # Note: Following returns full covariance matrices
    # TODO: Add support for tied, diag, and spherical covariances
    for comp in range(n_components):
        sufficient_stats_div_nk = sufficient_stats[comp] / nk
        augmented_covariance = _sweep(sufficient_stats_div_nk, 0)
        means[comp] = augmented_covariance[comp, 0, 1:]
        covariances[comp] = augmented_covariance[comp, 1:, 1:]
    #
    # means = np.dot(resp.T, X) / nk[:, np.newaxis]
    # covariances = {"full": _estimate_gaussian_covariances_full,
    #                "tied": _estimate_gaussian_covariances_tied,
    #                "diag": _estimate_gaussian_covariances_diag,
    #                "spherical": _estimate_gaussian_covariances_spherical
    #                }[covariance_type](resp, X, nk, means, reg_covar)
    return nk, means, covariances

class MaskedGaussianMixture(GaussianMixture):
    def __init__(self, n_components=1, covariance_type='full', tol=1e-3,
                 reg_covar=1e-6, max_iter=100, n_init=1, init_params='kmeans',
                 weights_init=None, means_init=None, precisions_init=None,
                 random_state=None, warm_start=False,
                 verbose=0, verbose_interval=10):
        super(MaskedGaussianMixture, self).__init__(
            n_components=n_components, tol=tol, reg_covar=reg_covar,
            max_iter=max_iter, n_init=n_init, init_params=init_params,
            random_state=random_state, warm_start=warm_start,
            verbose=verbose, verbose_interval=verbose_interval,
            covariance_type=covariance_type, weights_init=weights_init,
            means_init=means_init, precisions_init=precisions_init)

    def _e_step(self, X):
        """E step.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        log_prob_norm : float
            Mean of the logarithms of the probabilities of each sample in X

        log_responsibility : array, shape (n_samples, n_components)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.
        """
        # Step 1: Calculate responsibilities
        log_prob_norm = np.empty(X.shape[0], self.n_components)
        log_resp = np.empty(X.shape[0], self.n_components)

        for s in range(self.n_missing_patterns):
            s_rows = _I(s, self.row_index_of_patterns)
            s_cols = _O(s, self.feature_status_in_patterns)
            log_prob_norm[s_rows, :], log_resp[s_rows, :] = \
                self._estimate_log_prob_resp(X[s_rows, :][:, s_cols])
        # log_prob_norm, log_resp = self._estimate_log_prob_resp(X)
        resp = np.exp(log_resp)

        # Step 2: Calculate E(z * x_mis | x_obs, theta)
        T = np.empty((self.n_components, X.shape[1]+1, X.shape[1]+1))
        for comp in range(self.n_components):
            X_weighted = np.copy(X)
            X_weighted = X_weighted * resp[:, comp]

            # Setup theta
            # Covariance and augmented covariance matrix
            cov = self.covariances_[comp]
            means = self.means_[comp]

            theta = np.concatenate((means, cov), axis=0)
            theta = np.concatenate((np.insert(means, 0, 1, axis=1).T, theta),
                                   axis=1)

            # Setup Tobs, or the T-matrix with sums of the observed data and
            # also sums of their products (See Schafer, p. 215)
            vecOne = np.ones(X_weighted.shape[0]).reshape((-1, 1))

            top = np.insert(np.ma.dot(vecOne.T, X_weighted).data, 0,
                            X.shape[0], axis=1)
            bottom = np.concatenate((np.ma.dot(X_weighted.T, vecOne).data,
                                     np.ma.dot(X_weighted.T, X_weighted).data),
                                    axis=1)
            Tobs = np.concatenate((top, bottom), axis=0)

            # Calculate the matrix of expected sufficient statistics
            T[comp] = _exp_sufficient_stats(X, theta, np.copy(Tobs),
                                      self.n_missing_patterns)
            # n_weighted = np.sum(resp[:, comp])
            # T[comp] = T_exp / n_weighted
        self.sufficient_stats = T
        return np.mean(log_prob_norm), log_resp

    def _m_step(self, X, log_resp):
        """M step.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        log_resp : array-like, shape (n_samples, n_components)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.
        """
        n_samples, _ = X.shape
        self.weights_, self.means_, self.covariances_ = (
            _estimate_masked_gaussian_parameters(X, np.exp(log_resp),
                                                 self.reg_covar,
                                                 self.covariance_type,
                                                 self.sufficient_stats))
        self.weights_ /= n_samples
        self.precisions_cholesky_ = _compute_precision_cholesky(
            self.covariances_, self.covariance_type)

    def fit(self, X, y=None):
        """Estimate model parameters with the EM algorithm.

        The method fit the model `n_init` times and set the parameters with
        which the model has the largest likelihood or lower bound. Within each
        trial, the method iterates between E-step and M-step for `max_iter`
        times until the change of likelihood or lower bound is less than
        `tol`, otherwise, a `ConvergenceWarning` is raised.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        self
        """
        # Check input data and initial parameters
        X = _check_inf(_check_X(X, self.n_components, force_all_finite=False))
        self._check_initial_parameters(X)
        mask = _check_empty_axis(_get_mask(X, "NaN"))

        # if we enable warm_start, we will have a unique initialisation
        do_init = not(self.warm_start and hasattr(self, 'converged_'))
        n_init = self.n_init if do_init else 1

        max_lower_bound = -np.infty
        self.converged_ = False

        random_state = check_random_state(self.random_state)

        # Create a matrix of unique missingness patterns
        unique_anti_nan_pattern = np.unique(~mask, return_inverse=True, axis=0)
        R = unique_anti_nan_pattern[0]

        # Index for a particular missing pattern 's'
        I_all = unique_anti_nan_pattern[1]
        # To get "I" for say s=2, we do np.where(I==2)

        # Number of missingness patterns
        # S = R.shape[0]
        self.n_missing_patterns = R.shape[0]
        s_unique = np.unique(I_all)

        # X by missingness pattern
        sort_index = np.argsort(I_all)
        X = X[sort_index, :]
        mask = mask[sort_index, :]
        X = np.ma.array(X, mask=mask)
        unique_anti_nan_pattern = np.unique(~mask, return_inverse=True,
                                            axis=0)
        R = unique_anti_nan_pattern[0]
        I_all = unique_anti_nan_pattern[1]

        self.feature_status_in_patterns = unique_anti_nan_pattern[0]
        self.row_index_of_patterns = unique_anti_nan_pattern[1]
        # Functions to get patterns and column labels of specific 's'

        # Returns index of those rows belonging to pattern s
        # # Returns index of those rows belonging to pattern s
        # def I(s):
        #     return np.where(I_all == s)[0]
        #
        # # Returns observed columns in pattern s
        # def O(s):
        #     return np.where(R[s, ])[0]
        #
        # # Returns missing columns in pattern s
        # def M(s):
        #     return np.where(~R[s, ])[0]

        n_samples, _ = X.shape
        for init in range(n_init):
            self._print_verbose_msg_init_beg(init)

            if do_init:
                self._initialize_parameters(X, random_state)
                self.lower_bound_ = -np.infty

            for n_iter in range(self.max_iter):
                prev_lower_bound = self.lower_bound_
                log_prob_norm, log_resp = self._e_step(X)
                self._m_step(X, log_resp)
                self.lower_bound_ = self._compute_lower_bound(
                    log_resp, log_prob_norm)

                change = self.lower_bound_ - prev_lower_bound
                self._print_verbose_msg_iter_end(n_iter, change)

                if abs(change) < self.tol:
                    self.converged_ = True
                    break

            self._print_verbose_msg_init_end(self.lower_bound_)

            if self.lower_bound_ > max_lower_bound:
                max_lower_bound = self.lower_bound_
                best_params = self._get_parameters()
                best_n_iter = n_iter

        if not self.converged_:
            warnings.warn('Initialization %d did not converge. '
                          'Try different init parameters, '
                          'or increase max_iter, tol '
                          'or check for degenerate data.'
                          % (init + 1), ConvergenceWarning)

        self._set_parameters(best_params)
        self.n_iter_ = best_n_iter

        return self
