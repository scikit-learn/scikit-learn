import numpy as np
from scipy.stats import multivariate_normal

from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.mixture._base import BaseMixture
from sklearn.utils._array_api import _logsumexp
from sklearn.mixture._gaussian_mixture import _compute_precision_cholesky

class MaskedGaussianMixture(GaussianMixture):

    def __init__(
        self,
        n_components=1,
        *,
        covariance_type="full",
        tol=1e-3,
        reg_covar=1e-6,
        max_iter=100,
        n_init=1,
        init_params="kmeans",
        weights_init=None,
        means_init=None,
        precisions_init=None,
        random_state=None,
        warm_start=False,
        verbose=0,
        verbose_interval=10,
    ):
        super().__init__(
            n_components=n_components,
            tol=tol,
            reg_covar=reg_covar,
            max_iter=max_iter,
            n_init=n_init,
            init_params=init_params,
            covariance_type=covariance_type,
            weights_init = weights_init,
            means_init = means_init,
            precisions_init = precisions_init,
            random_state=random_state,
            warm_start=warm_start,
            verbose=verbose,
            verbose_interval=verbose_interval,
        )
        self.mean_cond = None
        self.cov_cond = None
        self.means_ = means_init
        self.covariances_ = None
        self.weights_ = None
        self.precisions_cholesky_ = None

    def _initialize_parameters(self, X, random_state, xp=None):
        n_samples, n_features = X.shape
        valid_rows = X[~np.isnan(X).any(axis=1)]
        idx = np.random.choice(len(valid_rows), self.n_components, replace=False)
        self.means_ = valid_rows[idx].copy()

        self.covariances_ = np.zeros((self.n_components, n_features, n_features))
        for i in range(self.n_components):
            self.covariances_[i] = np.identity(n_features)
        self.weights_ = np.ones(self.n_components) / self.n_components
        self.precisions_cholesky_ = _compute_precision_cholesky(
            self.covariances_, self.covariance_type, xp=xp
        )

    # def _e_step(self, X, xp=None):
    #     log_weights = self._estimate_log_weights()
    #     log_prob = self._estimate_log_prob(X)
    #     weighted_log_prob = log_prob + log_weights
    #     log_prob_norm = _logsumexp(weighted_log_prob, axis=1)
    #     log_resp = weighted_log_prob - log_prob_norm[:, None]
    #     self.compute_conditionals(X)
    #     return log_prob_norm.mean() , log_resp


    def _m_step(self, X, log_resp, xp=None):
        n_samples, n_features = X.shape
        resp = np.exp(log_resp)
        n_k = np.sum(resp, axis=0) + 1e-5
        for k in range(self.n_components):
          self.means_[k] = np.sum(resp[:, k][:, None] * self.mean_cond[:, k, :], axis=0)  / n_k[k]

        for k in range(self.n_components):
            cov = np.zeros((n_features, n_features))
            for i in range(n_samples):
                x_mean_cond = self.mean_cond[i, k, :]
                x_cov_cond = self.cov_cond[i, k, :, :] # shape (n, d, d)
                diff = x_mean_cond - self.means_[k]
                cov += resp[i, k] * (np.outer(diff, diff) + x_cov_cond)
            cov /= (n_k[k])
            cov += self.reg_covar * np.eye(n_features)
            self.covariances_[k] = cov

        self.weights_ = n_k / n_samples

        self.precisions_cholesky_ = _compute_precision_cholesky(
            self.covariances_, self.covariance_type, xp=xp
        )

    def compute_conditionals(self, X):
        n_samples, n_features = X.shape
        mask = ~np.isnan(X)
        self.mean_cond = np.zeros((n_samples, self.n_components, n_features))
        self.cov_cond = np.zeros((n_samples, self.n_components, n_features, n_features))
        for i in range(n_samples):
            for k in range(self.n_components):
                obs_idx = mask[i]
                miss_idx = ~mask[i]
                sigma_oo = self.covariances_[k][np.ix_(obs_idx, obs_idx)]
                sigma_om = self.covariances_[k][np.ix_(obs_idx, miss_idx)]
                sigma_mo = self.covariances_[k][np.ix_(miss_idx, obs_idx)]
                sigma_mm = self.covariances_[k][np.ix_(miss_idx, miss_idx)]
                mean_obs = self.means_[k, obs_idx]
                mean_miss = self.means_[k, miss_idx]
                x_obs = X[i, obs_idx]

                if miss_idx.any():
                    #x_mean = (mean_miss + sigma_mo @ np.linalg.inv(sigma_oo) @ (x_obs - mean_obs))
                    x_mean = (mean_miss + sigma_mo @ np.linalg.solve(sigma_oo , (x_obs - mean_obs)))
                    x_full = X[i].copy()
                    x_full[miss_idx] = x_mean
                    self.mean_cond[i, k] = x_full
                    cov_full = np.zeros((n_features, n_features))
                    #cov_full[np.ix_(miss_idx, miss_idx)] = sigma_mm - sigma_mo @ np.linalg.inv(sigma_oo) @ sigma_om
                    cov_full[np.ix_(miss_idx, miss_idx)] = sigma_mm - sigma_mo @ np.linalg.solve(sigma_oo, sigma_om)
                    self.cov_cond[i, k] = cov_full
                else:
                    self.mean_cond[i, k] = X[i].copy()
                    self.cov_cond[i, k] = np.zeros((n_features, n_features))


    def _estimate_log_prob(self, X, xp=None):
        n_samples, n_features = X.shape
        self.compute_conditionals(X)
        mask = ~np.isnan(X)
        probs = np.zeros((n_samples, self.n_components))
        for i in range(n_samples):
            obs_idx = mask[i]
            for k in range(self.n_components):
                x_obs = X[i, obs_idx]
                mean = self.means_[k, obs_idx]
                sigma_oo = self.covariances_[k][np.ix_(obs_idx, obs_idx)]
                sigma_oo = sigma_oo + self.reg_covar * np.eye(obs_idx.sum())  # stabilize covariance
                probs[i, k] = multivariate_normal.logpdf(x_obs, mean, sigma_oo, allow_singular=True)
        return probs

    def _get_parameters(self):
        return self.weights_, self.means_, self.covariances_, self.precisions_cholesky_

    def _set_parameters(self, params, xp=None):
        self.weights_, self.means_, self.covariances_, self.precisions_cholesky_ = params
