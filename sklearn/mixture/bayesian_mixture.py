"""Bayesian Gaussian Mixture Model."""

# TODO
# - Check if we can pass the covariance directly/precision
# - Remove the gaussian with small weights
# - check parameters XXX
# - change _alpha_prior -> alpha_prior_
# X Function in private
# - not do an EM but an ME
# - remove trace smw
# - check if this is precision or covariance
# - Modifier pour ne pas avoir deux fois le meme code.

# Author: Wei Xue <xuewei4d@gmail.com>
#         Thierry Guillemot <thierry.guillemot.work@gmail.com>

import numpy as np
from scipy import linalg
from scipy.special import digamma, gammaln

from .base import BaseMixture, _check_shape
from .gaussian_mixture import _check_precision_matrix
from .gaussian_mixture import _check_precision_positivity
from .gaussian_mixture import _compute_precision_cholesky
from .gaussian_mixture import _estimate_gaussian_parameters
from ..utils import check_array
from ..utils.validation import check_is_fitted


def _log_dirichlet_norm(alpha):
    """Estimate the log of the Dirichlet distribution normalization term.

    Parameters
    ----------
    alpha : array-like, shape (n_samples,)
        The parameters values of the Dirichlet distribution.

    Returns
    -------
    log_dirichlet_norm : float
        The log normalization of the Dirichlet distribution.
    """
    return gammaln(np.sum(alpha)) - np.sum(gammaln(alpha))


def _log_wishart_norm(nu, precision_chol, n_features):
    """Estimate the log of the Wishart distribution normalization term.

    Parameters
    ----------
    nu : float
        The parameters values of the Whishart distribution.

    precision_chol : array-like, shape (n_features, n_features)
        The Cholesky decomposition of the precision matrix.

    n_features : int
        The number of features.

    Return
    ------
    log_wishart_norm : float
        The log normalization of the Wishart distribution.
    """
    # XXX see if we remove the np.log(np.pi)
    # The determinant of the cholesky decomposition of the precision
    # is half of the determinant of the precision
    return -(nu * np.sum(np.log(np.diag(precision_chol))) +
             nu * n_features * .5 * np.log(2.) +
             n_features * (n_features - 1.) * .25 * np.log(np.pi) +
             np.sum(gammaln(.5 * (nu - np.arange(0, n_features)))))


# FIXME To remove
def _estimate_wishart_entropy(nu, precision_chol, log_lambda, n_features):
    """Estimate the entropy of the Wishart distribution.

    Parameters
    ----------
    nu : float
        The parameters values of the Whishart distribution.

    precision_chol : float
        The Cholesky decomposition of the precision matrix.

    log_lambda : array-like, shape (n_components,)
        The log values of the Wishart lambda.

    n_features : int
        The number of features.

    Return
    ------
    wishart_entropy : float
        The entropy of the Wishart distribution.
    """
    return (- _log_wishart_norm(nu, precision_chol, n_features) -
            .5 * (nu - n_features - 1.) * log_lambda + .5 * nu * n_features)


def _gamma_entropy_spherical(a, b):
    """Estimate the Gamma distribution entropy for 'spherical' precision.

    Parameters
    ----------
    a : array-like, shape (n_components,)

    b : array-like, shape (n_components,)

    Returns
    -------
    spherical_gamma_entropy : array-like, shape (n_components,)
    """
    return gammaln(a) - (a - 1.) * digamma(a) + np.log(b) + a


def _gamma_entropy_diag(a, b):
    """Estimate the Gamma distribution entropy for 'spherical' precision.

    Parameters
    ----------
    a : array-like, shape (n_components,)

    b : array-like, shape (n_components, n_features)

    Returns
    -------
    diag_gamma_entropy : array-like, shape (n_components,)
    """
    return ((gammaln(a) - (a - 1.) * digamma(a) + a) * len(b) +
            np.sum(np.log(b)))


class BayesianGaussianMixture(BaseMixture):
    """Variational estimation of a Gaussian mixture.

    Representation of a variational inference for a Bayesian Gaussian mixture
    model probability distribution. This class allows to do inference of an
    approximate posterior distribution over the parameters of a Gaussian
    mixture distribution. The number of components can be inferred from the
    data.

    Read more in the :ref:`User Guide <bgmm>`.

    Parameters
    ----------
    n_components: int, default to 1.
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

    max_iter : int, default to 100.
        The number of EM iterations to perform.

    n_init : int, default to 1.
        The number of initializations to perform. The best results is kept.

    init_params : {'kmeans', 'random'}, defaults to 'kmeans'.
        The method used to initialize the weights, the means and the
        covariances.
        Must be one of::
        'kmeans' : responsibilities are initialized using kmeans.
        'random' : responsibilities are initialized randomly.

    alpha_init : float, optional.
        The user-provided alpha prior parameter of the Dirichlet distribution.
        If is None, the alpha prior is set to `1. / n_components`.

    beta_init : float, optional.
        The user-provided beta prior parameter of the Gaussian
        distribution. If it is None, beta prior is set to 1.

    mean_init : array-like, shape (`n_features`,), optional
        The user-provided mean prior of the Gaussian distribution.
        If it is None, the mean prior is set to the mean of X.

    nu_init : float, optional.
        The user-provided nu prior parameter of the precision distribution.
        If it is None, the nu prior is set to `n_features`.

    covariance_init : float or array-like, optional
        The user-provided covariance prior of the precision distribution.
        If it is None, the covariance prior is initialized using the covariance
        of X. The shape depends on `covariance_type`::
            (`n_features`, `n_features`) if 'full',
            (`n_features`, `n_features`) if 'tied',
            (`n_features`)               if 'diag',
            float                        if 'spherical'

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
    weights_ : array-like, shape (`n_components`,)
        The weights of each mixture components.

    beta_ : array-like, shape (`n_components`, )
        The beta parameters of the Gaussian distributions of the means.

    means_ : array-like, shape (`n_components`, `n_features`)
        The mean of each mixture component.

    nu_ : array-like, shape (`n_components`,)
        The nu parameters of the precision distribution.

    covariances_ : array-like
        The covariance of each mixture component.
        The shape depends on `covariance_type`::
            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'

    precisions_ : array-like
        The precision matrices for each component in the mixture. A precision
        matrix is the inverse of a covariance matrix. A covariance matrix is
        symmetric positive definite so the mixture of Gaussian can be
        equivalently parameterized by the precision matrices. Storing the
        precision matrices instead of the covariance matrices makes it more
        efficient to compute the log-likelihood of new samples at test time.
        The shape depends on `covariance_type`::
            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'

    precisions_cholesky_ : array-like
        The cholesky decomposition of the precision matrices of each mixture
        component. A precision matrix is the inverse of a covariance matrix.
        A covariance matrix is symmetric positive definite so the mixture of
        Gaussian can be equivalently parameterized by the precision matrices.
        Storing the precision matrices instead of the covariance matrices makes
        it more efficient to compute the log-likelihood of new samples at test
        time. The shape depends on `covariance_type`::
            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'

    converged_ : bool
        True when convergence was reached in fit(), False otherwise.

    n_iter : int
        Number of step used by the best fit of EM to reach the convergence.

    See Also
    --------
    GaussianMixture : Finite Gaussian mixture fit with EM.
    """

    def __init__(self, n_components=1, covariance_type='full', tol=1e-3,
                 reg_covar=1e-6, max_iter=100, n_init=1, init_params='kmeans',
                 alpha_init=None, beta_init=None, mean_init=None,
                 nu_init=None, covariance_init=None, random_state=None,
                 warm_start=False, verbose=0, verbose_interval=20):
        super(BayesianGaussianMixture, self).__init__(
            n_components=n_components, tol=tol, reg_covar=reg_covar,
            max_iter=max_iter, n_init=n_init, init_params=init_params,
            random_state=random_state, warm_start=warm_start,
            verbose=verbose, verbose_interval=verbose_interval)

        self.covariance_type = covariance_type
        self.alpha_init = alpha_init
        self.beta_init = beta_init
        self.mean_init = mean_init
        self.nu_init = nu_init
        self.covariance_init = covariance_init

    def _check_parameters(self, X):
        """Check the VBGMM parameters are well defined.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        """
        if self.covariance_type not in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError("Invalid value for 'covariance_type': %s "
                             "'covariance_type' should be in "
                             "['spherical', 'tied', 'diag', 'full']"
                             % self.covariance_type)
        self._check_weights_parameters()
        self._check_means_parameters(X)
        self._check_precision_parameters(X)
        self._check_covariance_prior_parameter(X)

    def _check_weights_parameters(self):
        """Check the parameter of the Dirichlet distribution."""
        if self.alpha_init is None:
            self._alpha_prior = 1. / self.n_components
        elif self.alpha_init > 0.:
            self._alpha_prior = self.alpha_init
        else:
            raise ValueError("The parameter 'alpha_init' should be "
                             "greater than 0., but got %.3f."
                             % self.alpha_init)

    def _check_means_parameters(self, X):
        """Check the parameters of the Gaussian distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        """
        _, n_features = X.shape

        if self.beta_init is None:
            self._beta_prior = 1.
        elif self.beta_init > 0.:
            self._beta_prior = self.beta_init
        else:
            raise ValueError("The parameter 'beta_init' should be "
                             "greater than 0., but got %.3f."
                             % self.beta_init)

        if self.mean_init is None:
            self._mean_prior = X.mean(axis=0)
        else:
            self._mean_prior = check_array(self.mean_init,
                                           dtype=[np.float64, np.float32],
                                           ensure_2d=False)
            _check_shape(self._mean_prior, (n_features, ), 'means')

    def _check_precision_parameters(self, X):
        """Check the prior parameters of the precision distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        """
        _, n_features = X.shape

        if self.nu_init is None:
            self._nu_prior = n_features
        elif self.nu_init > n_features - 1.:
            self._nu_prior = self.nu_init
        else:
            raise ValueError("The parameter 'nu_init' "
                             "should be greater than %d, but got %.3f."
                             % (n_features - 1, self.nu_init))

    def _check_covariance_prior_parameter(self, X):
        """Check the `_covariance_prior` parameters depending of
        `covariance_type`.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
        """
        _, n_features = X.shape

        if self.covariance_init is None:
            self._covariance_prior = {
                'full': np.eye(n_features),
                'tied': np.eye(n_features),
                'diag': np.diag(np.atleast_2d(np.cov(X.T))),
                'spherical': np.var(X, axis=0).mean()
            }[self.covariance_type]

        elif self.covariance_type in ['full', 'tied']:
            self._covariance_prior = check_array(
                self.covariance_init, dtype=[np.float64, np.float32],
                ensure_2d=False)
            _check_shape(self._covariance_prior, (n_features, n_features),
                         '%s covariance_init' % self.covariance_type)
            _check_precision_matrix(self._covariance_prior,
                                    self.covariance_type)
        elif self.covariance_type is 'diag':
            self._covariance_prior = check_array(
                self.covariance_init, dtype=[np.float64, np.float32],
                ensure_2d=False)
            _check_shape(self._covariance_prior, (n_features,),
                         '%s covariance_init' % self.covariance_type)
            _check_precision_positivity(self._covariance_prior,
                                        self.covariance_type)
        # spherical case
        elif self.covariance_init > 0.:
            self._covariance_prior = self.covariance_init
        else:
            raise ValueError("The parameter 'spherical covariance_init' "
                             "should be greater than 0., but got %.3f."
                             % self.covariance_init)

    def _initialize(self, X, resp):
        """Initialization of the mixture parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        resp : array-like, shape (n_samples, n_components)
        """
        nk, xk, sk = _estimate_gaussian_parameters(X, resp, self.reg_covar,
                                                   self.covariance_type)

        self._estimate_weights(nk)
        self._estimate_means(nk, xk)
        self._estimate_precisions(nk, xk, sk)

        # XXX Remove constant term
        self._estimate_distribution_norms()

    def _estimate_weights(self, nk):
        """Estimate the parameters of the Dirichlet distribution.

        Parameters
        ----------
        nk : array-like, shape (n_components,)
        """
        self.alpha_ = self._alpha_prior + nk

    def _estimate_means(self, nk, xk):
        """Estimate the parameters of the Gaussian distribution.

        Parameters
        ----------
        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)
        """
        self.beta_ = self._beta_prior + nk
        self.means_ = (self._beta_prior * self._mean_prior +
                       nk[:, np.newaxis] * xk) / self.beta_[:, np.newaxis]

    def _estimate_precisions(self, nk, xk, sk):
        """Estimate the precisions parameters of the precision distribution.

        Parameters
        ----------
        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like
            The shape depends of `covariance_type`:
            'full' : (n_components, n_features, n_features)
            'tied' : (n_features, n_features)
            'diag' : (n_components, n_features)
            'spherical' : (n_components,)
        """
        {"full": self._estimate_wishart_full,
         "tied": self._estimate_wishart_tied,
         "diag": self._estimate_gamma_diag,
         "spherical": self._estimate_gamma_spherical
         }[self.covariance_type](nk, xk, sk)

        self.precisions_cholesky_, self.cov_chol_ = _compute_precision_cholesky(
            self.covariances_, self.covariance_type)

        if self.covariance_type is 'full':
            self.precisions_ = np.array([
                np.dot(prec_chol, prec_chol.T)
                for prec_chol in self.precisions_cholesky_])

        elif self.covariance_type is 'tied':
            self.precisions_ = np.dot(self.precisions_cholesky_,
                                      self.precisions_cholesky_.T)
        else:
            self.precisions_ = self.precisions_cholesky_ ** 2

    def _estimate_wishart_full(self, nk, xk, sk):
        """Estimate the Wishart distribution parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_components, n_features, n_features)
        """
        _, n_features = xk.shape

        # FIXME Typo error in certain version of Bishop for nu
        self.nu_ = self._nu_prior + nk

        self.covariances_ = np.empty((self.n_components, n_features,
                                      n_features))

        for k in range(self.n_components):
            diff = xk[k] - self._mean_prior
            self.covariances_[k] = (self._covariance_prior + nk[k] * sk[k] +
                                    nk[k] * self._beta_prior / self.beta_[k] *
                                    np.outer(diff, diff))

    def _estimate_wishart_tied(self, nk, xk, sk):
        """Estimate the Wishart distribution parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_features, n_features)
        """
        _, n_features = xk.shape

        # Typo error in certain version of Bishop for nu
        # FIXME We add n_features because we do a sum_n_features of 1 ???
        self.nu_ = self._nu_prior + nk.sum() / self.n_components

        diff = xk - self._mean_prior
        self.covariances_ = (self._covariance_prior +
                             sk * nk.sum() / self.n_components +
                             self._beta_prior / self.n_components *
                             (np.dot((nk / self.beta_) * diff.T, diff)))

    def _estimate_gamma_diag(self, nk, xk, sk):
        """Estimate the Gamma distribution parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_components, n_features)
        """
        _, n_features = xk.shape

        # FIXME Typo error in certain version of Bishop for nu
        self.nu_ = self._nu_prior + nk

        diff = xk - self._mean_prior
        self.covariances_ = (
            self._covariance_prior +
            .5 * (nk[:, np.newaxis] * sk +
                  (nk * self._beta_prior / self.beta_)[:, np.newaxis] *
                  np.square(diff)))

    def _estimate_gamma_spherical(self, nk, xk, sk):
        """Estimate the Gamma distribution parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_components,)
        """
        _, n_features = xk.shape

        # FIXME Typo error in certain version of Bishop for nu
        self.nu_ = self._nu_prior + nk

        diff = xk - self._mean_prior
        self.covariances_ = (self._covariance_prior + .5 / n_features *
                             (nk * sk + (nk * self._beta_prior / self.beta_) *
                              np.mean(np.square(diff), 1)))

    def _estimate_distribution_norms(self):
        """Estimate the distributions norm used to define the lowerbounds."""
        n_features, = self._mean_prior.shape

        # Distribution norm for the weights
        self._log_dirichlet_norm_prior = _log_dirichlet_norm(
            self._alpha_prior * np.ones(self.n_components))

        # Distribution norm for the means
        self._log_gaussian_norm_prior = (
            .5 * n_features * np.log(self._beta_prior / (2. * np.pi)))

        # Distribution norm for the covariance
        if self.covariance_type in ['full', 'tied']:
            # Computation of the cholesky decomposition of the precision matrix
            try:
                covariance_chol = linalg.cholesky(self._covariance_prior,
                                                  lower=True)
            except linalg.LinAlgError:
                raise ValueError("Invalid value for 'covariance_init'. The "
                                 "'covariance_init' should be a full rank.")

            precision_prior_chol = linalg.solve_triangular(
                covariance_chol, np.eye(n_features), lower=True).T
            self._log_wishart_norm_prior = (
                _log_wishart_norm(self._nu_prior, precision_prior_chol,
                                  n_features))

        elif self.covariance_type == 'diag':
            # lambda_inv_W_prior has n_feature Gamma distribution
            self._log_gamma_norm_prior = (
                self._nu_prior * np.sum(np.log(self._covariance_prior)) -
                len(self._covariance_prior) * gammaln(self._nu_prior))
        elif self.covariance_type == 'spherical':
            # lambda_inv_W_prior has only 1 Gamma distribution
            self._log_gamma_norm_prior = (
                self._nu_prior * np.log(self._covariance_prior) -
                gammaln(self._nu_prior))

    def _check_is_fitted(self):
        check_is_fitted(self, ['alpha_', 'beta_', 'means_', 'nu_',
                               'covariances_', 'precisions_',
                               'precisions_cholesky_'])

    def _m_step(self, X, log_resp, resp):
        """M step.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        resp : array-like, shape (n_samples, n_components)
        """
        n_samples, _ = X.shape

        nk, xk, sk = _estimate_gaussian_parameters(X, resp, self.reg_covar,
                                                   self.covariance_type)
        self._estimate_weights(nk)
        self._estimate_means(nk, xk)
        self._estimate_precisions(nk, xk, sk)

        return self._compute_lower_bound(log_resp, resp) / n_samples

    def _e_step(self, X):
        """E step.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        Returns
        -------
        log-responsibility : scalar

        responsibility : array, shape (n_samples, n_components)
        """
        n_features, _ = X.shape
        _, _, log_resp = self._estimate_log_prob_resp(X)
        return log_resp, np.exp(log_resp)

    def _estimate_log_weights(self):
        return digamma(self.alpha_) - digamma(np.sum(self.alpha_))

    def _estimate_log_prob(self, X):
        return {"full": self._estimate_log_prob_full,
                "tied": self._estimate_log_prob_tied,
                "diag": self._estimate_log_prob_diag,
                "spherical": self._estimate_log_prob_spherical
                }[self.covariance_type](X)

    def _estimate_log_prob_full(self, X):
        # second item in Equation 3.10
        n_samples, n_features = X.shape

        log_prob = np.empty((n_samples, self.n_components))
        self._log_lambda = np.empty((self.n_components, ))

        # SOME ERROR HERE !!!!
        for k in range(self.n_components):
            log_det_precisions = 2. * np.sum(np.log(np.diag(
                self.precisions_cholesky_[k])))
            # Equation 3.43 Corrected here XXX
            self._log_lambda[k] = (
                np.sum(digamma(.5 * (self.nu_[k] -
                                     np.arange(0, n_features)))) +
                n_features * np.log(2.) + log_det_precisions)

            # Equation 3.48
            # XXXX ERROR WITH THE TRANSPOSE !!!!!!!
            y = np.dot(X - self.means_[k], self.precisions_cholesky_[k].T)
            mahala_dist = (n_features / self.beta_[k] +
                           self.nu_[k] * np.sum(np.square(y), axis=1))

            log_prob[:, k] = .5 * (self._log_lambda[k] - mahala_dist)
        log_prob -= .5 * n_features * np.log(2. * np.pi)

        return log_prob

    def _estimate_log_prob_tied(self, X):
        # second item in Equation 3.10
        n_samples, n_features = X.shape

        log_prob = np.empty((n_samples, self.n_components))

        log_det_precisions = -2. * np.sum(np.log(np.diag(
            self.precisions_cholesky_)))

        # Equation 3.43
        self._log_lambda = (
            np.sum(digamma(.5 * (self.nu_ + 1 - np.arange(0, n_features)))) +
            n_features * np.log(2) - log_det_precisions)

        # Equation 3.48
        for k in range(self.n_components):
            y = np.dot(X - self.means_[k], self.precisions_cholesky_)
            mahala_dist = np.sum(np.square(y), axis=1)

            log_prob[:, k] = -.5 * (- self._log_lambda +
                                    n_features / self.beta_[k] +
                                    self.nu_ * mahala_dist)
        log_prob -= .5 * n_features * np.log(2. * np.pi)

        return log_prob

    def _estimate_log_prob_diag(self, X):
        _, n_features = X.shape

        # second item in Equation 4.10
        self._log_lambda = (n_features * digamma(self.nu_) +
                            np.sum(np.log(self.precisions_), axis=1))

        # Equation 4.35
        log_gauss = self.nu_ * (
            np.sum((self.means_ * self.precisions_cholesky_) ** 2, axis=1) -
            2. * np.dot(X, (self.means_ * self.precisions_).T) +
            np.dot(X ** 2, self.precisions_.T))

        return -.5 * (n_features * np.log(2. * np.pi) - self._log_lambda +
                      n_features / self.beta_ + log_gauss)

    def _estimate_log_prob_spherical(self, X):
        _, n_features = X.shape

        # second item in Equation 5.3
        self._log_lambda = n_features * (digamma(self.nu_) +
                                         np.log(self.precisions_))

        # Equation 4.35
        log_gauss = self.nu_ * self.precisions_ * (
            np.sum((self.means_ ** 2), 1) - 2. * np.dot(X, self.means_.T) +
            np.sum(X ** 2, 1)[:, np.newaxis])

        return -.5 * (n_features * np.log(2. * np.pi) - self._log_lambda +
                      n_features / self.beta_ + log_gauss)

    def _compute_lower_bound(self, log_resp, resp):
        """Estimate the lower bound of the model to check the convergence."""
        # XXX REMOVE CONSTANT TERMS
        n_features, = self._mean_prior.shape

        # Prob with that
        log_wishart = np.empty(self.n_components)
        for k, (nu, prec_chol) in enumerate(zip(self.nu_,
                                                self.precisions_cholesky_)):
            log_wishart[k] = _log_wishart_norm(nu, prec_chol, n_features)

        # We have removed the constant term
        lowe = (-np.sum(resp * log_resp) + self._log_dirichlet_norm_prior -
                _log_dirichlet_norm(self.alpha_) +
                self.n_components * self._log_wishart_norm_prior -
                np.sum(log_wishart) +
                0.5 * n_features *
                (self.n_components * np.log(self._beta_prior) -
                 np.sum(np.log(self.beta_))))


# 0.5 * n_features * self._log_wishart_norm_prior
# -0.5 * np.sum(log_wishart)
# 0.5 * n_features * self.n_components * np.log(self._beta_prior)
# -0.5 * n_features * np.sum(np.log(self.beta_))

        return lowe
        # log_q_z = np.sum(resp * log_resp)
        # log_p_Zw = np.sum(resp * log_prob)

        # # Equation 7.5
        # log_pi = self._estimate_log_weights()
        # log_p_XZ = np.sum(resp * (log_prob - log_pi)) - .5 * self._trace_sW
        # # Equation 7.6
        # # XXX (vérifier si log_pi ou log_prob)
        # # Equation 7.7
        # log_p_weight = ((self._alpha_prior - 1.) * np.sum(log_pi) +
        #                 self._log_dirichlet_norm_prior)
        # log_p_lambda = self._estimate_p_lambda()

        # # Equation 7.10
        # # Equation 7.11
        # log_q_weight = (np.sum((self.alpha_ - 1.) * log_pi) +
        #                 _log_dirichlet_norm(self.alpha_))
        # log_q_lambda = self._estimate_q_lambda()

        # return (log_p_XZ + log_p_Zw + log_p_weight + log_p_lambda -
        #         log_q_z - log_q_weight - log_q_lambda)

    def _estimate_p_lambda(self):
        return {'full': self._estimate_p_lambda_full,
                'tied': self._estimate_p_lambda_tied,
                'diag': self._estimate_p_lambda_diag,
                'spherical': self._estimate_p_lambda_spherical
                }[self.covariance_type]()

    def _estimate_p_lambda_full(self):
        n_features, = self._mean_prior.shape
        # Equation 7.9 :
        # m^t W m term
        sum_y = np.empty(self.n_components)
        for k in range(self.n_components):
            y = np.dot(self.means_[k] - self._mean_prior,
                       self.precisions_cholesky_[k])
            sum_y[k] = np.sum(np.square(y))

        # First sum of 7.9
        temp1 = (self.n_components * self._log_gaussian_norm_prior +
                 .5 * np.sum(self._log_lambda -
                             self._beta_prior * self.nu_ * sum_y -
                             n_features * self._beta_prior / self.beta_))

        # Wishart + sum(log_lambda)
        temp2 = (self.n_components * self._log_wishart_norm_prior +
                 .5 * (self._nu_prior - n_features - 1.) *
                 np.sum(self._log_lambda))

        # Final sum of 7.9
        trace_W0invW = np.empty(self.n_components)
        for k, precision in enumerate(self.precisions_):
            trace_W0invW[k] = np.sum((self._covariance_prior * precision))
        temp3 = -.5 * np.sum(self.nu_ * trace_W0invW)

        return temp1 + temp2 + temp3

    def _estimate_p_lambda_tied(self):
        n_features, = self._mean_prior.shape
        # Equation 7.9 :
        # (m^t W m) term
        sum_y = np.empty(self.n_components)
        for k in range(self.n_components):
            y = np.dot(self.means_[k] - self._mean_prior,
                       self.precisions_cholesky_)
            sum_y[k] = np.sum(np.square(y))

        # First sum of 7.9
        temp1 = (self.n_components * self._log_gaussian_norm_prior +
                 .5 * (self.n_components * self._log_lambda -
                       np.sum(self._beta_prior * self.nu_ * sum_y +
                              n_features * self._beta_prior / self.beta_)))

        # Wishart + sum(log_lambda)
        temp2 = (self.n_components * self._log_wishart_norm_prior +
                 .5 * (self._nu_prior - n_features - 1.) *
                 self.n_components * self._log_lambda)

        # Final sum of 7.9
        trace_W0inv_W = np.sum(self._covariance_prior * self.precisions_)
        temp3 = -.5 * self.n_components * self.nu_ * trace_W0inv_W

        return temp1 + temp2 + temp3

    def _estimate_p_lambda_diag(self):
        n_features, = self._mean_prior.shape
        # Equation 7.16 :
        # (m^t b m) term
        sum_y = np.sum(np.square(self.means_ - self._mean_prior) *
                       self.precisions_, axis=1)

        # First sum of 7.16
        temp1 = (self.n_components * self._log_gaussian_norm_prior +
                 .5 * np.sum(self._log_lambda - self._beta_prior *
                             (n_features / self.beta_ + self.nu_ * sum_y)))

        # Gaussian + sum(log_lambda) of 7.16
        temp2 = (self.n_components * self._log_gamma_norm_prior +
                 (self._nu_prior - 1.) * np.sum(self._log_lambda))

        # Final sum of 7.16
        temp3 = -np.sum(self.nu_ * np.sum(self._covariance_prior *
                                          self.precisions_, axis=1))
        return temp1 + temp2 + temp3

    def _estimate_p_lambda_spherical(self):
        n_features, = self._mean_prior.shape
        # Equation 7.16 :
        # (m^t W m) term
        sum_y = self.precisions_ * np.sum(np.square(self.means_ -
                                                    self._mean_prior), axis=1)

        # First sum of 7.16
        temp1 = (self.n_components * self._log_gaussian_norm_prior +
                 .5 * np.sum(self._log_lambda - self._beta_prior *
                             (n_features / self.beta_ + self.nu_ * sum_y)))

        # Wishart + sum(log_lambda) of 7.16
        temp2 = (self.n_components * self._log_gamma_norm_prior +
                 (self._nu_prior - 1.) * np.sum(self._log_lambda))

        # Final sum of 7.16
        temp3 = np.sum(- self.nu_ * self._covariance_prior * self.precisions_)

        return temp1 + temp2 + temp3

    def _estimate_q_lambda(self):
        return {'full': self._estimate_q_lambda_full,
                'tied': self._estimate_q_lambda_tied,
                'diag': self._estimate_q_lambda_diag,
                'spherical': self._estimate_q_lambda_spherical
                }[self.covariance_type]()

    def _estimate_q_lambda_full(self):
        # Equation 7.14
        n_features, = self._mean_prior.shape
        wishart_entropy = np.empty(self.n_components)
        for k in range(self.n_components):
            wishart_entropy[k] = _estimate_wishart_entropy(
                self.nu_[k], self.precisions_cholesky_[k],
                self._log_lambda[k], n_features)
        return np.sum(.5 * self._log_lambda +
                      .5 * n_features * np.log(self.beta_ / (2. * np.pi)) -
                      .5 * n_features - wishart_entropy)

    def _estimate_q_lambda_tied(self):
        # Equation 7.14
        n_features, = self._mean_prior.shape
        wishart_entropy = _estimate_wishart_entropy(
            self.nu_, self.precisions_cholesky_, self._log_lambda, n_features)
        return (.5 * self.n_components * self._log_lambda +
                .5 * n_features * np.sum(np.log(self.beta_ / (2. * np.pi))) -
                .5 * n_features * self.n_components -
                self.n_components * wishart_entropy)

    def _estimate_q_lambda_diag(self):
        # Equation 7.16
        n_features, = self._mean_prior.shape
        return np.sum(
            .5 * self._log_lambda +
            .5 * n_features * (np.log(self.beta_ / (2. * np.pi)) - 1.) -
            _gamma_entropy_diag(self.nu_, self.precisions_))

    def _estimate_q_lambda_spherical(self):
        # Equation 7.16
        n_features = self._mean_prior.shape[0]
        return np.sum(
            .5 * self._log_lambda +
            .5 * n_features * (np.log(self.beta_ / (2. * np.pi)) - 1.) -
            n_features * _gamma_entropy_spherical(self.nu_, self.precisions_))

    def _get_parameters(self):
        return (self.alpha_, self.beta_, self.means_, self.nu_,
                self.covariances_, self.precisions_, self.precisions_cholesky_)

    def _set_parameters(self, params):
        (self.alpha_, self.beta_, self.means_, self.nu_, self.covariances_,
         self.precisions_, self.precisions_cholesky_) = params

        # Attributes computation
        self. weights_ = self.alpha_ / np.sum(self.alpha_)
        if self.covariance_type is 'full':
            self.covariances_ /= self.nu_[:, np.newaxis, np.newaxis]
        elif self.covariance_type is 'diag':
            self.covariances_ /= self.nu_[:, np.newaxis]
        else:
            self.covariances_ /= self.nu_
