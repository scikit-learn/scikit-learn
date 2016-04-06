"""Bayesian Gaussian Mixture Model."""

# Author: Wei Xue <xuewei4d@gmail.com>
#         Thierry Guillemot <thierry.guillemot.work@gmail.com>

import numpy as np
from scipy import linalg
from scipy.special import digamma, gammaln

from .base import BaseMixture, _check_shape
from .gaussian_mixture import _check_covariance_matrix
from .gaussian_mixture import _check_covariance_positivity
from .gaussian_mixture import _estimate_gaussian_parameters
from ..utils import check_array
from ..utils.validation import check_is_fitted


def log_dirichlet_norm(alpha):
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


def log_wishart_norm(nu, precision_chol, n_features):
    """Estimate the log of the Wishart distribution normalization term.

    Parameters
    ----------
    nu : float
        The parameters values of the Whishart distribution.

    precision_chol : array-like, shape (n_features, n_features)
        The Cholesky decomposition of the inverse precision inverse matrix.

    n_features : int
        The number of features.

    Return
    ------
    log_wishart_norm : float
        The log normalization of the Wishart distribution.
    """
    return (nu * np.sum(np.log(np.diag(precision_chol))) -
            nu * n_features * .5 * np.log(2.) -
            n_features * (n_features - 1.) * .25 * np.log(np.pi) -
            np.sum(gammaln(.5 * (nu + 1. - np.arange(1, n_features + 1.)))))


def estimate_wishart_entropy(nu, precision_chol, log_lambda, n_features):
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
    return (- log_wishart_norm(nu, precision_chol, n_features) -
            .5 * (nu - n_features - 1.) * log_lambda + .5 * nu * n_features)


def gamma_entropy_spherical(a, inv_b):
    """Estimate the entropy of the Gamma distribution with 'diag' precision.

    Parameters
    ----------
    a : array-like, shape (n_components,)

    inv_b : array-like, shape (n_components,)

    Returns
    -------
    spherical_gamma_entropy : array-like, shape (n_components,)
    """
    return gammaln(a) - (a - 1.) * digamma(a) - np.log(inv_b) + a


def gamma_entropy_diag(a, inv_b):
    """The entropy of the Gamma distribution with 'diag' precision.

    Parameters
    ----------
    a : array-like, shape (n_components,)

    inv_b : array-like, shape (n_components, n_features)

    Returns
    -------
    diag_gamma_entropy : array-like, shape (n_components,)
    """
    print(a.shape, inv_b.shape)
    return ((gammaln(a) - (a - 1.) * digamma(a) + a) * len(inv_b) -
            np.sum(np.log(inv_b)))


class BayesianGaussianMixture(BaseMixture):
    """Variational estimation of a Gaussian mixture.

    Variational inference for a Bayesian Gaussian mixture model probability
    distribution. This class allows for easy and efficient inference
    of an approximate posterior distribution over the parameters of a
    Gaussian mixture model. The number of components can be inferred from the
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

    tol : float, defaults to 1e-6.
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

    alpha_prior_init : float, optional.
        The user-provided alpha prior parameter of the Dirichlet distribution.
        If is None, `_alpha_prior` is set to 1. / n_components.

    beta_prior_init : float, optional.
        The user-provided beta prior parameter of the Gaussian
        distribution. If it is None, `_beta_prior` is set to 1.

    m_prior_init : array-like, shape (`n_features`,), optional
        The user-provided mean prior of the Gaussian distribution.
        If it is None, `_m_prior` is set to the mean of X.

    nu_prior_init : float, optional.
        The user-provided nu prior parameter of the precision distribution.
        If it is None, `_nu_prior` is set to `n_features`.

    precision_prior_init : float or array-like, optional
        The user-provided precision prior of the precision distribution.
        If it is None, `precision_prior` is initialized using the covariance of
        X. The shape depends on `covariance_type`::
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
        `weights_` will not exist before a call to fit.

    means_ : array-like, shape (`n_components`, `n_features`)
        The mean of each mixture component.
        `means_` will not exist before a call to fit.

    covariances_ : array-like
        The covariance of each mixture component.
        The shape depends on `covariance_type`::
            (`n_components`,)                            if 'spherical',
            (`n_features`, `n_features`)                 if 'tied',
            (`n_components`, `n_features`)               if 'diag',
            (`n_components`, `n_features`, `n_features`) if 'full'
        `covariances_` will not exist before a call to fit.

    alpha_ : array-like, shape (`n_components`,)
        The alpha parameters of the Dirichlet distribution (weights).
        `alpha_` will not exist before a call to fit.

    beta_ : array-like, shape (`n_components`, )
        The beta parameters of the Gaussian distribution (means).
        `beta_` will not exist before a call to fit.

    m_ : array-like, shape (`n_components`, `n_features`)
        The means parameter of the Gaussian distribution (means).
        `m_` will not exist before a call to fit.

    nu_ : array-like, shape (`n_components`,)
        The nu parameters of the precision distribution.
        `nu_` will not exist before a call to fit.

    precisions_ : array-like,
        The precisions parameters of the precision distribution
        The shape depends on `precision_type`::
            (`n_components`,)                            if 'spherical',
            (`n_features`, `n_features`)                 if 'tied',
            (`n_components`, `n_features`)               if 'diag',
            (`n_components`, `n_features`, `n_features`) if 'full'
        `precisions_` will not exist before a call to fit.

    converged_ : bool
        True when convergence was reached in fit(), False otherwise.
        `converged_` will not exist before a call to fit.

    best_n_iter : int
        Number of step used by the best fit of EM to reach the convergence.
        `best_n_iter` will not exist before a call to fit.

    See Also
    --------
    GaussianMixture : Finite Gaussian mixture fit with EM.
    """

    def __init__(self, n_components=1, covariance_type='full', tol=1e-6,
                 reg_covar=0, max_iter=100, n_init=1, init_params='kmeans',
                 alpha_prior_init=None, beta_prior_init=None,
                 m_prior_init=None, nu_prior_init=None,
                 precision_prior_init=None, random_state=None,
                 warm_start=False, verbose=0, verbose_interval=10):
        super(BayesianGaussianMixture, self).__init__(
            n_components=n_components, tol=tol, reg_covar=reg_covar,
            max_iter=max_iter, n_init=n_init, init_params=init_params,
            random_state=random_state, warm_start=warm_start,
            verbose=verbose, verbose_interval=verbose_interval)

        self.covariance_type = covariance_type
        self.alpha_prior_init = alpha_prior_init
        self.beta_prior_init = beta_prior_init
        self.m_prior_init = m_prior_init
        self.nu_prior_init = nu_prior_init
        self.precision_prior_init = precision_prior_init

    def _check_parameters(self, X):
        """Check the Gaussian mixture parameters are well defined."""
        if self.covariance_type not in ['spherical', 'tied', 'diag', 'full']:
            raise ValueError("Invalid value for 'covariance_type': %s "
                             "'covariance_type' should be in "
                             "['spherical', 'tied', 'diag', 'full']"
                             % self.covariance_type)

    def _initialize(self, X, resp):
        """Initialization of the mixture parameters.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        resp : array-like, shape (n_samples, n_components)
        """
        nk, xk, sk = _estimate_gaussian_parameters(X, resp, self.reg_covar,
                                                   self.covariance_type)

        self._initialize_weights(nk)
        self._initialize_means(X, nk, xk)
        self._initialize_precisions(X, nk, xk, sk)
        self._estimate_distribution_norms()

    def _initialize_weights(self, nk):
        """Initialize the prior parameter of the Dirichlet distribution.

        Parameters
        ----------
        nk : array-like, shape (n_components,)
        """
        if self.alpha_prior_init is None:
            self._alpha_prior = 1. / self.n_components
        elif self.alpha_prior_init <= 0.:
            raise ValueError("The parameter 'alpha_prior_init' should be "
                             "greater than 0., but got %.3f."
                             % self.alpha_prior_init)
        else:
            self._alpha_prior = self.alpha_prior_init
        self._estimate_weights(nk)

    def _estimate_weights(self, nk):
        """Estimate the parameters of the Dirichlet distribution.

        Parameters
        ----------
        nk : array-like, shape (n_components,)
        """
        self.alpha_ = self._alpha_prior + nk

    def _initialize_means(self, X, nk, xk):
        """Initialize the prior parameters of the Gaussian distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)
        """
        n_features = X.shape[1]

        if self.beta_prior_init is None:
            self._beta_prior = 1.
        elif self.beta_prior_init <= 0.:
            raise ValueError("The parameter 'beta_prior_init' should be "
                             "greater than 0., but got %.3f."
                             % self.beta_prior_init)
        else:
            self._beta_prior = self.beta_prior_init

        if self.m_prior_init is None:
            self._m_prior = X.mean(axis=0)
        else:
            self._m_prior = check_array(self.m_prior_init,
                                        dtype=[np.float64, np.float32],
                                        ensure_2d=False)
            _check_shape(self.m_prior_init, (n_features, ), 'means')
        self._estimate_means(nk, xk)

    def _estimate_means(self, nk, xk):
        """Estimate the parameters of the Gaussian distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)
        """
        self.beta_ = self._beta_prior + nk
        self.m_ = (self._beta_prior * self._m_prior +
                   nk[:, np.newaxis] * xk) / self.beta_[:, np.newaxis]

    def _initialize_precisions(self, X, nk, xk, sk):
        """Initialize the prior parameters of the precision distribution.

        The precision distribution is the Wishart distribution for 'full' or
        'tied' covariance models, or the Gamma distribution for 'diagonal' and
        'spherical' covariance models.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)
        """
        n_features = X.shape[1]

        if self.nu_prior_init is None:
            self._nu_prior = n_features
        elif self.nu_prior_init <= n_features - 1.:
            raise ValueError("The parameter 'nu_prior_init' "
                             "should be greater than %d, but got %.3f."
                             % (n_features - 1, self.nu_prior_init))
        else:
            self._nu_prior = self.nu_prior_init

        self._initialize_precision_prior(X)

        self._estimate_precisions(nk, xk, sk)

    def _initialize_precision_prior(self, X):
        """Initialize `_precision_prior` depending of `covariance_type`.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like
            The shape depends of `covariance_type`:
            'full' : (n_components, n_features, n_features)
            'tied' : (n_features, n_features)
            'diag' : (n_components, n_features)
            'spherical' : (n_components,)
        """
        n_features = X.shape[1]

        if self.precision_prior_init is None:
            self._precision_prior = {
                'full': np.eye(X.shape[1]),
                'tied': np.eye(X.shape[1]),
                'diag': .5 * np.diag(np.atleast_2d(np.cov(X.T, bias=1))),
                'spherical': .5 * np.var(X, axis=0).mean()
            }[self.covariance_type]

        elif self.covariance_type in ['full', 'tied']:
            self._precision_prior = check_array(self.precision_prior_init,
                                                dtype=[np.float64, np.float32],
                                                ensure_2d=False)
            _check_shape(self._precision_prior, (n_features, n_features),
                         '%s precision_prior_init' % self.covariance_type)
            _check_covariance_matrix(self._precision_prior,
                                     self.covariance_type)
        elif self.covariance_type is 'diag':
            self._precision_prior = check_array(self.precision_prior_init,
                                                dtype=[np.float64, np.float32],
                                                ensure_2d=False)
            _check_shape(self._precision_prior, (n_features,),
                         '%s precision_prior_init' % self.covariance_type)
            _check_covariance_positivity(self._precision_prior,
                                         self.covariance_type)
        # spherical case
        elif self.precision_prior_init <= 0.:
            raise ValueError("The parameter 'spherical precision_prior_init' "
                             "should be greater than 0., but got %.3f."
                             % self.precision_prior_init)
        else:
            self._precision_prior = self.precision_prior_init

    def _estimate_precisions(self, nk, xk, Sk):
        """Estimate the precisions parameters of the precision distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like
            The shape depends of `covariance_type`:
            'full' : (n_components, n_features, n_features)
            'tied' : (n_features, n_features)
            'diag' : (n_components, n_features)
            'spherical' : (n_components,)
        """
        return {"full": self._estimate_precisions_full,
                "tied": self._estimate_precisions_tied,
                "diag": self._estimate_precisions_diag,
                "spherical": self._estimate_precisions_spherical
                }[self.covariance_type](nk, xk, Sk)

    def _estimate_precisions_full(self, nk, xk, Sk):
        """Estimate the precisions parameters of the Wishart distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_components, n_features, n_features)
        """
        n_features = xk.shape[1]

        self.nu_ = self._nu_prior + nk

        self.precisions_ = np.empty((self.n_components,
                                     n_features, n_features))
        for k in range(self.n_components):
            diff = xk[k] - self._m_prior
            self.precisions_[k] = (self._precision_prior + nk[k] * Sk[k] +
                                   (nk[k] * self._beta_prior / self.beta_[k]) *
                                   np.outer(diff, diff))

    def _estimate_precisions_tied(self, nk, xk, Sk):
        """Estimate the precisions parameters of the Wishart distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_features, n_features)
        """
        self.nu_ = self._nu_prior + nk.sum() / self.n_components

        diff = xk - self._m_prior
        self.precisions_ = (
            self._precision_prior +
            Sk * nk.sum() / self.n_components +
            self._beta_prior / self.n_components *
            (np.dot((nk / self.beta_) * diff.T, diff)))

    def _estimate_precisions_diag(self, nk, xk, Sk):
        """Estimate the precisions parameters of the Gamma distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_components, n_features)
        """
        self.nu_ = self._nu_prior + .5 * nk

        diff = xk - self._m_prior
        self.precisions_ = (
            self._precision_prior +
            .5 * (nk[:, np.newaxis] * Sk +
                  (nk * self._beta_prior / self.beta_)[:, np.newaxis] *
                  np.square(diff)))

    def _estimate_precisions_spherical(self, nk, xk, Sk):
        """Estimate the precisions parameters of the Gamma distribution.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        nk : array-like, shape (n_components,)

        xk : array-like, shape (n_components, n_features)

        sk : array-like, shape (n_components,)
        """
        n_features = xk.shape[1]

        self.nu_ = self._nu_prior + .5 * nk

        diff = xk - self._m_prior
        self.precisions_ = (
            self._precision_prior + .5 / n_features *
            (nk * Sk + (nk * self._beta_prior / self.beta_) *
             np.mean(np.square(diff), 1)))

    def _estimate_distribution_norms(self):
        """Estimate the distributions norm used to define the lowerbounds."""
        n_features = self._m_prior.shape[0]

        alpha = self._alpha_prior * np.ones(self.n_components)
        self._log_dirichlet_norm_prior = log_dirichlet_norm(alpha)

        self._log_gaussian_norm_prior = (
            .5 * n_features * np.log(self._beta_prior / (2. * np.pi)))

        if self.covariance_type in ['full', 'tied']:
            self._log_wishart_norm_prior = (
                log_wishart_norm(self._nu_prior, self._precision_prior,
                                 n_features))
        elif self.covariance_type == 'diag':
            # lambda_inv_W_prior has n_feature Gamma distribution
            self._log_gamma_norm_prior = (
                self._nu_prior * np.sum(np.log(self._precision_prior)) -
                len(self._precision_prior) * gammaln(self._nu_prior))
        elif self.covariance_type == 'spherical':
            # lambda_inv_W_prior has only 1 Gamma distribution
            self._log_gamma_norm_prior = (
                self._nu_prior * np.log(self._precision_prior) -
                gammaln(self._nu_prior))

    def _check_is_fitted(self):
        check_is_fitted(self, ['alpha_', 'beta_', 'm_', 'nu_', 'precisions_'])

    def _m_step(self, X, resp):
        nk, xk, sk = _estimate_gaussian_parameters(X, resp, self.reg_covar,
                                                   self.covariance_type)

        self._estimate_weights(nk)
        self._estimate_means(nk, xk)
        self._estimate_precisions(nk, xk, sk)

    def _e_step(self, X):
        _, log_prob, log_resp = self._estimate_log_prob_resp(X)
        resp = np.exp(log_resp)
        self._lower_bound = self._estimate_lower_bound(log_prob, resp,
                                                       log_resp)
        return self._lower_bound, resp

    def _estimate_log_weights(self):
        # save the value for computing the lower bound
        self._log_pi = digamma(self.alpha_) - digamma(np.sum(self.alpha_))
        return self._log_pi

    def _estimate_log_prob(self, X):
        return {"full": self._estimate_log_prob_full,
                "tied": self._estimate_log_prob_tied,
                "diag": self._estimate_log_prob_diag,
                "spherical": self._estimate_log_prob_spherical
                }[self.covariance_type](X)

    def _estimate_log_prob_full(self, X):
        # second item in Equation 3.8
        n_samples, n_features = X.shape[0], X.shape[1]

        ln_W_digamma = np.arange(1, n_features + 1)
        log_prob = np.empty((n_samples, self.n_components))
        log_lambda = np.empty((self.n_components, ))
        precs_chol = np.empty(self.precisions_.shape)

        for k in range(self.n_components):
            try:
                precs_chol[k] = linalg.cholesky(self.precisions_[k],
                                                lower=True)
            except linalg.LinAlgError:
                raise ValueError("'precisions_' must be symmetric, "
                                 "positive-definite")

            log_det_precisions = 2 * np.sum(np.log(np.diagonal(precs_chol[k])))
            log_lambda[k] = (
                np.sum(digamma(.5 * (self.nu_[k] + 1. - ln_W_digamma))) +
                n_features * np.log(2.) - log_det_precisions)

            W_sol = linalg.solve_triangular(precs_chol[k],
                                            (X - self.m_[k]).T,
                                            lower=True).T
            mahala_dist = np.sum(np.square(W_sol), axis=1)

            log_prob[:, k] = (
                -.5 * (log_lambda[k] + n_features / self.beta_[k] +
                       self.nu_[k] * mahala_dist))
        log_prob -= .5 * n_features * np.log(2. * np.pi)

        # save precision_chol, log_lambda for computing lower bound
        self._precisions_chol = precs_chol
        self._log_lambda = log_lambda
        return log_prob

    def _estimate_log_prob_tied(self, X):
        n_samples, n_features = X.shape[0], X.shape[1]

        ln_W_digamma = np.arange(1, n_features + 1)
        log_prob = np.empty((n_samples, self.n_components))

        try:
            precs_chol = linalg.cholesky(self.precisions_, lower=True)
        except linalg.LinAlgError:
            raise ValueError("'precisions_' must be symmetric, "
                             "positive-definite")
        log_det_precisions = 2. * np.sum(np.log(np.diagonal(precs_chol)))
        log_lambda = (
            np.sum(digamma(.5 * (self.nu_ + 1. - ln_W_digamma))) +
            n_features * np.log(2) - log_det_precisions)

        for k in range(self.n_components):
            W_sol = linalg.solve_triangular(precs_chol, (X - self.m_[k]).T,
                                            lower=True).T
            mahala_dist = np.sum(np.square(W_sol), axis=1)

            log_prob[:, k] = -.5 * (- log_lambda +
                                    n_features / self.beta_[k] +
                                    self.nu_ * mahala_dist)
        log_prob -= .5 * n_features * np.log(2. * np.pi)

        # save precision_chol, log_lambda for computing lower bound
        self._precisions_chol = precs_chol
        self._log_lambda = log_lambda
        return log_prob

    def _estimate_log_prob_diag(self, X):
        n_features = X.shape[1]
        log_lambda = (n_features * digamma(self.nu_) -
                      np.sum(np.log(self.precisions_), axis=1))

        log_gauss = self.nu_ * (np.sum((self.m_ ** 2 / self.precisions_), 1) -
                                2 * np.dot(X, (self.m_ / self.precisions_).T) +
                                np.dot(X ** 2, (1. / self.precisions_).T))
        log_prob = -.5 * (n_features * np.log(2. * np.pi) - log_lambda +
                          n_features / self.beta_ + log_gauss)

        # save log_lambda for computing lower bound
        self._log_lambda = log_lambda
        return log_prob

    def _estimate_log_prob_spherical(self, X):
        n_features = X.shape[1]

        log_lambda = n_features * (digamma(self.nu_) -
                                   np.log(self.precisions_))

        log_gauss = self.nu_ / self.precisions_ * (
            np.sum((self.m_ ** 2), 1) - 2 * np.dot(X, self.m_.T) +
            np.sum(X ** 2, 1)[:, np.newaxis])
        log_prob = -.5 * (n_features * np.log(2. * np.pi) - log_lambda +
                          n_features / self.beta_ + log_gauss)

        # save log_lambda for computing lower bound
        self._log_lambda = log_lambda
        return log_prob

    def _estimate_lower_bound(self, log_prob, resp, log_resp):
        """Estimate the lower bound of the model to check the convergence."""
        # Equation 7.5, 7.6
        log_p_XZ = np.sum(log_prob * resp)
        # Equation 7.7
        log_p_weight = (self._log_dirichlet_norm_prior +
                        (self._alpha_prior - 1.) * np.sum(self._log_pi))
        log_p_lambda = self._estimate_p_lambda()

        # Equation 7.10
        log_q_z = np.sum(resp * log_resp)
        # Equation 7.11
        log_q_weight = (np.sum((self.alpha_ - 1.) * self._log_pi) +
                        log_dirichlet_norm(self.alpha_))
        log_q_lambda = self._estimate_q_lambda()

        return (log_p_XZ + log_p_weight + log_p_lambda -
                log_q_z - log_q_weight - log_q_lambda)

    def _estimate_p_lambda(self):
        return {'full': self._estimate_p_lambda_full,
                'tied': self._estimate_p_lambda_tied,
                'diag': self._estimate_p_lambda_diag,
                'spherical': self._estimate_p_lambda_spherical
                }[self.covariance_type]()

    def _estimate_p_lambda_full(self):
        n_features = self._m_prior.shape[0]
        # Equation 7.9
        temp1 = (self._log_lambda -
                 n_features * self._beta_prior / self.beta_)

        mk_sol = np.empty(self.n_components)
        for k in range(self.n_components):
            sol = linalg.solve_triangular(
                self._precisions_chol[k],
                (self.m_[k] - self._m_prior).T, lower=True).T
            mk_sol[k] = np.sum(np.square(sol))
        temp1 -= self._beta_prior * self.nu_ * mk_sol

        # XXX @xuewei4d I think you forgot n_component in your code ?
        temp1 = (.5 * np.sum(temp1) +
                 self.n_components * self._log_gaussian_norm_prior)

        temp2 = (self.n_components * self._log_wishart_norm_prior +
                 .5 * (self._nu_prior - n_features - 1.) *
                 np.sum(self._log_lambda))

        trace_W0inv_Wk = np.empty(self.n_components)
        for k in range(self.n_components):
            chol_sol = linalg.inv(self._precisions_chol[k])
            trace_W0inv_Wk[k] = np.sum(self._precision_prior.T *
                                       np.dot(chol_sol.T, chol_sol))
        temp3 = -.5 * np.sum(self.nu_ * trace_W0inv_Wk)
        return temp1 + temp2 + temp3

    def _estimate_p_lambda_tied(self):
        n_features = self._m_prior.shape[0]

        temp1 = -n_features * self._beta_prior / self.beta_

        mk_sol = np.empty(self.n_components)
        for k in range(self.n_components):
            sol = linalg.solve_triangular(
                self._precisions_chol,
                (self.m_[k] - self._m_prior).T, lower=True).T
            mk_sol[k] = np.sum(np.square(sol))
        temp1 -= self._beta_prior * self.nu_ * mk_sol

        # XXX @xuewei4d I think you forgot n_component in your code ?
        temp1 = (.5 * (np.sum(temp1) + self.n_components * self._log_lambda) +
                 self.n_components * self._log_gaussian_norm_prior)

        temp2 = (self.n_components * self._log_wishart_norm_prior +
                 .5 * (self._nu_prior - n_features - 1.) *
                 self.n_components * self._log_lambda)

        chol_sol = linalg.inv(self._precisions_chol)
        trace_W0inv_W = np.sum(self._precision_prior.T *
                               np.dot(chol_sol.T, chol_sol))
        temp3 = -.5 * self.n_components * self.nu_ * trace_W0inv_W
        return temp1 + temp2 + temp3

    def _estimate_p_lambda_diag(self):
        n_features = self._m_prior.shape[0]

        temp1 = (self._log_lambda -
                 n_features * self._beta_prior / self.beta_)

        diff = self.m_ - self._m_prior
        temp1 -= (self._beta_prior * self.nu_ *
                  np.sum(np.square(diff) / self.precisions_, axis=1))

        # XXX @xuewei4d I think you forgot n_component in your code ?
        temp1 = (.5 * np.sum(temp1) +
                 self.n_components * self._log_gaussian_norm_prior)

        temp2 = (self.n_components * self._log_gamma_norm_prior +
                 (self._nu_prior - 1.) * np.sum(self._log_lambda))
        temp3 = (np.sum(- self.nu_ * np.sum(self._precision_prior /
                        self.precisions_, axis=1)))
        return temp1 + temp2 + temp3

    def _estimate_p_lambda_spherical(self):
        n_features = self._m_prior.shape[0]

        temp1 = (self._log_lambda -
                 n_features * self._beta_prior / self.beta_)

        diff = self.m_ - self._m_prior
        temp1 -= (self._beta_prior * self.nu_ / self.precisions_ *
                  np.sum(np.square(diff), axis=1))

        # XXX @xuewei4d I think you forgot n_component in your code ?
        temp1 = (.5 * np.sum(temp1) +
                 self.n_components * self._log_gaussian_norm_prior)

        temp2 = (self.n_components * self._log_gamma_norm_prior +
                 (self._nu_prior - 1.) * np.sum(self._log_lambda))

        temp3 = np.sum(- self.nu_ * self._precision_prior /
                       self.precisions_)
        return temp1 + temp2 + temp3

    def _estimate_q_lambda(self):
        return {'full': self._estimate_q_lambda_full,
                'tied': self._estimate_q_lambda_tied,
                'diag': self._estimate_q_lambda_diag,
                'spherical': self._estimate_q_lambda_spherical
                }[self.covariance_type]()

    def _estimate_q_lambda_full(self):
        n_features = self._m_prior.shape[0]
        wishart_entropy = np.empty(self.n_components)
        for k in range(self.n_components):
            wishart_entropy[k] = estimate_wishart_entropy(
                self.nu_[k], self._precisions_chol[k],
                self._log_lambda[k], n_features)
        return np.sum(.5 * self._log_lambda +
                      .5 * n_features * np.log(self.beta_ / (2. * np.pi)) -
                      .5 * n_features - wishart_entropy)

    def _estimate_q_lambda_tied(self):
        n_features = self._m_prior.shape[0]
        wishart_entropy = estimate_wishart_entropy(
            self.nu_, self._precisions_chol, self._log_lambda, n_features)
        return (.5 * self.n_components * self._log_lambda +
                .5 * n_features * np.sum(np.log(self.beta_ / (2. * np.pi))) -
                .5 * n_features * self.n_components -
                self.n_components * wishart_entropy)

    def _estimate_q_lambda_diag(self):
        n_features = self._m_prior.shape[0]
        return np.sum(
            .5 * self._log_lambda +
            .5 * n_features * (np.log(self.beta_ / (2. * np.pi)) - 1.) -
            gamma_entropy_diag(self.nu_, self.precisions_))

    def _estimate_q_lambda_spherical(self):
        n_features = self._m_prior.shape[0]
        return np.sum(
            .5 * self._log_lambda +
            .5 * n_features * (np.log(self.beta_ / (2. * np.pi)) - 1.) -
            n_features * gamma_entropy_spherical(self.nu_, self.precisions_))

    def _get_parameters(self):
        return (self.alpha_, self.beta_, self.m_, self.nu_, self.precisions_)

    def _set_parameters(self, params):
        (self.alpha_, self.beta_, self.m_, self.nu_, self.precisions_) = params

    @property
    def weights_(self):
        """The expected weight of each component in the mixture."""
        self._check_is_fitted()
        return self.alpha_ / np.sum(self.alpha_)

    @property
    def means_(self):
        """The expected mean of each component in the mixture."""
        self._check_is_fitted()
        return self.m_

    @property
    def covariances_(self):
        """The expected covariance of each component in the mixture."""
        self._check_is_fitted()

        if self.covariance_type is 'full':
            covars = self.precisions_ / self.nu_[:, np.newaxis, np.newaxis]
        elif self.covariance_type is 'diag':
            covars = self.precisions_ / self.nu_[:, np.newaxis]
        else:
            covars = self.precisions_ / self.nu_
        return covars
