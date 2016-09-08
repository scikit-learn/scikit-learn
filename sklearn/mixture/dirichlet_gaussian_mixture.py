"""Dirichlet Process Bayesian Gaussian Mixture Model."""
# Author: Wei Xue <xuewei4d@gmail.com>
#         Thierry Guillemot <thierry.guillemot.work@gmail.com>
# License: BSD 3 clause

import numpy as np
from scipy.special import digamma, betaln

from ..utils.validation import check_is_fitted
from .gaussian_mixture import _compute_log_det_cholesky
from .bayesian_mixture import BayesianGaussianMixture, _log_wishart_norm


class DirichletGaussianMixture(BayesianGaussianMixture):
    """Dirichlet Process Bayesian Gaussian Mixture.

    An infinite mixture model with the Dirichlet Process as a prior
    distribution on the number of clusters. In practice the approximate
    inference algorithm uses a truncated distribution with a fixed maximum
    number of components (called the Stick-breaking representation), but almost
    always the number of components actually used depends on the data.

    This class allows for easy and efficient inference of an approximate
    posterior distribution over the parameters of a Gaussian mixture model with
    a variable number of components (smaller than the truncation parameter
    n_components).

    Read more in the :ref:`User Guide <dpgmm>`.

    Parameters
    ----------
    n_components: int, defaults to 1.
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
        lower bound average gain on the likelihood (of the training data with
        respect to the model) is below this threshold.

    reg_covar : float, defaults to 1e-6.
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

    weight_concentration_prior : float, default to 1.
        The user-provided concentration prior parameter of the
        Beta distribution. The higher concentration puts more mass in the
        center and will lead to more components being active, while a lower
        concentration parameter will lead to more mass at the edge of the
        simplex. The value of the parameter must be greater than 0.

    mean_precision_prior : float, optional.
        The user-provided mean precision prior parameter of the Gaussian
        distribution. Controls the extend to where means can be placed. Smaller
        values concentrate the means of each clusters around `mean_prior`.
        The value of the parameter must be greater than 0.
        If it is None, mean precision prior is set to 1.

    mean_prior : array-like, shape (`n_features`,), optional
        The user-provided mean prior of the Gaussian distribution.
        If it is None, the mean prior is set to the mean of X.

    degrees_of_freedom_prior : float, optional.
        The user-provided number of degrees of freedom prior parameter of the
        covariance distribution. If it is None, the prior of the number of
        degrees of freedom is set to `n_features`.

    covariance_prior : float or array-like, optional
        The user-provided covariance prior of the covariance distribution.
        If it is None, the emiprical covariance prior is initialized using the
        covariance of X. The shape depends on `covariance_type`::
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

    verbose_interval : int, default to 10.
        Number of iteration done before the next print.

    Attributes
    ----------
    weights_ : array-like, shape (`n_components`,)
        The weights of each mixture components.

    means_ : array-like, shape (`n_components`, `n_features`)
        The mean of each mixture component.

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

    n_iter_ : int
        Number of step used by the best fit of EM to reach the convergence.

    lower_bound_ : float
        Lower bound value on the likelihood (of the training data with
        respect to the model) of the best fit of EM.

    beta_concentration_alpha_ : array-like, shape (`n_components`, )
        The alpha concentration parameters of the Beta distribution of each
        component.

    beta_concentration_beta_ : array-like, shape (`n_components`, )
        The beta concentration parameters of the Beta distribution of each
        component.

    mean_precision_prior : float
        The mean precision prior parameters of the Gaussian distributions of
        the means used during the fit process. Controls the extend to where
        means can be placed. Smaller values concentrate the means of each
        clusters around `mean_prior`.

    mean_precision_ : array-like, shape (`n_components`, )
        The mean precision parameters of the Gaussian distributions of
        the means.

    means_prior_ : array-like, shape (`n_features`,)
        The mean prior of each mixture component.

    freedom_degrees_prior_ : float
        The prior of the number of degrees of freedom parameters of the
        covariance distribution.

    freedom_degrees_ : array-like, shape (`n_components`,)
        The number of degrees of freedom parameters of the covariance
        distribution.

    covariance_prior_ : float or array-like
        The covariance prior of the covariance distribution.
        The shape depends on `covariance_type`::
            (`n_features`, `n_features`) if 'full',
            (`n_features`, `n_features`) if 'tied',
            (`n_features`)               if 'diag',
            float                        if 'spherical'

    See Also
    --------
    GaussianMixture : Finite Gaussian mixture fit with EM.

    BayesianGaussianMixture : Finite gaussian mixture model fit with a
        variational algorithm.
    """

    def __init__(self, n_components=1, covariance_type='full', tol=1e-3,
                 reg_covar=1e-6, max_iter=100, n_init=1, init_params='kmeans',
                 weight_concentration_prior=1,
                 mean_precision_prior=None, mean_prior=None,
                 degrees_of_freedom_prior=None, covariance_prior=None,
                 random_state=None, warm_start=False, verbose=0,
                 verbose_interval=10):
        super(DirichletGaussianMixture, self).__init__(
            n_components=n_components, covariance_type=covariance_type,
            tol=tol, reg_covar=reg_covar, max_iter=max_iter, n_init=n_init,
            init_params=init_params, weight_concentration_prior=None,
            mean_precision_prior=mean_precision_prior, mean_prior=mean_prior,
            degrees_of_freedom_prior=degrees_of_freedom_prior,
            covariance_prior=covariance_prior, random_state=random_state,
            warm_start=warm_start, verbose=verbose,
            verbose_interval=verbose_interval)
        self.weight_concentration_prior = weight_concentration_prior

    def _check_weights_parameters(self):
        """Check the parameter of the Beta distribution."""
        if self.weight_concentration_prior <= 0:
            raise ValueError("The parameter 'weight_concentration_prior' "
                             "should be greater than 0., but got %.3f."
                             % self.weight_concentration_prior)

    def _estimate_weights(self, nk):
        """Estimate the parameters of the Dirichlet distribution.

        Parameters
        ----------
        nk : array-like, shape (n_components,)
        """
        self.beta_concentration_alpha_ = 1. + nk
        self.beta_concentration_beta_ = (
            self.weight_concentration_prior +
            np.hstack((np.cumsum(nk[::-1])[-2::-1], 0)))

    def _estimate_log_weights(self):
        digamma_sum = digamma(self.beta_concentration_alpha_ +
                              self.beta_concentration_beta_)
        digamma_a = digamma(self.beta_concentration_alpha_)
        digamma_b = digamma(self.beta_concentration_beta_)
        return (digamma_a - digamma_sum +
                np.hstack((0, np.cumsum(digamma_b - digamma_sum)[:-1])))

    def _compute_lower_bound(self, log_resp, log_prob_norm):
        """Estimate the lower bound of the model.

        The lower bound is used to detect the convergence and has to decrease
        at each iteration.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        log_resp : array, shape (n_samples, n_components)
            Logarithm of the posterior probabilities (or responsibilities) of
            the point of each sample in X.

        log_prob_norm : float
            Logarithm of the probability of each sample in X.

        Returns
        -------
        lower_bound : float
        """
        # Contrary to the original formula, we have done some simplification
        # and removed all the constant terms.
        n_features, = self.mean_prior_.shape

        # We removed `.5 * n_features * np.log(self.degrees_of_freedom_)`
        # because the precision matrix is normalized.
        log_det_precisions_chol = (_compute_log_det_cholesky(
            self.precisions_cholesky_, self.covariance_type, n_features) -
            .5 * n_features * np.log(self.degrees_of_freedom_))

        if self.covariance_type == 'tied':
            log_wishart = self.n_components * np.float64(_log_wishart_norm(
                self.degrees_of_freedom_, log_det_precisions_chol, n_features))
        else:
            log_wishart = np.sum(_log_wishart_norm(
                self.degrees_of_freedom_, log_det_precisions_chol, n_features))

        return -(np.sum(np.exp(log_resp) * log_resp) + log_wishart -
                 np.sum(betaln(self.beta_concentration_alpha_,
                               self.beta_concentration_beta_)) +
                 0.5 * n_features * np.sum(np.log(self.mean_precision_)))

    def _check_is_fitted(self):
        check_is_fitted(self, ['beta_concentration_alpha_',
                               'beta_concentration_beta_',
                               'mean_precision_',
                               'means_', 'degrees_of_freedom_',
                               'covariances_', 'precisions_cholesky_'])

    def _get_parameters(self):
        return (self.beta_concentration_alpha_,
                self.beta_concentration_beta_,
                self.mean_precision_, self.means_,
                self.degrees_of_freedom_, self.covariances_,
                self.precisions_cholesky_)

    def _set_parameters(self, params):
        (self.beta_concentration_alpha_,
         self.beta_concentration_beta_, self.mean_precision_,
         self.means_, self.degrees_of_freedom_, self.covariances_,
         self.precisions_cholesky_) = params

        weight_dirichlet_sum = (self.beta_concentration_alpha_ +
                                self.beta_concentration_beta_)
        tmp = (self.beta_concentration_beta_ / weight_dirichlet_sum)
        self.weights_ = (
            self.beta_concentration_alpha_ / weight_dirichlet_sum *
            np.hstack((1, np.cumprod(tmp[:-1]))))
        self.weights_ /= np.sum(self.weights_)

        if self.covariance_type == 'full':
            self.precisions_ = np.array([
                np.dot(prec_chol, prec_chol.T)
                for prec_chol in self.precisions_cholesky_])

        elif self.covariance_type == 'tied':
            self.precisions_ = np.dot(self.precisions_cholesky_,
                                      self.precisions_cholesky_.T)
        else:
            self.precisions_ = self.precisions_cholesky_ ** 2
