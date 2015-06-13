"""Dirichlet Process Bayesian Gaussian Mixture Model."""
import numpy as np
from scipy.special import digamma, gammaln

from ..utils.validation import check_is_fitted
from .bayesianmixture import BayesianGaussianMixture, check_shape


def _log_beta_norm(a, b):
    """The log of the normalization of Beta distribution."""
    return gammaln(a+b) - gammaln(a) - gammaln(b)


class DirichletProcessGaussianMixture(BayesianGaussianMixture):

    """Dirichlet Process Bayesian Gaussian Mixture Model.

    An infinite mixture model with the Dirichlet Process as a prior
    distribution on the number of clusters. In practice the
    approximate inference algorithm uses a truncated distribution with
    a fixed maximum number of components, but almost always the number
    of components actually used depends on the data.

    Stick-breaking Representation of a Gaussian mixture model
    probability distribution. This class allows for easy and efficient
    inference of an approximate posterior distribution over the
    parameters of a Gaussian mixture model with a variable number of
    components (smaller than the truncation parameter n_components).

    Read more in the :ref:`User Guide <dpgmm>`.

    Parameters
    ----------
    n_components: int, default 1.
        Number of mixture components.

    precision_type: string, default to 'full'.
        String describing the type of covariance parameters to
        use.  Must be one of 'spherical', 'tied', 'diag', 'full'.

    tol : float, default 1e-6
        Convergence threshold.

    reg_covar : float, defaults to 0.
        Non-negative regularization to the diagonal of covariance.

    max_iter : int, default to 100.
        Maximum number of iterations to perform before convergence.

    n_init : int, default to 1.
        Number of initializations to perform. The best results is kept.

    params : string, default to None.
        Controls which parameters are updated in the initialization
        process.  Can contain any combination of 'w' for weights,
        'm' for means, and 'c' for covars.  Defaults to 'wmc'.

    init_params : string, defaults to 'kmeans'.
        Controls how parameters are initialized unless the prior
        parameters are provided by users. It should be one of 'kmeans',
        'random', None. If it is not None, the responsibilities and model
        parameters are initialized by the chosen method.

    gamma_prior: float, default to 0.1.
        Real number representing the concentration parameter of
        the dirichlet process. Intuitively, the Dirichlet Process
        is as likely to start a new cluster for a point as it is
        to add that point to a cluster with gamma_prior elements. A
        higher alpha means more clusters, as the expected number
        of clusters is ``gamma_prior*log(N)``.

    m_prior : array-like, shape (1, `n_components`), defaults to None.
        User-provided prior parameter of 'm_0'. If it None, m_prior_
        is set to the mean of X.

    beta_prior : float, defaults to None.
        User-provided prior parameter of 'beta_0'. If it is None,
        beta_prior_ is set to 1.

    nu_prior : float, defaults to None.
        User-provided prior parameter of 'nu_0'. If it is None,
        nu_prior is set to `n_features` for 'full' and 'tied' precision,
        and set to 0.5 for 'diag' and 'spherical' precision.

    inv_W_prior : array-like, defaults to None.
        User-provided prior parameter of 'inv(W_0)'. If it is None,
        inv_W_prior is set to cov(X) * n_features for 'full' and 'tied'
        precision, and set to '.5 * np.diag(np.cov(X.T, bias=1))'

    random_state: RandomState or an int seed, defaults to None.
        A random number generator instance.

    verbose : int, default 0
        Controls output verbosity.

    Attributes
    ----------
    covariance_type : string
        String describing the type of covariance parameters used by
        the DP-GMM.  Must be one of 'spherical', 'tied', 'diag', 'full'.

    n_components : int
        Number of mixture components.

    weights_ : array, shape (`n_components`,)
        Mixing weights for each mixture component.

    means_ : array, shape (`n_components`, `n_features`)
        Mean parameters for each mixture component.

    precs_ : array
        Precision (inverse covariance) parameters for each mixture
        component.  The shape depends on `covariance_type`::

            (`n_components`, 'n_features')                if 'spherical',
            (`n_features`, `n_features`)                  if 'tied',
            (`n_components`, `n_features`)                if 'diag',
            (`n_components`, `n_features`, `n_features`)  if 'full'

    converged_ : bool
        True when convergence was reached in fit(), False otherwise.

    See Also
    --------
    GMM : Finite Gaussian mixture model fit with EM

    VBGMM : Finite Gaussian mixture model fit with a variational
        algorithm, better for situations where there might be too little
        data to get a good estimate of the covariance matrix.
    """

    def __init__(self, n_components=1, precision_type='full',
                 tol=1e-6, reg_covar=0,
                 max_iter=100, n_init=1, params=None, init_params='kmeans',
                 gamma_prior=0.1,
                 m_prior=None, beta_prior=None,
                 nu_prior=None, inv_W_prior=None,
                 random_state=None, verbose=0, verbose_interval=10):
        super(DirichletProcessGaussianMixture, self).__init__(
            n_components=n_components, precision_type=precision_type, tol=tol,
            reg_covar=reg_covar, max_iter=max_iter, n_init=n_init,
            init_params=init_params, alpha_prior=None, m_prior=m_prior,
            beta_prior=beta_prior, nu_prior=nu_prior, inv_W_prior=inv_W_prior,
            random_state=random_state, verbose=verbose,
            verbose_interval=verbose_interval)
        self.gamma_prior = gamma_prior
        self.weight_gamma_a_ = None
        self.weight_gamma_b_ = None

    def _check_weight_prior(self, weight_gamma_prior, desired_shape):
        """Validate the prior parameter  of the weight Beta distribution.

        Parameters
        ----------
        weight_gamma_prior : float

        desired_shape : tuple

        Returns
        -------
        weight_gamma_prior : float
        """
        check_shape(weight_gamma_prior, desired_shape, 'gamma_prior')
        if weight_gamma_prior <= 0:
            raise ValueError("The parameter 'gamma_prior' should be "
                             "greater than 0, but got %.3f ."
                             % weight_gamma_prior)
        return weight_gamma_prior

    def _initialize_weight(self, nk):
        """Initialize the prior parameter of weight gamma distribution."""
        self.weight_gamma_a_, self.weight_gamma_b_ = self._estimate_weights(nk)

    def _initialize_weight_norm_prior(self):
        self._log_beta_norm_gamma_prior = _log_beta_norm(
            1, self.gamma_prior)

    def _estimate_weights(self, nk):
        tmp = np.hstack((np.cumsum(nk[::-1])[-2::-1], 0))
        return 1 + nk, self.gamma_prior + tmp

    def _m_step(self, X, resp):
        nk, xk, Sk = self._estimate_suffstat(X, resp)
        self.weight_gamma_a_, self.weight_gamma_b_ = self._estimate_weights(nk)
        self.mu_beta_, self.mu_m_ = self._estimate_mu(nk, xk)
        self.lambda_nu_, self.lambda_inv_W_ = self._estimate_lambda(nk, xk, Sk)

    def _estimate_log_weights(self):
        """Equation 9.22."""
        digamma_sum = digamma(self.weight_gamma_a_ + self.weight_gamma_b_)
        digamma_a = digamma(self.weight_gamma_a_)
        digamma_b = digamma(self.weight_gamma_b_)
        self._log_pi_a = digamma_a - digamma_sum
        self._log_pi_b = digamma_b - digamma_sum
        return (digamma_a - digamma_sum +
                np.hstack((0, np.cumsum(digamma_b - digamma_sum)[:-1])))

    def _estimate_p_weight(self):
        """Equation 9.31."""
        # c = np.arange(self.n_components - 1, 0, -1)
        # tmp = self._log_beta_norm_gamma_prior + \
        #     (self.weight_gamma_prior - 1)*self._log_pi_a
        return (self.n_components * self._log_beta_norm_gamma_prior +
                (self.gamma_prior - 1) * np.sum(self._log_pi_b))
        # + np.sum(c * tmp[:-1])

    def _estimate_q_weight(self):
        """Equation 9.32."""
        log_norm = np.empty(self.n_components)
        for k in range(self.n_components):
            log_norm[k] = _log_beta_norm(self.weight_gamma_a_[k],
                                         self.weight_gamma_b_[k])
        # c = np.arange(self.n_components-1, 0, -1)

        tmp1 = (log_norm + (self.weight_gamma_a_ - 1) * self._log_pi_a +
                (self.weight_gamma_b_ - 1) * self._log_pi_b)
        # tmp2 = log_norm + (self.weight_gamma_b_ - 1) * self._log_pi_a + \
        #     (self.weight_gamma_a_ - 1) * self._log_pi_b
        return np.sum(tmp1)
        # + np.sum(c * tmp2[:-1])

    def _check_is_fitted(self):
        check_is_fitted(self, 'weight_gamma_a_')
        check_is_fitted(self, 'weight_gamma_b_')
        check_is_fitted(self, 'mu_beta_')
        check_is_fitted(self, 'mu_m_')
        check_is_fitted(self, 'lambda_nu_')
        check_is_fitted(self, 'lambda_inv_W_')

    def _get_parameters(self):
        return (self.weight_gamma_a_, self.weight_gamma_b_, self.mu_beta_,
                self.mu_m_, self.lambda_nu_, self.lambda_inv_W_)

    def _set_parameters(self, params):
        (self.weight_gamma_a_, self.weight_gamma_b_, self.mu_beta_, self.mu_m_,
         self.lambda_nu_, self.lambda_inv_W_) = params

    @property
    def weights_(self):
        """The expected weight of each component in the mixture."""
        tmp1 = (self.weight_gamma_a_ /
                (self.weight_gamma_a_ + self.weight_gamma_b_))
        tmp2 = (self.weight_gamma_b_ /
                (self.weight_gamma_a_ + self.weight_gamma_b_))
        return tmp1 * np.hstack((1, np.cumprod(tmp2[:-1])))
