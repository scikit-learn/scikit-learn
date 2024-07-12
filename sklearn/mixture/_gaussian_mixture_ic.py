"""GaussianMixtureIC"""

# Authors: Tingshan Liu <tliu68@jhmi.edu>
#          Thomas Athey <tathey1@jhmi.edu>
#          Benjamin Pedigo <bpedigo@jhu.edu>


import numpy as np

from ..base import BaseEstimator, ClusterMixin
from ..model_selection import GridSearchCV
from ..utils._param_validation import (
    Integral,
    Interval,
    InvalidParameterError,
    StrOptions,
)
from ..utils.validation import check_is_fitted
from . import GaussianMixture


def _check_multi_comp_inputs(input, name, default):
    if isinstance(input, (np.ndarray, list)):
        input = list(np.unique(input))
    elif isinstance(input, str):
        if input not in default:
            raise InvalidParameterError(
                f"The '{name}' parameter of GaussianMixtureIC must be one of {default}."
                f" Got {input} instead."
            )
        if input != "all":
            input = [input]
        else:
            input = default.copy()
            input.remove("all")
    else:
        raise InvalidParameterError(
            f"The '{name}' parameter of GaussianMixtureIC must be one of {default}. "
            f"Got {input} instead."
        )
    return input


class GaussianMixtureIC(ClusterMixin, BaseEstimator):
    """Gaussian mixture with BIC/AIC.

    Automatic Gaussian Mixture Model (GMM) selection via the
    Bayesian Information Criterion (BIC)
    or the Akaike Information Criterion (AIC).

    Such criteria are useful to select the value
    of the gaussian mixture parameters by making a trade-off
    between the goodness of fit and the complexity of the model.

    Parameters
    ----------
    min_components : int, default=2
        The minimum number of mixture components to consider.
        If ``max_components`` is not None, ``min_components`` must be
        less than or equal to ``max_components``.

    max_components : int or None, default=10
        The maximum number of mixture components to consider.
        Must be greater than or equal to ``min_components``.

    covariance_type : {'full' (default), 'tied', 'diag', 'spherical', 'all'},
            optional
        String or list/array describing the type of covariance parameters
        to use.
        If a string, it must be one of:

        - 'full'
            each component has its own general covariance matrix
        - 'tied'
            all components share the same general covariance matrix
        - 'diag'
            each component has its own diagonal covariance matrix
        - 'spherical'
            each component has its own single variance
        - 'all'
            considers all covariance structures in
            ['spherical', 'diag', 'tied', 'full']

        If a list/array, it must be a list/array of strings containing only
        'spherical', 'tied', 'diag', and/or 'spherical'.

    n_init : int, optional (default = 1)
        The number of initializations to perform.

    init_params : {'kmeans' (default), 'k-means++', 'random', 'random_from_data'}
        The method used to initialize the weights, the means and the precisions
        for Gaussian mixture modeling.

    criterion : str {"bic" or "aic"}, optional, (default = "bic")
        Select the best model based on Bayesian Information Criterion (bic) or
        Aikake Information Criterion (aic).

    n_jobs : int
        The number of jobs to use for the computation
        This works by computing each of the n_init runs in parallel.

    tol : float, default=1e-3
        The convergence threshold. EM iterations will stop when the
        lower bound average gain is below this threshold.

    reg_covar : float, default=1e-6
        Non-negative regularization added to the diagonal of covariance.
        Allows to assure that the covariance matrices are all positive.

    weights_init : array-like of shape (n_components, ), default=None
        The user-provided initial weights.
        If it is None, weights are initialized using the `init_params` method.

    means_init : array-like of shape (n_components, n_features), default=None
        The user-provided initial means,
        If it is None, means are initialized using the `init_params` method.

    precisions_init : array-like, default=None
        The user-provided initial precisions (inverse of the covariance
        matrices).
        If it is None, precisions are initialized using the 'init_params'
        method.
        The shape depends on 'covariance_type'::

            (n_components,)                        if 'spherical',
            (n_features, n_features)               if 'tied',
            (n_components, n_features)             if 'diag',
            (n_components, n_features, n_features) if 'full'

    random_state : int, RandomState instance or None, default=None
        Controls the random seed given to the method chosen to initialize the
        parameters (see `init_params`).
        In addition, it controls the generation of random samples from the
        fitted distribution (see the method `sample`).
        Pass an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    warm_start : bool, default=False
        If 'warm_start' is True, the solution of the last fitting is used as
        initialization for the next call of fit(). This can speed up
        convergence when fit is called several times on similar problems.
        In that case, 'n_init' is ignored and only a single initialization
        occurs upon the first call.
        See :term:`the Glossary <warm_start>`.

    max_iter : int, optional (default = 100)
        The maximum number of EM iterations to perform.

    verbose : int, default=0
        Enable verbose output. If 1 then it prints the current
        initialization and each iteration step. If greater than 1 then
        it prints also the log probability and the time needed
        for each step.

    verbose_interval : int, default=10
        Number of iteration done before the next print.

    Attributes
    ----------
    criterion_ : array-like
        The value of the information criteria ('aic', 'bic') across all
        numbers of components. The number of component which has the smallest
        information criterion is chosen.

    n_components_ : int
        Number of clusters for the model with the best bic/aic.

    covariance_type_ : str
        Covariance type for the model with the best bic/aic.

    best_estimator_ : :class:`sklearn.mixture.GaussianMixture`
        Object with the best bic/aic.

    weights_ : array-like of shape (n_components,)
        The weights of each mixture components for the model with the best bic/aic.

    means_ : array-like of shape (n_components, n_features)
        The mean of each mixture component for the model with the best bic/aic.

    covariances_ : array-like
        The covariance of each mixture component for the model with the best bic/aic.
        The shape depends on `covariance_type_`. See
        :class:`~sklearn.mixture.GaussianMixture` for details.

    precisions_ : array-like
        The precision matrices for each component in the mixture for the model with
        the best bic/aic. See :class:`~sklearn.mixture.GaussianMixture` for details.

    precisions_cholesky_ : array-like
        The cholesky decomposition of the precision matrices of each mixture component
        for the model with the best bic/aic.
        See :class:`~sklearn.mixture.GaussianMixture` for details.

    converged_ : bool
        True only when convergence was reached in :term:`fit` for the model
        with the best bic/aic, False otherwise.

    n_iter_ : int
        Number of step used by the best fit of EM for the best model
        to reach the convergence.

    lower_bound_ : float
        Lower bound value on the log-likelihood (of the training data with
        respect to the model) of the best fit of EM.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point.

    See Also
    --------
    GaussianMixture : Fit Gaussian mixture model.
    BayesianGaussianMixture : Gaussian mixture model fit with a variational
        inference.

    Notes
    -----
    This algorithm was strongly inspired by mclust [3]_,
    a clustering package for R.

    References
    ----------
    .. [1] `Fraley, C., & Raftery, A. E. (2002). Model-based clustering,
        discriminant analysis, and density estimation.
        Journal of the American statistical Association, 97(458), 611-631.
        <https://doi.org/10.1198/016214502760047131>_`

    .. [2] `Athey, T. L., Pedigo, B. D., Liu, T., & Vogelstein, J. T. (2019).
        AutoGMM: Automatic and Hierarchical Gaussian Mixture Modeling
        in Python. arXiv preprint arXiv:1909.02688.
        <https://arxiv.org/abs/1909.02688>_`

    .. [3] `Scrucca, L., Fop, M., Murphy, T. B., & Raftery, A. E. (2016).
        mclust 5: Clustering, Classification and Density Estimation Using
        Gaussian Finite Mixture Models. The R journal, 8(1), 289-317.
        <https://doi.org/10.32614/RJ-2016-021>_`

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.mixture import GaussianMixtureIC
    >>> X = np.array([[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]])
    >>> gmIC = GaussianMixtureIC(max_components=4)
    >>> print(np.sort(gmIC.fit_predict(X)))
    [0 0 0 1 1 1]
    >>> print(gmIC.n_components_)
    2
    """

    _parameter_constraints: dict = {
        **GaussianMixture._parameter_constraints,
        "criterion": [StrOptions({"aic", "bic"})],
        "min_components": [Interval(Integral, 1, None, closed="left")],
        "max_components": [Interval(Integral, 1, None, closed="left")],
        "n_jobs": [Integral, None],
        "covariance_type": [
            StrOptions({"spherical", "diag", "tied", "full", "all"}),
            list,
            np.ndarray,
        ],
    }
    _parameter_constraints.pop("n_components")

    def __init__(
        self,
        *,
        min_components=2,
        max_components=10,
        covariance_type="full",
        n_init=1,
        init_params="kmeans",
        criterion="bic",
        n_jobs=None,
        tol=1e-3,
        reg_covar=1e-6,
        weights_init=None,
        means_init=None,
        precisions_init=None,
        random_state=None,
        warm_start=False,
        max_iter=100,
        verbose=0,
        verbose_interval=10,
    ):
        super().__init__()
        self.covariance_type = covariance_type
        self.min_components = min_components
        self.max_components = max_components
        self.criterion = criterion
        self.n_jobs = n_jobs
        self.n_init = n_init
        self.init_params = init_params
        self.tol = tol
        self.reg_covar = reg_covar
        self.weights_init = weights_init
        self.means_init = means_init
        self.precisions_init = precisions_init
        self.random_state = random_state
        self.warm_start = warm_start
        self.max_iter = max_iter
        self.verbose = verbose
        self.verbose_interval = verbose_interval

    def _check_parameters(self):
        covariance_type = _check_multi_comp_inputs(
            self.covariance_type,
            "covariance_type",
            ["spherical", "diag", "tied", "full", "all"],
        )

        return covariance_type

    def criterion_score(self, estimator, X):
        """Callable to pass to GridSearchCV that will use the BIC score.

        Parameters
        ----------
        estimator : estimator object
            A score function to calculate either BIC or AIC.

        X : array-like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        score : float
            The BIC or AIC score.
        """
        if self.criterion == "bic":
            return -estimator.bic(X)
        else:
            return -estimator.aic(X)

    def fit(self, X, y=None):
        """Fit several Gaussian mixture models to the data.

        Initialize with agglomerative clustering then
        estimate model parameters with EM algorithm.
        Select the best model according to the chosen
        information criterion.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        self._validate_params()
        covariance_type = self._check_parameters()
        X = self._validate_data(X, dtype=[np.float64, np.float32], ensure_min_samples=1)

        # check n_components against sample size
        if self.max_components > X.shape[0]:
            msg = "max_components must be <= n_samples, but max_components"
            msg += "= {}, n_samples = {}".format(self.max_components, X.shape[0])
            raise ValueError(msg)

        # Ensure reproducibility
        if self.random_state is not None:
            np.random.seed(self.random_state)

        param_grid = {
            "covariance_type": covariance_type,
            "n_components": range(self.min_components, self.max_components + 1),
        }

        grid_search = GridSearchCV(
            GaussianMixture(
                init_params=self.init_params, max_iter=self.max_iter, n_init=self.n_init
            ),
            param_grid=param_grid,
            scoring=self.criterion_score,
        )
        grid_search.fit(X)

        self.criterion_ = -grid_search.cv_results_["mean_test_score"]
        self.n_components_ = grid_search.best_params_["n_components"]
        self.covariance_type_ = grid_search.best_params_["covariance_type"]

        best_estimator = grid_search.best_estimator_
        self.best_estimator_ = best_estimator
        self.weights_ = best_estimator.weights_
        self.means_ = best_estimator.means_
        self.covariances_ = best_estimator.covariances_
        self.precisions_ = best_estimator.precisions_
        self.precisions_cholesky_ = best_estimator.precisions_cholesky_
        self.converged_ = best_estimator.converged_
        self.n_iter_ = best_estimator.n_iter_
        self.lower_bound_ = best_estimator.lower_bound_
        self.n_features_in_ = X.shape[1]
        self.labels_ = best_estimator.predict(X)

        return self

    def predict(self, X):
        """Predict clusters based on the best Gaussian mixture model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        labels : array, shape (n_samples,)
            Component labels.
        """
        check_is_fitted(self, ["best_estimator_"], all_or_any=all)
        X = self._validate_data(X, reset=False)
        labels = self.best_estimator_.predict(X)

        return labels

    def fit_predict(self, X, y=None):
        """Fit the models and predict clusters based on the best model.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        labels : array, shape (n_samples,)
            Component labels.
        """
        self.fit(X, y)

        labels = self.predict(X)
        return labels
