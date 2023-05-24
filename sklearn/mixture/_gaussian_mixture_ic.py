"""GaussianMixtureIC"""

# Author: Thomas Athey <tathey1@jhmi.edu>
# Modified by: Benjamin Pedigo <bpedigo@jhu.edu>
#              Tingshan Liu <tliu68@jhmi.edu>


import numpy as np
import joblib
from ..utils.fixes import parse_version
from ..utils.parallel import delayed, Parallel
from ..utils import check_scalar
from ..utils.validation import check_is_fitted, check_random_state

from . import GaussianMixture
from ..base import BaseEstimator, ClusterMixin
from ..model_selection import ParameterGrid




class GaussianMixtureIC(ClusterMixin, BaseEstimator):
    """Gaussian mixture with BIC/AIC.

    Automatic Gaussian Mixture Model (GMM) selection via the
    Bayesian Information Criterion (BIC)
    or the Akaike Information Criterion (AIC).

    Different combinations of initialization, GMM,
    and cluster numbers are used and the clustering
    with the best selection criterion (BIC or AIC) is chosen.

    Parameters
    ----------
    min_components : int, default=2
        The minimum number of mixture components to consider.
        If ``max_components`` is not None, ``min_components`` must be
        less than or equal to ``max_components``.

    max_components : int or None, default=10
        The maximum number of mixture components to consider.
        Must be greater than or equal to ``min_components``.

    covariance_type : {'full', 'tied', 'diag', 'spherical', 'all' (default)},
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

    random_state : int, RandomState instance or None, optional (default=None)
        There is randomness in k-means initialization of
        :class:`sklearn.mixture.GaussianMixture`. This parameter is passed to
        :class:`~sklearn.mixture.GaussianMixture` to control the random state.
        If int, ``random_state`` is used as the random number generator seed.
        If RandomState instance, ``random_state`` is the random number
        generator; If None, the random number generator is the
        RandomState instance used by ``np.random``.

    n_init : int, optional (default = 1)
        If ``n_init`` is larger than 1, additional
        ``n_init``-1 runs of :class:`sklearn.mixture.GaussianMixture`
        initialized with k-means will be performed
        for all covariance parameters in ``covariance_type``.

    init_params : {‘kmeans’ (default), ‘k-means++’, ‘random’, ‘random_from_data’}
        The method used to initialize the weights, the means and the precisions
        for Gaussian mixture modeling.

    max_iter : int, optional (default = 100)
        The maximum number of EM iterations to perform.

    verbose : int, optional (default = 0)
        Enable verbose output. If 1 then it prints the current initialization
        and each iteration step. If greater than 1 then it prints also
        the log probability and the time needed for each step.

    criterion : str {"bic" or "aic"}, optional, (default = "bic")
        Select the best model based on Bayesian Information Criterion (bic) or
        Aikake Information Criterion (aic).


    Attributes
    ----------
    best_criterion_ : float
        The best (lowest) Bayesian or Aikake Information Criterion.

    n_components_ : int
        Number of clusters for the model with the best bic/aic.

    covariance_type_ : str
        Covariance type for the model with the best bic/aic.

    best_model_ : :class:`sklearn.mixture.GaussianMixture`
        Object with the best bic/aic.

    labels_ : array-like, shape (n_samples,)
        Labels of each point predicted by the best model.

    n_iter_ : int
        Number of step used by the best fit of EM for the best model
        to reach the convergence.

    results_ : list
        Contains exhaustive information about all the clustering runs.
        Each item represents a class object storing the results
        for a single run with the following attributes:

        model : :class:`~sklearn.mixture.GaussianMixture` object
            GMM clustering fit to the data.
        criterion_score : float
            Bayesian or Aikake Information Criterion score.
        n_components : int
            Number of clusters.
        covariance_type : {'full', 'tied', 'diag', 'spherical'}
            Covariance type used for the GMM.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

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
        Gaussian Finite Mixture Models. The R journal, 8(1), 289–317.
        <https://doi.org/10.32614/RJ-2016-021>_`

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.mixture import GaussianMixtureIC
    >>> X = np.array([[1, 2], [1, 4], [1, 0], [10, 2], [10, 4], [10, 0]])
    >>> gmIC = GaussianMixtureIC(max_components=4, random_state=0)
    >>> gmIC.fit_predict(X)
    array([0, 0, 0, 1, 1, 1])
    >>> print(gmIC.n_components_)
    2
    """

    def __init__(
        self,
        min_components=2,
        max_components=10,
        covariance_type="all",
        random_state=None,
        n_init=1,
        init_params="kmeans",
        max_iter=100,
        verbose=0,
        criterion="bic",
        n_jobs=None,
    ):

        self.min_components = min_components
        self.max_components = max_components
        self.covariance_type = covariance_type
        self.random_state = random_state
        self.n_init = n_init
        self.init_params = init_params
        self.max_iter = max_iter
        self.verbose = verbose
        self.criterion = criterion
        self.n_jobs = n_jobs

    def _check_multi_comp_inputs(self, input, name, default):
        if isinstance(input, (np.ndarray, list)):
            input = list(np.unique(input))
        elif isinstance(input, str):
            if input not in default:
                raise ValueError(f"{name} is {input} but must be one of {default}.")
            if input != "all":
                input = [input]
            else:
                input = default.copy()
                input.remove("all")
        else:
            raise TypeError(
                f"{name} is a {type(input)} but must be a numpy array, "
                "a list, or a string."
            )
        return input

    def _check_parameters(self):
        check_scalar(
            self.min_components,
            min_val=1,
            max_val=self.max_components,
            name="min_components",
            target_type=int,
        )
        check_scalar(
            self.max_components, min_val=1, name="max_components", target_type=int
        )

        covariance_type = self._check_multi_comp_inputs(
            self.covariance_type,
            "covariance_type",
            ["spherical", "diag", "tied", "full", "all"],
        )

        check_scalar(self.n_init, name="n_init", target_type=int, min_val=1)

        if self.criterion not in ["aic", "bic"]:
            raise ValueError(
                f'criterion is {self.criterion} but must be "aic" or "bic".'
            )

        return covariance_type


    def _fit_cluster(self, X, gm_params, seed):
        gm_params["init_params"] = self.init_params
        gm_params["max_iter"] = self.max_iter
        gm_params["n_init"] = self.n_init
        gm_params["random_state"] = seed

        model = GaussianMixture(**gm_params)
        model.fit(X)

        if self.criterion == "bic":
            criterion_value = model.bic(X)
        else:
            criterion_value = model.aic(X)

        # change the precision of "criterion_value" based on sample size
        criterion_value = round(criterion_value, int(np.log10(X.shape[0])))
        results = _CollectResults(model, criterion_value, gm_params)
        return results

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

        covariance_type = self._check_parameters()
        X = self._validate_data(X, dtype=[np.float64, np.float32], ensure_min_samples=1)

        random_state = check_random_state(self.random_state)

        # check n_components against sample size
        n_comps = [self.max_components, self.min_components]
        names = ["max_components", "min_components"]
        for i in range(len(names)):
            if n_comps[i] > X.shape[0]:
                msg = names[i] + "must be <= n_samples, but" + names[i]
                msg += "= {}, n_samples = {}".format(n_comps[i], X.shape[0])
                raise ValueError(msg)

        param_grid = dict(
            covariance_type=covariance_type,
            n_components=range(self.min_components, self.max_components + 1),
        )
        param_grid = list(ParameterGrid(param_grid))

        seeds = random_state.randint(np.iinfo(np.int32).max, size=len(param_grid))

        if parse_version(joblib.__version__) < parse_version("0.12"):
            parallel_kwargs = {"backend": "threading"}
        else:
            parallel_kwargs = {"prefer": "threads"}

        results = Parallel(n_jobs=self.n_jobs, verbose=self.verbose, **parallel_kwargs)(
            delayed(self._fit_cluster)(X, gm_params, seed)
            for gm_params, seed in zip(param_grid, seeds)
        )
        best_criter = [result.criterion for result in results]

        if sum(best_criter == np.min(best_criter)) == 1:
            best_idx = np.argmin(best_criter)
        else:
            # in case there is a tie,
            # select the model with the least number of parameters
            ties = np.where(best_criter == np.min(best_criter))[0]
            n_params = [results[tie].model._n_parameters() for tie in ties]
            best_idx = ties[np.argmin(n_params)]

        self.best_criterion_ = results[best_idx].criterion
        self.n_components_ = results[best_idx].n_components
        self.covariance_type_ = results[best_idx].covariance_type
        self.best_model_ = results[best_idx].model
        self.n_iter_ = results[best_idx].model.n_iter_
        self.labels_ = results[best_idx].model.predict(X)
        self.results_ = results
        self.n_features_in_ = X.shape[1]

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
        check_is_fitted(self, ["best_model_"], all_or_any=all)
        X = self._validate_data(X, reset=False)
        labels = self.best_model_.predict(X)

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
    


class _CollectResults:
    """Collect intermediary results.

    Represent the intermediary results for a single GMM clustering run
    of :class:`sklearn.mixture.GaussianMixtureIC`

    Attributes
    ----------

    model : GaussianMixture object
        GMM clustering fit to the data.

    criterion : float
        Bayesian or Aikake Information Criterion.

    n_components : int
        Number of components.

    covariance_type : {'full', 'tied', 'diag', 'spherical'}
        Covariance type used for the GMM.

    """

    def __init__(self, model, criter, gm_params):
        self.model = model
        self.criterion = criter
        self.n_components = gm_params["n_components"]
        self.covariance_type = gm_params["covariance_type"]