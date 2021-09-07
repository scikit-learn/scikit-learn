"""GaussianMixtureIC"""

# Author: Thomas Athey <tathey1@jhmi.edu>
# Modified by: Benjamin Pedigo <bpedigo@jhu.edu>
#              Tingshan Liu <tliu68@jhmi.edu>
#              Giorgio Di Salvo <gdisalv1@jhu.edu>


import numpy as np
import joblib
from joblib import Parallel
from ..utils.fixes import delayed, parse_version
import warnings

from . import GaussianMixture
from ._gaussian_mixture import (
    _compute_precision_cholesky,
    _estimate_gaussian_parameters,
)
from ..cluster import AgglomerativeClustering
from ..model_selection import ParameterGrid
from ..utils import check_scalar
from ..utils.validation import check_is_fitted, check_random_state
from ..preprocessing import OneHotEncoder
from ..exceptions import ConvergenceWarning
from ..base import BaseEstimator, ClusterMixin


class GaussianMixtureIC(ClusterMixin, BaseEstimator):
    """Gaussian mixture with BIC/AIC.

    Automatic Gaussian Mixture Model (GMM) selection via the
    Bayesian Information Criterion (BIC)
    or the Akaike Information Criterion (AIC).

    Clustering algorithm using a hierarchical agglomerative clustering
    or k-means clustering prior to a Gaussian mixture model (GMM) fitting.
    Different combinations of initialization, agglomeration, GMM,
    and cluster numbers are used and the clustering
    with the best selection criterion (BIC or AIC) is chosen.

    Parameters
    ----------
    min_components : int, default=2
        The minimum number of mixture components to consider.
        If ``max_components`` is not None, ``min_components`` must be
        less than or equal to ``max_components``. If ``label_init`` is given,
        ``min_components`` must match the number of unique labels
        in ``label_init``.

    max_components : int or None, default=10
        The maximum number of mixture components to consider.
        Must be greater than or equal to ``min_components``.
        If ``label_init`` is given, ``max_components`` must match
        the number of unique labels in ``label_init``.

    affinity : {'euclidean', 'manhattan', 'cosine', 'all', 'none' (default)},
            optional
        String or list/array describing the type of affinities to use
        in agglomeration. If a string, it must be one of:

        - 'euclidean'
            L2 norm
        - 'manhattan'
            L1 norm
        - 'cosine'
            cosine similarity
        - 'all'
            considers all affinities in
            ['euclidean', 'manhattan', 'cosine', 'none']
        - 'none'
            no agglomeration - GMM is initialized with k-means

        If a list/array, it must be a list/array of strings containing only
        'euclidean', 'manhattan', 'cosine', and/or 'none'.

        Note that cosine similarity can only work when all of the rows are not
        the zero vector. If the input matrix has a zero row, cosine similarity
        will be skipped and a warning will be thrown.

    linkage : {'ward', 'complete', 'average', 'single', 'all' (default)},
            optional
        String or list/array describing the type of linkages to use
        in agglomeration. Not used if ``affinity`` is 'none'.
        If a string, it must be one of:

        - 'ward'
            ward's clustering, can only be used with euclidean affinity
        - 'complete'
            complete linkage
        - 'average'
            average linkage
        - 'single'
            single linkage
        - 'all'
            considers all linkages in ['ward', 'complete', 'average', 'single']

        If a list/array, it must be a list/array of strings containing only
        'ward', 'complete', 'average', and/or 'single'.

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

    label_init : array-like, shape (n_samples,), optional (default=None)
        List of labels for samples if available. Used to initialize the model.
        If provided, ``min_components`` and ``max_components`` must match
        the number of unique labels given here.

    n_init : int, optional (default = 1)
        If ``n_init`` is larger than 1 and ``label_init`` is None, additional
        ``n_init``-1 runs of :class:`sklearn.mixture.GaussianMixture`
        initialized with k-means will be performed
        for all covariance parameters in ``covariance_type``.

    max_iter : int, optional (default = 100)
        The maximum number of EM iterations to perform.

    verbose : int, optional (default = 0)
        Enable verbose output. If 1 then it prints the current initialization
        and each iteration step. If greater than 1 then it prints also
        the log probability and the time needed for each step.

    criterion : str {"bic" or "aic"}, optional, (default = "bic")
        Select the best model based on Bayesian Information Criterion (bic) or
        Aikake Information Criterion (aic).

    max_agglom_size : int or None, optional (default = 2000)
        The maximum number of datapoints on which to do agglomerative
        clustering as the initialization to GMM. If the number of datapoints
        is larger than this value, a random subset of the data is used
        for agglomerative initialization. If None, all data is used
        for agglomerative clustering for initialization.

    n_jobs : int or None, optional (default = None)
        The number of jobs to use for the computation. This works by computing
        each of the initialization runs in parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    best_criterion_ : float
        The best (lowest) Bayesian or Aikake Information Criterion.

    n_components_ : int
        Number of clusters for the model with the best bic/aic.

    covariance_type_ : str
        Covariance type for the model with the best bic/aic.

    affinity_ : str
        Affinity used for the model with the best bic/aic.

    linkage_ : str
        Linkage used for the model with the best bic/aic.

    reg_covar_ : float
        Regularization used in the model with the best bic/aic.

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
        affinity : {'euclidean', 'manhattan', 'cosine', 'none'}
            Affinity used for the Agglomerative Clustering.
        linkage : {'ward', 'complete', 'average', 'single', 'none'}
            Linkage used for the Agglomerative Clustering.
        covariance_type : {'full', 'tied', 'diag', 'spherical'}
            Covariance type used for the GMM.
        reg_covar : float
            Regularization used for the GMM.

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
        affinity="none",
        linkage="all",
        covariance_type="all",
        random_state=None,
        label_init=None,
        n_init=1,
        max_iter=100,
        verbose=0,
        criterion="bic",
        max_agglom_size=2000,
        n_jobs=None,
    ):

        self.min_components = min_components
        self.max_components = max_components
        self.affinity = affinity
        self.linkage = linkage
        self.covariance_type = covariance_type
        self.random_state = random_state
        self.label_init = label_init
        self.n_init = n_init
        self.max_iter = max_iter
        self.verbose = verbose
        self.criterion = criterion
        self.max_agglom_size = max_agglom_size
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

        affinity = self._check_multi_comp_inputs(
            self.affinity,
            "affinity",
            ["euclidean", "manhattan", "cosine", "none", "all"],
        )

        if ("ward" in self.linkage) and ("euclidean" not in self.affinity):
            raise ValueError(
                'If "ward" is a linkage option, "euclidean" must be an affinity option.'
            )

        linkage = self._check_multi_comp_inputs(
            self.linkage, "linkage", ["ward", "complete", "average", "single", "all"]
        )

        covariance_type = self._check_multi_comp_inputs(
            self.covariance_type,
            "covariance_type",
            ["spherical", "diag", "tied", "full", "all"],
        )

        if isinstance(self.label_init, list):
            self.label_init = np.array(self.label_init)
        elif isinstance(self.label_init, np.ndarray):
            if self.label_init.ndim > 2 or (
                self.label_init.ndim == 2 and 1 not in self.label_init.shape
            ):
                raise TypeError("label_init must be a one dimension array.")
        elif self.label_init is not None:
            raise TypeError(
                f"label_init is a {type(self.label_init)} but must be a "
                "1-D number array, a list, or None."
            )

        # Adjust elements in label_init to range(n_components of label_init)
        if self.label_init is not None:
            check_scalar(
                np.size(np.unique(self.label_init)),
                target_type=int,
                name="the number of init labels",
                min_val=self.min_components,
                max_val=self.max_components,
            )
            _, labels_init = np.unique(self.label_init, return_inverse=True)
        else:
            labels_init = None
        self.label_init = labels_init

        check_scalar(self.n_init, name="n_init", target_type=int, min_val=1)

        if self.criterion not in ["aic", "bic"]:
            raise ValueError(
                f'criterion is {self.criterion} but must be "aic" or "bic".'
            )

        check_scalar(
            self.max_agglom_size, name="max_agglom_size", target_type=int, min_val=2
        )

        return affinity, linkage, covariance_type

    def _init_gm_params(self, X, labels, gm_params):
        labels = labels.reshape(-1, 1)
        enc = OneHotEncoder()
        enc.fit(labels)
        onehot = enc.fit_transform(labels).toarray()
        weights_init, means_init, precisions_init = _onehot_to_initial_params(
            X, onehot, gm_params["covariance_type"]
        )
        gm_params["weights_init"] = weights_init
        gm_params["means_init"] = means_init
        gm_params["precisions_init"] = precisions_init
        return gm_params

    def _fit_cluster(self, X, X_subset, y, ag_params, gm_params, agg_clustering, seed):
        label_init = self.label_init
        n_samples = X.shape[0]
        if label_init is not None:
            gm_params = self._init_gm_params(X, label_init, gm_params)
        elif ag_params["affinity"] != "none":
            gm_params = self._init_gm_params(X_subset, agg_clustering, gm_params)
            n_samples = X_subset.shape[0]
        else:
            gm_params["init_params"] = "kmeans"
        gm_params["max_iter"] = self.max_iter
        gm_params["random_state"] = seed

        # if none of the iterations converge, bic/aic is set to inf
        criter = np.inf
        # below is the regularization scheme
        reg_covar = 0
        while reg_covar <= 1 and criter == np.inf:
            model = GaussianMixture(**gm_params)
            model.reg_covar = reg_covar
            try:
                # ignoring warning here because if convergence is not reached,
                # the regularization is automatically increased
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ConvergenceWarning)
                    model.fit(X)
                predictions = model.predict(X)
                counts = [
                    sum(predictions == i) for i in range(gm_params["n_components"])
                ]
                # singleton clusters not allowed
                assert not any([count <= 1 for count in counts])

            except ValueError:
                reg_covar = _increase_reg(reg_covar)
                continue
            except AssertionError:
                reg_covar = _increase_reg(reg_covar)
                continue
            # if the code gets here, then the model has been fit with
            # no errors or singleton clusters
            if self.criterion == "bic":
                criter = model.bic(X)
            else:
                criter = model.aic(X)
            break

        gm_params["reg_covar"] = reg_covar
        # change the precision of "criter" based on sample size
        criter = round(criter, int(np.log10(n_samples)))
        results = _CollectResults(model, criter, gm_params, ag_params)
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

        affinity, linkage, covariance_type = self._check_parameters()
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

        # check if X contains the 0 vector
        if np.any(~X.any(axis=1)):
            if affinity == ["cosine"]:
                raise ValueError(
                    "X contains a zero vector, cannot run cosine affinity."
                )
            elif "cosine" in affinity:
                affinity.remove("cosine")
                warnings.warn("X contains a zero vector, will not run cosine affinity.")

        label_init = self.label_init
        if label_init is not None:
            if label_init.size != X.shape[0]:
                raise ValueError(
                    "n_samples must be the same as the length of label_init"
                )

        param_grid = dict(
            affinity=affinity,
            linkage=linkage,
            covariance_type=covariance_type,
            n_components=range(self.min_components, self.max_components + 1),
        )

        param_grid = list(ParameterGrid(param_grid))
        param_grid_ag, param_grid = _process_paramgrid(
            param_grid, self.n_init, self.label_init
        )

        seeds = random_state.randint(np.iinfo(np.int32).max, size=len(param_grid))

        n_samples = X.shape[0]
        if self.max_agglom_size is None or n_samples <= self.max_agglom_size:
            X_subset = X
        else:  # if dataset is huge, agglomerate a subset
            subset_idxs = random_state.choice(
                np.arange(0, n_samples), self.max_agglom_size
            )
            X_subset = X[subset_idxs, :]

        ag_labels = []
        if self.label_init is None:
            for p_ag in param_grid_ag:
                if p_ag["affinity"] != "none":
                    agg = AgglomerativeClustering(
                        n_clusters=self.min_components, **p_ag
                    )
                    agg.fit(X_subset)
                    hierarchical_labels = _hierarchical_labels(
                        agg.children_, self.min_components, self.max_components
                    )
                    ag_labels.append(hierarchical_labels)

        def _fit_for_data(ag_params, gm_params, seed):
            n_clusters = gm_params["n_components"]
            if (ag_params["affinity"] != "none") and (self.label_init is None):
                index = param_grid_ag.index(ag_params)
                agg_clustering = ag_labels[index][:, n_clusters - self.min_components]
            else:
                agg_clustering = []
            return self._fit_cluster(
                X, X_subset, y, ag_params, gm_params, agg_clustering, seed
            )

        if parse_version(joblib.__version__) < parse_version("0.12"):
            parallel_kwargs = {"backend": "threading"}
        else:
            parallel_kwargs = {"prefer": "threads"}

        results = Parallel(n_jobs=self.n_jobs, verbose=self.verbose, **parallel_kwargs)(
            delayed(_fit_for_data)(ag_params, gm_params, seed)
            for (ag_params, gm_params), seed in zip(param_grid, seeds)
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
        self.affinity_ = results[best_idx].affinity
        self.linkage_ = results[best_idx].linkage
        self.reg_covar_ = results[best_idx].reg_covar
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


def _increase_reg(reg):
    """Scale the regularization factor by 10.

    Parameters
    ----------
    reg: float
        Current regularization factor.

    Returns
    -------
    reg : float
        Increased regularization
    """
    if reg == 0:
        reg = 1e-06
    else:
        reg *= 10
    return reg


def _onehot_to_initial_params(X, onehot, cov_type):
    """Compute initial parameters.

    Compute cluster weights, cluster means and cluster precisions from
    a given clustering.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row
        corresponds to a single data point.
    onehot : ndarray, shape (n_samples, n_clusters)
        Each row has a 1 indicating cluster membership, other entries are 0.
    cov_type : {'full', 'tied', 'diag', 'spherical'}
        Covariance type for :class:`sklearn.mixture.GaussianMixture`.
    """
    n_samples = X.shape[0]
    weights, means, covariances = _estimate_gaussian_parameters(
        X, onehot, 1e-06, cov_type
    )
    weights /= n_samples

    precisions_cholesky_ = _compute_precision_cholesky(covariances, cov_type)

    if cov_type == "tied":
        c = precisions_cholesky_
        precisions = np.dot(c, c.T)
    elif cov_type == "diag":
        precisions = precisions_cholesky_
    else:  # "spherical" or "full"
        precisions = [np.dot(c, c.T) for c in precisions_cholesky_]

    return weights, means, precisions


def _process_paramgrid(paramgrid, n_init, label_init):
    """Remove combinations of affinity and linkage that are not possible.

    Parameters
    ----------
    paramgrid : list of dicts
        Each dict has the keys 'affinity', 'covariance_type', 'linkage',
        'n_components', and 'random_state'

    Returns
    -------
    paramgrid_processed : list pairs of dicts
        For each pair, the first dict are the options
        for :class:`sklearn.cluster.AgglomerativeClustering`.
        The second dict include the options
        for :class:`sklearn.mixture.GaussianMixture`.
    ag_paramgrid_processed : list of dicts
        Options for :class:`~sklearn.cluster.AgglomerativeClustering`.
    """
    gm_keys = ["covariance_type", "n_components"]
    ag_keys = ["affinity", "linkage"]
    ag_params_processed = []
    paramgrid_processed = []

    for params in paramgrid:
        if (
            params["affinity"] == "none"
            and params["linkage"] != paramgrid[0]["linkage"]
        ):
            continue
        elif (
            params["linkage"] == "ward"
            and params["affinity"] != "euclidean"
            and params["affinity"] != "none"
        ):
            continue
        else:
            gm_params = {key: params[key] for key in gm_keys}
            ag_params = {key: params[key] for key in ag_keys}
            paramgrid_processed.append([ag_params, gm_params])

            if ag_params not in ag_params_processed:
                ag_params_processed.append(ag_params)

            if ag_params["affinity"] == "none" and n_init > 1 and label_init is None:
                more_kmeans_init = gm_params.copy()
                more_kmeans_init.update({"n_init": 1})
                paramgrid_processed += [
                    [{"affinity": "none", "linkage": "none"}, more_kmeans_init]
                ] * (n_init - 1)

    return ag_params_processed, paramgrid_processed


def _hierarchical_labels(children, min_components, max_components):
    n_samples = len(children) + 1
    hierarchical_labels = np.arange(n_samples).reshape((-1, 1))
    merge_start = n_samples - max_components - 1
    merge_end = n_samples - min_components - 1

    for n in range(merge_end + 1):
        inds = np.where(np.isin(hierarchical_labels[:, n], children[n, :]))[0]
        hierarchical_labels[inds, -1] = n_samples + n
        if n < merge_end:
            hierarchical_labels = np.hstack(
                (hierarchical_labels, hierarchical_labels[:, -1].reshape((-1, 1)))
            )

    if n_samples == max_components:
        hierarchical_labels = np.hstack(
            (np.arange(max_components).reshape(-1, 1), hierarchical_labels)
        )
    else:
        hierarchical_labels = hierarchical_labels[:, merge_start:]
    for i in range(hierarchical_labels.shape[1]):
        _, hierarchical_labels[:, i] = np.unique(
            hierarchical_labels[:, i], return_inverse=True
        )

    return hierarchical_labels[:, ::-1]


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

    affinity : {'euclidean', 'manhattan', 'cosine', 'none'}
        Affinity used for the Agglomerative Clustering.

    linkage : {'ward', 'complete', 'average', 'single'}
        Linkage used for the Agglomerative Clustering.

    covariance_type : {'full', 'tied', 'diag', 'spherical'}
        Covariance type used for the GMM.

    reg_covar : float
        Regularization factor used for the regularization of the
        GMM covariance matrices.
    """

    def __init__(self, model, criter, gm_params, ag_params):
        self.model = model
        self.criterion = criter
        self.n_components = gm_params["n_components"]
        self.affinity = ag_params["affinity"]
        self.linkage = ag_params["linkage"]
        self.covariance_type = gm_params["covariance_type"]
        self.reg_covar = gm_params["reg_covar"]
