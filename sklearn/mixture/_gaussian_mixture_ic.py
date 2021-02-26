"""GaussianMixtureIC"""

# Author: Thomas Athey <tathey1@jhmi.edu>
# Modified by Benjamin Pedigo <bpedigo@jhu,edu>, Tingshan Liu <tliu68@jhmi.edu>


import numpy as np
from joblib import Parallel
from ..utils.fixes import delayed
import warnings

from . import GaussianMixture
from ._gaussian_mixture import (
    _compute_precision_cholesky,
    _estimate_gaussian_parameters,
)
from ..cluster import AgglomerativeClustering
from ..model_selection import ParameterGrid
from ..utils import check_scalar, check_array
from ..utils.validation import check_is_fitted
from ..exceptions import ConvergenceWarning
from ..base import BaseEstimator, ClusterMixin


class GaussianMixtureIC(ClusterMixin, BaseEstimator):
    """
    Automatic Gaussian Mixture Model (GMM) selection via BIC/AIC.

    Clustering algorithm using a hierarchical agglomerative clustering
    then Gaussian mixtured model (GMM) fitting. Different combinations
    of agglomeration, GMM, and cluster numbers are used and the clustering
    with the best selection criterion (bic/aic) is chosen.

    Parameters
    ----------
    min_components : int, default=2.
        The minimum number of mixture components to consider.
        If ``max_components`` is not None, ``min_components`` must be
        less than or equal to ``max_components``. If ``label_init`` is given,
        min_components must match number of unique labels in ``label_init``.

    max_components : int or None, default=10.
        The maximum number of mixture components to consider. Must be greater
        than or equal to min_components.
        If ``label_init`` is given, ``min_components`` must match number of
        unique labels in ``label_init``.

    affinity : {'euclidean','manhattan','cosine','none', 'all' (default)},
        optional
        String or list/array describing the type of affinities to use
        in agglomeration. If a string, it must be one of:

        - 'euclidean'
            L2 norm
        - 'manhattan'
            L1 norm
        - 'cosine'
            cosine similarity
        - 'none'
            no agglomeration - GMM is initialized with k-means
        - 'all'
            considers all affinities in
            ['euclidean','manhattan','cosine','none']

        If a list/array, it must be a list/array of strings containing only
        'euclidean', 'manhattan', 'cosine', and/or 'none'.

        Note that cosine similarity can only work when all of the rows are not
        the zero vector. If the input matrix has a zero row, cosine similarity
        will be skipped and a warning will be thrown.

    linkage : {'ward','complete','average','single', 'all' (default)}, optional
        String or list/array describing the type of linkages to use
        in agglomeration. If a string, it must be one of:

        - 'ward'
            ward's clustering, can only be used with euclidean affinity
        - 'complete'
            complete linkage
        - 'average'
            average linkage
        - 'single'
            single linkage
        - 'all'
            considers all linkages in ['ward','complete','average','single']

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
        If int, ``random_state`` is used as the random number generator seed;
        If RandomState instance, ``random_state`` is the random number
        generator; If None, the random number generator is the
        RandomState instance used by ``np.random``.

    label_init : array-like, shape (n_samples,), optional (default=None)
        List of labels for samples if available. Used to initialize the model.
        If provided, ``min_components`` and ``max_components`` must match
        the number of unique labels given here.

    max_iter : int, optional (default = 100).
        The maximum number of EM iterations to perform.

    criterion : str {"bic" or "aic"}, optional, (default="bic")
        select the best model based on Bayesian Information Criterion (bic) or
        Aikake Information Criterion (aic)

    verbose : int, optional (default = 0)
        Enable verbose output. If 1 then it prints the current initialization
        and each iteration step. If greater than 1 then it prints also
        the log probability and the time needed for each step.

    max_agglom_size : int or None, optional (default = 2000)
        The maximum number of datapoints on which to do agglomerative
        clustering as the initialization to GMM. If the number of datapoints
        is larger than this value, a random subset of the data is used
        for agglomerative initialization. If None, all data is used
        for agglomerative clustering for initialization.

    n_jobs : int or None, optional (default = None)
        The number of jobs to use for the computation. This works by computing
        each of the initialization runs in parallel. None means 1 unless in a
        ``joblib.parallel_backend context``. -1 means using all processors.
        See https://scikit-learn.org/stable/glossary.html#term-n-jobs
        for more details.

    Attributes
    ----------
    results_ : dict
        Contains exhaustive information about all the clustering runs.
        Keys are:

        'model' : GaussianMixture object
            GMM clustering fit to the data
        'bic/aic' : float
            Bayesian Information Criterion
        'n_components' : int
            number of clusters
        'affinity' : {'euclidean','manhattan','cosine','none'}
            affinity used in Agglomerative Clustering
        'linkage' : {'ward','complete','average','single'}
            linkage used in Agglomerative Clustering
        'covariance_type' : {'full', 'tied', 'diag', 'spherical'}
            covariance type used in GMM
        'reg_covar' : float
            regularization used in GMM

    best_criterion_ : float
        the best (lowest) Bayesian or Aikake Information Criterion

    n_components_ : int
        number of clusters in the model with the best bic/aic

    covariance_type_ : str
        covariance type in the model with the best bic/aic

    affinity_ : str
        affinity used in the model with the best bic/aic

    linkage_ : str
        linkage used in the model with the best bic/aic

    reg_covar_ : float
        regularization used in the model with the best bic/aic

    best_model_ : :class:`sklearn.mixture.GaussianMixture`
        object with the best bic/aic

    labels_ : array-like, shape (n_samples,)
        labels of each point predicted by the best model

    n_iter_ : int
        number of step used by the best fit of EM for the best model
        to reach the convergence

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

    See Also
    --------
    GaussianMixture : fit Gaussian mixture model
    BayesianGaussianMixture : Gaussian mixture model fit with a variational
    inference

    Notes
    -----
    This algorithm was strongly inspired by mclust, a clustering package in R

    References
    ----------
    .. [1] Jeffrey D. Banfield and Adrian E. Raftery. Model-based gaussian and
       non-gaussian clustering. Biometrics, 49:803–821, 1993.

    .. [2] Abhijit Dasgupta and Adrian E. Raftery. Detecting features in
       spatial point processes with clutter via model-based clustering.
       Journal of the American Statistical Association, 93(441):294–302, 1998.
    """

    def __init__(
        self,
        min_components=2,
        max_components=10,
        affinity="all",
        linkage="all",
        covariance_type="all",
        random_state=None,
        label_init=None,
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
        self.max_iter = max_iter
        self.label_init = label_init
        self.verbose = verbose
        self.criterion = criterion
        self.max_agglom_size = max_agglom_size
        self.n_jobs = n_jobs

    def _check_multi_comp_inputs(self, input, name, default):
        if isinstance(input, (np.ndarray, list)):
            input = np.unique(input)
        elif isinstance(input, str):
            if input not in default:
                msg = (
                    name
                    + " must be one of "
                    + str(default)
                    + " not {}".format(input)
                )
                raise ValueError(msg)
            if input != "all":
                input = [input]
        else:
            msg = name + "must be a numpy array, a list, or "
            msg += "string, not {}".format(type(input))
            raise TypeError(msg)
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
            self.max_components,
            min_val=1,
            name="max_components",
            target_type=int,
        )

        self.affinity = self._check_multi_comp_inputs(
            self.affinity,
            "affinity",
            ["euclidean", "manhattan", "cosine", "none", "all"],
        )

        if ("ward" in self.linkage) and not ("euclidean" in self.affinity):
            msg = (
                'if "ward" is a linkage option, '
                + '"euclidean" must be an affinity option'
            )
            raise ValueError(msg)

        self.linkage = self._check_multi_comp_inputs(
            self.linkage,
            "linkage",
            ["ward", "complete", "average", "single", "all"],
        )

        self.covariance_type = self._check_multi_comp_inputs(
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
                msg = "label_init must be a one dimension array."
                raise TypeError(msg)
        elif self.label_init is not None:
            msg = "label_init must be a 1-D numpy array, a list, or None,"
            msg += "not {}".format(type(self.label_init))
            raise TypeError(msg)

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

        if self.criterion not in ["aic", "bic"]:
            msg = "criterion must be one of " + '["aic, "bic"]'
            msg += " not {}".format(self.criterion)
            raise ValueError(msg)

        check_scalar(
            self.max_agglom_size,
            name="max_agglom_size",
            target_type=int,
            min_val=2,
        )

    def _init_gm_params(self, X, labels, params):
        onehot = _labels_to_onehot(labels)
        weights_init, means_init, precisions_init = _onehot_to_initial_params(
            X, onehot, params[1]["covariance_type"]
        )
        gm_params = params[1]
        gm_params["weights_init"] = weights_init
        gm_params["means_init"] = means_init
        gm_params["precisions_init"] = precisions_init
        return gm_params

    def _fit_cluster(self, X, X_subset, y, params, agg_clustering):
        label_init = self.label_init
        if label_init is not None:
            gm_params = self._init_gm_params(X, label_init, params)
        elif params[0]["affinity"] != "none":
            gm_params = self._init_gm_params(X_subset, agg_clustering, params)
        else:
            gm_params = params[1]
            gm_params["init_params"] = "kmeans"
        gm_params["reg_covar"] = 0
        gm_params["max_iter"] = self.max_iter

        # if none of the iterations converge, bic/aic is set to inf
        criter = np.inf
        # below is the regularization scheme
        while gm_params["reg_covar"] <= 1 and criter == np.inf:
            model = GaussianMixture(**gm_params)
            try:
                # ignoring warning here because if convergence is not reached,
                # the regularization is automatically increased
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", ConvergenceWarning)
                    model.fit(X)
                predictions = model.predict(X)
                counts = [
                    sum(predictions == i)
                    for i in range(gm_params["n_components"])
                ]
                # singleton clusters not allowed
                assert not any([count <= 1 for count in counts])

            except ValueError:
                gm_params["reg_covar"] = _increase_reg(gm_params["reg_covar"])
                continue
            except AssertionError:
                gm_params["reg_covar"] = _increase_reg(gm_params["reg_covar"])
                continue
            # if the code gets here, then the model has been fit with
            # no errors or singleton clusters
            if self.criterion == "bic":
                criter = model.bic(X)
            else:
                criter = model.aic(X)
            break

        results = {
            "model": model,
            "bic/aic": criter,
            "n_components": gm_params["n_components"],
            "affinity": params[0]["affinity"],
            "linkage": params[0]["linkage"],
            "covariance_type": gm_params["covariance_type"],
            "reg_covar": gm_params["reg_covar"],
        }
        return results

    def fit(self, X, y=None):
        """
        Fits Gaussian mixture model to the data.
        Initialize with agglomerative clustering then
        estimate model parameters with EM algorithm.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            List of n_features-dimensional data points. Each row
            corresponds to a single data point.

        Returns
        -------
        self : object
            Returns an instance of self.
        """

        self._check_parameters()
        X = check_array(
            X, dtype=[np.float64, np.float32], ensure_min_samples=1
        )

        # check n_components against sample size
        n_comps = [self.max_components, self.min_components]
        names = ["max_components", "min_components"]
        for i in range(len(names)):
            if n_comps[i] > X.shape[0]:
                msg = names[i] + "must be <= n_samples, but" + names[i]
                msg += "= {}, n_samples = {}".format(n_comps[i], X.shape[0])
                raise ValueError(msg)

        # check if X contains the 0 vector
        if np.any(~X.any(axis=1)) and (
            "cosine" in self.affinity or self.affinity == "all"
        ):
            if isinstance(self.affinity, np.ndarray):
                self.affinity = np.delete(
                    self.affinity, np.argwhere(self.affinity == "cosine")
                )
            if isinstance(self.affinity, list):
                self.affinity.remove("cosine")
            if self.affinity == "all":
                self.affinity = ["euclidean", "manhattan", "none"]
            warnings.warn(
                "X contains a zero vector, will not run cosine affinity."
            )

        label_init = self.label_init
        if label_init is not None:
            if label_init.size != X.shape[0]:
                msg = "n_samples must be the same as the length of label_init"
                raise ValueError(msg)

        param_grid = dict(
            n_components=range(self.min_components, self.max_components + 1),
            random_state=[self.random_state],
        )
        param_grid["affinity"] = (
            self.affinity
            if self.affinity != "all"
            else ["euclidean", "manhattan", "cosine", "none"]
        )
        param_grid["linkage"] = (
            self.linkage
            if self.linkage != "all"
            else ["ward", "complete", "average", "single"]
        )
        param_grid["covariance_type"] = (
            self.covariance_type
            if self.covariance_type != "all"
            else ["spherical", "diag", "tied", "full"]
        )

        param_grid = list(ParameterGrid(param_grid))
        param_grid_ag, param_grid = _process_paramgrid(param_grid)

        n = X.shape[0]
        if self.max_agglom_size is None or n <= self.max_agglom_size:
            X_subset = X
        else:  # if dataset is huge, agglomerate a subset
            subset_idxs = np.random.choice(
                np.arange(0, n), self.max_agglom_size
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

        def _fit_for_data(p):
            n_clusters = p[1]["n_components"]
            if (p[0]["affinity"] != "none") and (self.label_init is None):
                index = param_grid_ag.index(p[0])
                agg_clustering = ag_labels[index][
                    :, n_clusters - self.min_components
                ]
            else:
                agg_clustering = []
            return self._fit_cluster(X, X_subset, y, p, agg_clustering)

        results = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
            delayed(_fit_for_data)(p) for p in param_grid
        )
        best_idx = np.argmin(
            [results[i]["bic/aic"] for i in range(len(results))]
        )

        self.best_criterion_ = results[best_idx]["bic/aic"]
        self.n_components_ = results[best_idx]["n_components"]
        self.covariance_type_ = results[best_idx]["covariance_type"]
        self.affinity_ = results[best_idx]["affinity"]
        self.linkage_ = results[best_idx]["linkage"]
        self.reg_covar_ = results[best_idx]["reg_covar"]
        self.best_model_ = results[best_idx]["model"]
        self.n_iter_ = results[best_idx]["model"].n_iter_
        self.labels_ = results[best_idx]["model"].predict(X)

        results_dict = results[0].copy()
        for key in results_dict:
            for sub in range(1, len(results)):
                results_dict[key] = np.append(
                    results_dict[key], results[sub][key]
                )
        self.results_ = results_dict

        return self

    def predict(self, X):
        """
        Predict clusters based on the best Gaussian mixture model.

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
        labels = self.best_model_.predict(X)

        return labels

    def fit_predict(self, X, y=None):
        """
        Fit the models and predict clusters based on the best model.

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
    """
    Increase regularization factor by factor of 10.

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
        reg = reg * 10
    return reg


def _onehot_to_initial_params(X, onehot, cov_type):
    """
    Computes cluster weights, cluster means and cluster precisions from
    a given clustering.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
        List of n_features-dimensional data points. Each row
        corresponds to a single data point.
    onehot : ndarray, shape (n_samples, n_clusters)
        Each row has a 1 indicating cluster membership, other entries are 0.
    cov_type : {'full', 'tied', 'diag', 'spherical'}
        Covariance type for Gaussian mixture model
    """
    n = X.shape[0]
    weights, means, covariances = _estimate_gaussian_parameters(
        X, onehot, 1e-06, cov_type
    )
    weights /= n

    precisions_cholesky_ = _compute_precision_cholesky(covariances, cov_type)

    if cov_type == "tied":
        c = precisions_cholesky_
        precisions = np.dot(c, c.T)
    elif cov_type == "diag":
        precisions = precisions_cholesky_
    else:
        precisions = [np.dot(c, c.T) for c in precisions_cholesky_]

    return weights, means, precisions


def _labels_to_onehot(labels):
    """
    Converts labels to one-hot format.

    Parameters
    ----------
    labels : ndarray, shape (n_samples,)
        Cluster labels

    Returns
    -------
    onehot : ndarray, shape (n_samples, n_clusters)
        Each row has a single one indicating cluster membership.
        All other entries are zero.
    """
    n = len(labels)
    k = max(labels) + 1
    onehot = np.zeros([n, k])
    onehot[np.arange(n), labels] = 1
    return onehot


def _process_paramgrid(paramgrid):
    """
    Removes combinations of affinity and linkage that are not possible.

    Parameters
    ----------
    paramgrid : list of dicts
        Each dict has the keys 'affinity', 'covariance_type', 'linkage',
        'n_components', and 'random_state'

    Returns
    -------
    paramgrid_processed : list pairs of dicts
        For each pair, the first dict are the options
        for AgglomerativeClustering.
        The second dict include the options for GaussianMixture.
    ag_paramgrid_processed : list of dicts
        options for AgglomerativeClustering
    """
    gm_keys = ["covariance_type", "n_components", "random_state"]
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
            if ag_params not in ag_params_processed:
                ag_params_processed.append(ag_params)

            paramgrid_processed.append([ag_params, gm_params])
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
                (
                    hierarchical_labels,
                    hierarchical_labels[:, -1].reshape((-1, 1)),
                )
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
