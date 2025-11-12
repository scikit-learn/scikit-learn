"""
Multi-dimensional Scaling (MDS).
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Integral, Real

import numpy as np
from joblib import effective_n_jobs

from sklearn.base import BaseEstimator, _fit_context
from sklearn.isotonic import IsotonicRegression
from sklearn.manifold import ClassicalMDS
from sklearn.metrics import euclidean_distances, pairwise_distances
from sklearn.utils import check_array, check_random_state, check_symmetric
from sklearn.utils._param_validation import (
    Hidden,
    Interval,
    StrOptions,
    validate_params,
)
from sklearn.utils.parallel import Parallel, delayed
from sklearn.utils.validation import validate_data


def _smacof_single(
    dissimilarities,
    metric=True,
    n_components=2,
    init=None,
    max_iter=300,
    verbose=0,
    eps=1e-6,
    random_state=None,
    normalized_stress=False,
):
    """Computes multidimensional scaling using SMACOF algorithm.

    Parameters
    ----------
    dissimilarities : ndarray of shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Must be symmetric.

    metric : bool, default=True
        Compute metric or nonmetric SMACOF algorithm.
        When ``False`` (i.e. non-metric MDS), dissimilarities with 0 are considered as
        missing values.

    n_components : int, default=2
        Number of dimensions in which to immerse the dissimilarities. If an
        ``init`` array is provided, this option is overridden and the shape of
        ``init`` is used to determine the dimensionality of the embedding
        space.

    init : ndarray of shape (n_samples, n_components), default=None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.

    max_iter : int, default=300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, default=0
        Level of verbosity.

    eps : float, default=1e-6
        The tolerance with respect to stress (normalized by the sum of squared
        embedding distances) at which to declare convergence.

        .. versionchanged:: 1.7
           The default value for `eps` has changed from 1e-3 to 1e-6, as a result
           of a bugfix in the computation of the convergence criterion.

    random_state : int, RandomState instance or None, default=None
        Determines the random number generator used to initialize the centers.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    normalized_stress : bool, default=False
        Whether to return normalized stress value (Stress-1) instead of raw
        stress.

        .. versionadded:: 1.2

        .. versionchanged:: 1.7
           Normalized stress is now supported for metric MDS as well.

    Returns
    -------
    X : ndarray of shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).
        If `normalized_stress=True`, returns Stress-1.
        A value of 0 indicates "perfect" fit, 0.025 excellent, 0.05 good,
        0.1 fair, and 0.2 poor [1]_.

    n_iter : int
        The number of iterations corresponding to the best stress.

    References
    ----------
    .. [1] "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
           Psychometrika, 29 (1964)

    .. [2] "Multidimensional scaling by optimizing goodness of fit to a nonmetric
           hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    .. [3] "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
           Groenen P. Springer Series in Statistics (1997)
    """
    dissimilarities = check_symmetric(dissimilarities, raise_exception=True)

    n_samples = dissimilarities.shape[0]
    random_state = check_random_state(random_state)

    dissimilarities_flat = ((1 - np.tri(n_samples)) * dissimilarities).ravel()
    dissimilarities_flat_w = dissimilarities_flat[dissimilarities_flat != 0]
    if init is None:
        # Randomly choose initial configuration
        X = random_state.uniform(size=n_samples * n_components)
        X = X.reshape((n_samples, n_components))
    else:
        # overrides the parameter p
        n_components = init.shape[1]
        if n_samples != init.shape[0]:
            raise ValueError(
                "init matrix should be of shape (%d, %d)" % (n_samples, n_components)
            )
        X = init
    distances = euclidean_distances(X)

    # Out of bounds condition cannot happen because we are transforming
    # the training set here, but does sometimes get triggered in
    # practice due to machine precision issues. Hence "clip".
    ir = IsotonicRegression(out_of_bounds="clip")

    old_stress = None
    for it in range(max_iter):
        # Compute distance and monotonic regression
        if metric:
            disparities = dissimilarities
        else:
            distances_flat = distances.ravel()
            # dissimilarities with 0 are considered as missing values
            distances_flat_w = distances_flat[dissimilarities_flat != 0]

            # Compute the disparities using isotonic regression.
            # For the first SMACOF iteration, use scaled original dissimilarities.
            # (This choice follows the R implementation described in this paper:
            # https://www.jstatsoft.org/article/view/v102i10)
            if it < 1:
                disparities_flat = dissimilarities_flat_w
            else:
                disparities_flat = ir.fit_transform(
                    dissimilarities_flat_w, distances_flat_w
                )
            disparities = np.zeros_like(distances_flat)
            disparities[dissimilarities_flat != 0] = disparities_flat
            disparities = disparities.reshape((n_samples, n_samples))
            disparities *= np.sqrt(
                (n_samples * (n_samples - 1) / 2) / (disparities**2).sum()
            )
            disparities = disparities + disparities.T

        # Update X using the Guttman transform
        distances[distances == 0] = 1e-5
        ratio = disparities / distances
        B = -ratio
        B[np.arange(len(B)), np.arange(len(B))] += ratio.sum(axis=1)
        X = 1.0 / n_samples * np.dot(B, X)

        # Compute stress
        distances = euclidean_distances(X)
        stress = ((distances.ravel() - disparities.ravel()) ** 2).sum() / 2

        if verbose >= 2:  # pragma: no cover
            print(f"Iteration {it}, stress {stress:.4f}")
        if old_stress is not None:
            sum_squared_distances = (distances.ravel() ** 2).sum()
            if ((old_stress - stress) / (sum_squared_distances / 2)) < eps:
                if verbose:  # pragma: no cover
                    print(f"Convergence criterion reached (iteration {it}).")
                break
        old_stress = stress

    if normalized_stress:
        sum_squared_distances = (distances.ravel() ** 2).sum()
        stress = np.sqrt(stress / (sum_squared_distances / 2))

    return X, stress, it + 1


# TODO(1.9): change default `n_init` to 1, see PR #31117
@validate_params(
    {
        "dissimilarities": ["array-like"],
        "metric": ["boolean"],
        "n_components": [Interval(Integral, 1, None, closed="left")],
        "init": ["array-like", None],
        "n_init": [Interval(Integral, 1, None, closed="left"), StrOptions({"warn"})],
        "n_jobs": [Integral, None],
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "verbose": ["verbose"],
        "eps": [Interval(Real, 0, None, closed="left")],
        "random_state": ["random_state"],
        "return_n_iter": ["boolean"],
        "normalized_stress": ["boolean", StrOptions({"auto"})],
    },
    prefer_skip_nested_validation=True,
)
def smacof(
    dissimilarities,
    *,
    metric=True,
    n_components=2,
    init=None,
    n_init="warn",
    n_jobs=None,
    max_iter=300,
    verbose=0,
    eps=1e-6,
    random_state=None,
    return_n_iter=False,
    normalized_stress="auto",
):
    """Compute multidimensional scaling using the SMACOF algorithm.

    The SMACOF (Scaling by MAjorizing a COmplicated Function) algorithm is a
    multidimensional scaling algorithm which minimizes an objective function
    (the *stress*) using a majorization technique. Stress majorization, also
    known as the Guttman Transform, guarantees a monotone convergence of
    stress, and is more powerful than traditional techniques such as gradient
    descent.

    The SMACOF algorithm for metric MDS can be summarized by the following
    steps:

    1. Set an initial start configuration, randomly or not.
    2. Compute the stress
    3. Compute the Guttman Transform
    4. Iterate 2 and 3 until convergence.

    The nonmetric algorithm adds a monotonic regression step before computing
    the stress.

    Parameters
    ----------
    dissimilarities : array-like of shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Must be symmetric.

    metric : bool, default=True
        Compute metric or nonmetric SMACOF algorithm.
        When ``False`` (i.e. non-metric MDS), dissimilarities with 0 are considered as
        missing values.

    n_components : int, default=2
        Number of dimensions in which to immerse the dissimilarities. If an
        ``init`` array is provided, this option is overridden and the shape of
        ``init`` is used to determine the dimensionality of the embedding
        space.

    init : array-like of shape (n_samples, n_components), default=None
        Starting configuration of the embedding to initialize the algorithm. By
        default, the algorithm is initialized with a randomly chosen array.

    n_init : int, default=8
        Number of times the SMACOF algorithm will be run with different
        initializations. The final results will be the best output of the runs,
        determined by the run with the smallest final stress. If ``init`` is
        provided, this option is overridden and a single run is performed.

        .. versionchanged:: 1.9
           The default value for `n_iter` will change from 8 to 1 in version 1.9.

    n_jobs : int, default=None
        The number of jobs to use for the computation. If multiple
        initializations are used (``n_init``), each run of the algorithm is
        computed in parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    max_iter : int, default=300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, default=0
        Level of verbosity.

    eps : float, default=1e-6
        The tolerance with respect to stress (normalized by the sum of squared
        embedding distances) at which to declare convergence.

        .. versionchanged:: 1.7
           The default value for `eps` has changed from 1e-3 to 1e-6, as a result
           of a bugfix in the computation of the convergence criterion.

    random_state : int, RandomState instance or None, default=None
        Determines the random number generator used to initialize the centers.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    return_n_iter : bool, default=False
        Whether or not to return the number of iterations.

    normalized_stress : bool or "auto", default="auto"
        Whether to return normalized stress value (Stress-1) instead of raw
        stress. By default, metric MDS returns raw stress while non-metric MDS
        returns normalized stress.

        .. versionadded:: 1.2

        .. versionchanged:: 1.4
           The default value changed from `False` to `"auto"` in version 1.4.

        .. versionchanged:: 1.7
           Normalized stress is now supported for metric MDS as well.

    Returns
    -------
    X : ndarray of shape (n_samples, n_components)
        Coordinates of the points in a ``n_components``-space.

    stress : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).
        If `normalized_stress=True`, returns Stress-1.
        A value of 0 indicates "perfect" fit, 0.025 excellent, 0.05 good,
        0.1 fair, and 0.2 poor [1]_.

    n_iter : int
        The number of iterations corresponding to the best stress. Returned
        only if ``return_n_iter`` is set to ``True``.

    References
    ----------
    .. [1] "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
           Psychometrika, 29 (1964)

    .. [2] "Multidimensional scaling by optimizing goodness of fit to a nonmetric
           hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    .. [3] "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
           Groenen P. Springer Series in Statistics (1997)

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.manifold import smacof
    >>> from sklearn.metrics import euclidean_distances
    >>> X = np.array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
    >>> dissimilarities = euclidean_distances(X)
    >>> Z, stress = smacof(
    ...     dissimilarities, n_components=2, n_init=1, eps=1e-6, random_state=42
    ... )
    >>> Z.shape
    (3, 2)
    >>> np.round(stress, 6).item()
    3.2e-05
    """

    if n_init == "warn":
        warnings.warn(
            "The default value of `n_init` will change from 8 to 1 in 1.9.",
            FutureWarning,
        )
        n_init = 8

    dissimilarities = check_array(dissimilarities)
    random_state = check_random_state(random_state)

    if normalized_stress == "auto":
        normalized_stress = not metric

    if hasattr(init, "__array__"):
        init = np.asarray(init).copy()
        if not n_init == 1:
            warnings.warn(
                "Explicit initial positions passed: "
                "performing only one init of the MDS instead of %d" % n_init
            )
            n_init = 1

    best_pos, best_stress = None, None

    if effective_n_jobs(n_jobs) == 1:
        for it in range(n_init):
            pos, stress, n_iter_ = _smacof_single(
                dissimilarities,
                metric=metric,
                n_components=n_components,
                init=init,
                max_iter=max_iter,
                verbose=verbose,
                eps=eps,
                random_state=random_state,
                normalized_stress=normalized_stress,
            )
            if best_stress is None or stress < best_stress:
                best_stress = stress
                best_pos = pos.copy()
                best_iter = n_iter_
    else:
        seeds = random_state.randint(np.iinfo(np.int32).max, size=n_init)
        results = Parallel(n_jobs=n_jobs, verbose=max(verbose - 1, 0))(
            delayed(_smacof_single)(
                dissimilarities,
                metric=metric,
                n_components=n_components,
                init=init,
                max_iter=max_iter,
                verbose=verbose,
                eps=eps,
                random_state=seed,
                normalized_stress=normalized_stress,
            )
            for seed in seeds
        )
        positions, stress, n_iters = zip(*results)
        best = np.argmin(stress)
        best_stress = stress[best]
        best_pos = positions[best]
        best_iter = n_iters[best]

    if return_n_iter:
        return best_pos, best_stress, best_iter
    else:
        return best_pos, best_stress


# TODO(1.9): change default `n_init` to 1, see PR #31117
# TODO(1.10): change default `init` to "classical_mds", see PR #32229
# TODO(1.10): drop support for boolean `metric`, see PR #32229
# TODO(1.10): drop support for `dissimilarity`, see PR #32229
class MDS(BaseEstimator):
    """Multidimensional scaling.

    Read more in the :ref:`User Guide <multidimensional_scaling>`.

    Parameters
    ----------
    n_components : int, default=2
        Number of dimensions in which to immerse the dissimilarities.

    metric_mds : bool, default=True
        If ``True``, perform metric MDS; otherwise, perform nonmetric MDS.
        When ``False`` (i.e. non-metric MDS), dissimilarities with 0 are considered as
        missing values.

        .. versionchanged:: 1.8
           The parameter `metric` was renamed into `metric_mds`.

    n_init : int, default=4
        Number of times the SMACOF algorithm will be run with different
        initializations. The final results will be the best output of the runs,
        determined by the run with the smallest final stress.

        .. versionchanged:: 1.9
           The default value for `n_init` will change from 4 to 1 in version 1.9.

    init : {'random', 'classical_mds'}, default='random'
        The initialization approach. If `random`, random initialization is used.
        If `classical_mds`, then classical MDS is run and used as initialization
        for MDS (in this case, the value of `n_init` is ignored).

        .. versionadded:: 1.8

        .. versionchanged:: 1.10
           The default value for `init` will change to `classical_mds`.

    max_iter : int, default=300
        Maximum number of iterations of the SMACOF algorithm for a single run.

    verbose : int, default=0
        Level of verbosity.

    eps : float, default=1e-6
        The tolerance with respect to stress (normalized by the sum of squared
        embedding distances) at which to declare convergence.

        .. versionchanged:: 1.7
           The default value for `eps` has changed from 1e-3 to 1e-6, as a result
           of a bugfix in the computation of the convergence criterion.

    n_jobs : int, default=None
        The number of jobs to use for the computation. If multiple
        initializations are used (``n_init``), each run of the algorithm is
        computed in parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    random_state : int, RandomState instance or None, default=None
        Determines the random number generator used to initialize the centers.
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    dissimilarity : {'euclidean', 'precomputed'}
        Dissimilarity measure to use:

        - 'euclidean':
            Pairwise Euclidean distances between points in the dataset.

        - 'precomputed':
            Pre-computed dissimilarities are passed directly to ``fit`` and
            ``fit_transform``.

        .. deprecated:: 1.8
           `dissimilarity` was renamed `metric` in 1.8 and will be removed in 1.10.

    metric : str or callable, default='euclidean'
        Metric to use for dissimilarity computation. Default is "euclidean".

        If metric is a string, it must be one of the options allowed by
        `scipy.spatial.distance.pdist` for its metric parameter, or a metric
        listed in :func:`sklearn.metrics.pairwise.distance_metrics`

        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square during fit.

        If metric is a callable function, it takes two arrays representing 1D
        vectors as inputs and must return one value indicating the distance
        between those vectors. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        .. versionchanged:: 1.8
           Prior to 1.8, `metric=True/False` was used to select metric/non-metric
           MDS, which is now the role of `metric_mds`.  The support for ``True``
           and ``False`` will be dropped in version 1.10, use `metric_mds` instead.

    metric_params : dict, default=None
        Additional keyword arguments for the dissimilarity computation.

        .. versionadded:: 1.8

    normalized_stress : bool or "auto" default="auto"
        Whether to return normalized stress value (Stress-1) instead of raw
        stress. By default, metric MDS returns raw stress while non-metric MDS
        returns normalized stress.

        .. versionadded:: 1.2

        .. versionchanged:: 1.4
           The default value changed from `False` to `"auto"` in version 1.4.

        .. versionchanged:: 1.7
           Normalized stress is now supported for metric MDS as well.

    Attributes
    ----------
    embedding_ : ndarray of shape (n_samples, n_components)
        Stores the position of the dataset in the embedding space.

    stress_ : float
        The final value of the stress (sum of squared distance of the
        disparities and the distances for all constrained points).
        If `normalized_stress=True`, returns Stress-1.
        A value of 0 indicates "perfect" fit, 0.025 excellent, 0.05 good,
        0.1 fair, and 0.2 poor [1]_.

    dissimilarity_matrix_ : ndarray of shape (n_samples, n_samples)
        Pairwise dissimilarities between the points. Symmetric matrix that:

        - either uses a custom dissimilarity matrix by setting `dissimilarity`
          to 'precomputed';
        - or constructs a dissimilarity matrix from data using
          Euclidean distances.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_iter_ : int
        The number of iterations corresponding to the best stress.

    See Also
    --------
    sklearn.decomposition.PCA : Principal component analysis that is a linear
        dimensionality reduction method.
    sklearn.decomposition.KernelPCA : Non-linear dimensionality reduction using
        kernels and PCA.
    TSNE : T-distributed Stochastic Neighbor Embedding.
    Isomap : Manifold learning based on Isometric Mapping.
    LocallyLinearEmbedding : Manifold learning using Locally Linear Embedding.
    SpectralEmbedding : Spectral embedding for non-linear dimensionality.

    References
    ----------
    .. [1] "Nonmetric multidimensional scaling: a numerical method" Kruskal, J.
       Psychometrika, 29 (1964)

    .. [2] "Multidimensional scaling by optimizing goodness of fit to a nonmetric
       hypothesis" Kruskal, J. Psychometrika, 29, (1964)

    .. [3] "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
       Groenen P. Springer Series in Statistics (1997)

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.manifold import MDS
    >>> X, _ = load_digits(return_X_y=True)
    >>> X.shape
    (1797, 64)
    >>> embedding = MDS(n_components=2, n_init=1, init="random")
    >>> X_transformed = embedding.fit_transform(X[:100])
    >>> X_transformed.shape
    (100, 2)

    For a more detailed example of usage, see
    :ref:`sphx_glr_auto_examples_manifold_plot_mds.py`.

    For a comparison of manifold learning techniques, see
    :ref:`sphx_glr_auto_examples_manifold_plot_compare_methods.py`.
    """

    _parameter_constraints: dict = {
        "n_components": [Interval(Integral, 1, None, closed="left")],
        "metric_mds": ["boolean"],
        "n_init": [
            Interval(Integral, 1, None, closed="left"),
            Hidden(StrOptions({"warn"})),
        ],
        "init": [StrOptions({"random", "classical_mds"}), Hidden(StrOptions({"warn"}))],
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "verbose": ["verbose"],
        "eps": [Interval(Real, 0.0, None, closed="left")],
        "n_jobs": [None, Integral],
        "random_state": ["random_state"],
        "dissimilarity": [
            StrOptions({"euclidean", "precomputed"}),
            Hidden(StrOptions({"deprecated"})),
        ],
        "metric": [str, callable, Hidden("boolean")],
        "metric_params": [dict, None],
        "normalized_stress": ["boolean", StrOptions({"auto"})],
    }

    def __init__(
        self,
        n_components=2,
        *,
        metric_mds=True,
        n_init="warn",
        init="warn",
        max_iter=300,
        verbose=0,
        eps=1e-6,
        n_jobs=None,
        random_state=None,
        dissimilarity="deprecated",
        metric="euclidean",
        metric_params=None,
        normalized_stress="auto",
    ):
        self.n_components = n_components
        self.dissimilarity = dissimilarity
        self.metric = metric
        self.metric_params = metric_params
        self.metric_mds = metric_mds
        self.n_init = n_init
        self.init = init
        self.max_iter = max_iter
        self.eps = eps
        self.verbose = verbose
        self.n_jobs = n_jobs
        self.random_state = random_state
        self.normalized_stress = normalized_stress

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.pairwise = (self.dissimilarity == "precomputed") | (
            self.metric == "precomputed"
        )
        return tags

    def fit(self, X, y=None, init=None):
        """
        Compute the position of the points in the embedding space.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or \
                (n_samples, n_samples)
            Input data. If ``metric=='precomputed'``, the input should
            be the dissimilarity matrix.

        y : Ignored
            Not used, present for API consistency by convention.

        init : ndarray of shape (n_samples, n_components), default=None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        self.fit_transform(X, init=init)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, X, y=None, init=None):
        """
        Fit the data from `X`, and returns the embedded coordinates.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or \
                (n_samples, n_samples)
            Input data. If ``metric=='precomputed'``, the input should
            be the dissimilarity matrix.

        y : Ignored
            Not used, present for API consistency by convention.

        init : ndarray of shape (n_samples, n_components), default=None
            Starting configuration of the embedding to initialize the SMACOF
            algorithm. By default, the algorithm is initialized with a randomly
            chosen array.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            X transformed in the new space.
        """

        if self.n_init == "warn":
            warnings.warn(
                "The default value of `n_init` will change from 4 to 1 in 1.9. "
                "To suppress this warning, provide some value of `n_init`.",
                FutureWarning,
            )
            self._n_init = 4
        else:
            self._n_init = self.n_init

        if self.init == "warn":
            warnings.warn(
                "The default value of `init` will change from 'random' to "
                "'classical_mds' in 1.10. To suppress this warning, provide "
                "some value of `init`.",
                FutureWarning,
            )
            self._init = "random"
        else:
            self._init = self.init

        if self.dissimilarity != "deprecated":
            if not isinstance(self.metric, bool) and self.metric != "euclidean":
                raise ValueError(
                    "You provided both `dissimilarity` and `metric`. Please use "
                    "only `metric`."
                )
            else:
                warnings.warn(
                    "The `dissimilarity` parameter is deprecated and will be "
                    "removed in 1.10. Use `metric` instead.",
                    FutureWarning,
                )
                self._metric = self.dissimilarity

        if isinstance(self.metric, bool):
            warnings.warn(
                f"Use metric_mds={self.metric} instead of metric={self.metric}. The "
                "support for metric={True/False} will be dropped in 1.10.",
                FutureWarning,
            )
            if self.dissimilarity == "deprecated":
                self._metric = "euclidean"
            self._metric_mds = self.metric
        else:
            if self.dissimilarity == "deprecated":
                self._metric = self.metric
            self._metric_mds = self.metric_mds

        X = validate_data(self, X)
        if X.shape[0] == X.shape[1] and self._metric != "precomputed":
            warnings.warn(
                "The provided input is a square matrix. Note that ``fit`` constructs "
                "a dissimilarity matrix from data and will treat rows as samples "
                "and columns as features. To use a pre-computed dissimilarity matrix, "
                "set ``metric='precomputed'``."
            )

        if self._metric == "precomputed":
            self.dissimilarity_matrix_ = X
            self.dissimilarity_matrix_ = check_symmetric(
                self.dissimilarity_matrix_, raise_exception=True
            )
        else:
            self.dissimilarity_matrix_ = pairwise_distances(
                X,
                metric=self._metric,
                **(self.metric_params if self.metric_params is not None else {}),
            )

        if init is not None:
            init_array = init
        elif self._init == "classical_mds":
            cmds = ClassicalMDS(metric="precomputed")
            init_array = cmds.fit_transform(self.dissimilarity_matrix_)
        else:
            init_array = None

        self.embedding_, self.stress_, self.n_iter_ = smacof(
            self.dissimilarity_matrix_,
            metric=self._metric_mds,
            n_components=self.n_components,
            init=init_array,
            n_init=self._n_init,
            n_jobs=self.n_jobs,
            max_iter=self.max_iter,
            verbose=self.verbose,
            eps=self.eps,
            random_state=self.random_state,
            return_n_iter=True,
            normalized_stress=self.normalized_stress,
        )

        return self.embedding_
