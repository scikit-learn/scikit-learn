"""
DBSCAN: Density-Based Spatial Clustering of Applications with Noise
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from numbers import Integral, Number, Real

import numpy as np
from scipy import sparse

from ..base import BaseEstimator, ClusterMixin, _fit_context
from ..metrics.pairwise import _VALID_METRICS
from ..neighbors import NearestNeighbors
from ..utils import check_random_state
from ..utils._param_validation import Interval, StrOptions, validate_params
from ..utils.validation import _check_sample_weight, validate_data
from ._dbscan_inner import dbscan_inner


@validate_params(
    {
        "X": ["array-like", "sparse matrix"],
        "sample_weight": ["array-like", None],
        "subsample": [Interval(Real, 0, 1, closed="neither"), None],
        "random_state": ["random_state"],
    },
    prefer_skip_nested_validation=False,
)
def dbscan(
    X,
    eps=0.5,
    *,
    min_samples=5,
    metric="minkowski",
    metric_params=None,
    algorithm="auto",
    leaf_size=30,
    p=2,
    sample_weight=None,
    n_jobs=None,
    subsample=None,
    random_state=None,
):
    """Perform DBSCAN clustering from vector array or distance matrix.

    Read more in the :ref:`User Guide <dbscan>`.

    Parameters
    ----------
    X : {array-like, sparse (CSR) matrix} of shape (n_samples, n_features) or \
            (n_samples, n_samples)
        A feature array, or array of distances between samples if
        ``metric='precomputed'``.

    eps : float, default=0.5
        The maximum distance between two samples for one to be considered
        as in the neighborhood of the other. This is not a maximum bound
        on the distances of points within a cluster. This is the most
        important DBSCAN parameter to choose appropriately for your data set
        and distance function.

    min_samples : int, default=5
        The number of samples (or total weight) in a neighborhood for a point
        to be considered as a core point. This includes the point itself.

    metric : str or callable, default='minkowski'
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by :func:`sklearn.metrics.pairwise_distances` for
        its metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square during fit.
        X may be a :term:`sparse graph <sparse graph>`,
        in which case only "nonzero" elements may be considered neighbors.

    metric_params : dict, default=None
        Additional keyword arguments for the metric function.

        .. versionadded:: 0.19

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, default='auto'
        The algorithm to be used by the NearestNeighbors module
        to compute pointwise distances and find nearest neighbors.
        See NearestNeighbors module documentation for details.

    leaf_size : int, default=30
        Leaf size passed to BallTree or cKDTree. This can affect the speed
        of the construction and query, as well as the memory required
        to store the tree. The optimal value depends
        on the nature of the problem.

    p : float, default=2
        The power of the Minkowski metric to be used to calculate distance
        between points.

    sample_weight : array-like of shape (n_samples,), default=None
        Weight of each sample, such that a sample with a weight of at least
        ``min_samples`` is by itself a core sample; a sample with negative
        weight may inhibit its eps-neighbor from being core.
        Note that weights are absolute, and default to 1.

    n_jobs : int, default=None
        The number of parallel jobs to run for neighbors search. ``None`` means
        1 unless in a :obj:`joblib.parallel_backend` context. ``-1`` means
        using all processors. See :term:`Glossary <n_jobs>` for more details.
        If precomputed distance are used, parallel execution is not available
        and thus n_jobs will have no effect.

    subsample : float, default=None
        Should be between [0, 1]. By default, no sampling is done.
        Sampling probability, representing the proportion of the dataset
        that can be labeled a core sample. The lower the subsample, the
        less memory and computation is used.
        See: Jang, J. and Jiang, H. "DBSCAN++: Towards fast and scalable
        density clustering". Proceedings of the 36th International Conference
        on Machine Learning, 2019.

    random_state : int, RandomState instance or None, default=None
        Only relevant when ``subsample`` is set. Controls the randomness
        of the subsampling. Pass an int for reproducible output across
        multiple function calls. See :term:`Glossary <random_state>`.

    Returns
    -------
    core_samples : ndarray of shape (n_core_samples,)
        Indices of core samples.

    labels : ndarray of shape (n_samples,)
        Cluster labels for each point.  Noisy samples are given the label -1.

    See Also
    --------
    DBSCAN : An estimator interface for this clustering algorithm.
    OPTICS : A similar estimator interface clustering at multiple values of
        eps. Our implementation is optimized for memory usage.

    Notes
    -----
    For an example, see :ref:`sphx_glr_auto_examples_cluster_plot_dbscan.py`.

    This implementation bulk-computes all neighborhood queries, which increases
    the memory complexity to O(n.d) where d is the average number of neighbors,
    while original DBSCAN had memory complexity O(n). It may attract a higher
    memory complexity when querying these nearest neighborhoods, depending
    on the ``algorithm``.

    One way to avoid the query complexity is to pre-compute sparse
    neighborhoods in chunks using
    :func:`NearestNeighbors.radius_neighbors_graph
    <sklearn.neighbors.NearestNeighbors.radius_neighbors_graph>` with
    ``mode='distance'``, then using ``metric='precomputed'`` here.

    Another way to reduce memory and computation time is to remove
    (near-)duplicate points and use ``sample_weight`` instead.

    :class:`~sklearn.cluster.OPTICS` provides a similar clustering with lower
    memory usage.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, `"A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise"
    <https://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD-96.final.frame.pdf>`_.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226-231. 1996

    Schubert, E., Sander, J., Ester, M., Kriegel, H. P., & Xu, X. (2017).
    :doi:`"DBSCAN revisited, revisited: why and how you should (still) use DBSCAN."
    <10.1145/3068335>`
    ACM Transactions on Database Systems (TODS), 42(3), 19.

    Examples
    --------
    >>> from sklearn.cluster import dbscan
    >>> X = [[1, 2], [2, 2], [2, 3], [8, 7], [8, 8], [25, 80]]
    >>> core_samples, labels = dbscan(X, eps=3, min_samples=2)
    >>> core_samples
    array([0, 1, 2, 3, 4])
    >>> labels
    array([ 0,  0,  0,  1,  1, -1])
    """

    est = DBSCAN(
        eps=eps,
        min_samples=min_samples,
        metric=metric,
        metric_params=metric_params,
        algorithm=algorithm,
        leaf_size=leaf_size,
        p=p,
        n_jobs=n_jobs,
    )
    est.fit(
        X, sample_weight=sample_weight, subsample=subsample, random_state=random_state
    )
    return est.core_sample_indices_, est.labels_


class DBSCAN(ClusterMixin, BaseEstimator):
    """Perform DBSCAN clustering from vector array or distance matrix.

    DBSCAN - Density-Based Spatial Clustering of Applications with Noise.
    Finds core samples of high density and expands clusters from them.
    Good for data which contains clusters of similar density.

    This implementation has a worst case memory complexity of :math:`O({n}^2)`,
    which can occur when the `eps` param is large and `min_samples` is low,
    while the original DBSCAN only uses linear memory.
    For further details, see the Notes below.

    Read more in the :ref:`User Guide <dbscan>`.

    Parameters
    ----------
    eps : float, default=0.5
        The maximum distance between two samples for one to be considered
        as in the neighborhood of the other. This is not a maximum bound
        on the distances of points within a cluster. This is the most
        important DBSCAN parameter to choose appropriately for your data set
        and distance function.

    min_samples : int, default=5
        The number of samples (or total weight) in a neighborhood for a point to
        be considered as a core point. This includes the point itself. If
        `min_samples` is set to a higher value, DBSCAN will find denser clusters,
        whereas if it is set to a lower value, the found clusters will be more
        sparse.

    metric : str, or callable, default='euclidean'
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string or callable, it must be one of
        the options allowed by :func:`sklearn.metrics.pairwise_distances` for
        its metric parameter.
        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square. X may be a :term:`sparse graph`, in which
        case only "nonzero" elements may be considered neighbors for DBSCAN.

        .. versionadded:: 0.17
           metric *precomputed* to accept precomputed sparse matrix.

    metric_params : dict, default=None
        Additional keyword arguments for the metric function.

        .. versionadded:: 0.19

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, default='auto'
        The algorithm to be used by the NearestNeighbors module
        to compute pointwise distances and find nearest neighbors.
        See NearestNeighbors module documentation for details.

    leaf_size : int, default=30
        Leaf size passed to BallTree or cKDTree. This can affect the speed
        of the construction and query, as well as the memory required
        to store the tree. The optimal value depends
        on the nature of the problem.

    p : float, default=None
        The power of the Minkowski metric to be used to calculate distance
        between points. If None, then ``p=2`` (equivalent to the Euclidean
        distance).

    n_jobs : int, default=None
        The number of parallel jobs to run.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    core_sample_indices_ : ndarray of shape (n_core_samples,)
        Indices of core samples.

    components_ : ndarray of shape (n_core_samples, n_features)
        Copy of each core sample found by training.

    labels_ : ndarray of shape (n_samples)
        Cluster labels for each point in the dataset given to fit().
        Noisy samples are given the label -1.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    See Also
    --------
    OPTICS : A similar clustering at multiple values of eps. Our implementation
        is optimized for memory usage.

    Notes
    -----
    This implementation bulk-computes all neighborhood queries, which increases
    the memory complexity to O(n.d) where d is the average number of neighbors,
    while original DBSCAN had memory complexity O(n). It may attract a higher
    memory complexity when querying these nearest neighborhoods, depending
    on the ``algorithm``.

    One way to avoid the query complexity is to pre-compute sparse
    neighborhoods in chunks using
    :func:`NearestNeighbors.radius_neighbors_graph
    <sklearn.neighbors.NearestNeighbors.radius_neighbors_graph>` with
    ``mode='distance'``, then using ``metric='precomputed'`` here.

    Another way to reduce memory and computation time is to remove
    (near-)duplicate points and use ``sample_weight`` instead.

    Yet another way is to use ``subsample`` in order to reduce the core
    samples search space.

    :class:`~sklearn.cluster.OPTICS` provides a similar clustering with lower memory
    usage.

    References
    ----------
    Ester, M., H. P. Kriegel, J. Sander, and X. Xu, `"A Density-Based
    Algorithm for Discovering Clusters in Large Spatial Databases with Noise"
    <https://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD-96.final.frame.pdf>`_.
    In: Proceedings of the 2nd International Conference on Knowledge Discovery
    and Data Mining, Portland, OR, AAAI Press, pp. 226-231. 1996

    Schubert, E., Sander, J., Ester, M., Kriegel, H. P., & Xu, X. (2017).
    :doi:`"DBSCAN revisited, revisited: why and how you should (still) use DBSCAN."
    <10.1145/3068335>`
    ACM Transactions on Database Systems (TODS), 42(3), 19.

    Examples
    --------
    >>> from sklearn.cluster import DBSCAN
    >>> import numpy as np
    >>> X = np.array([[1, 2], [2, 2], [2, 3],
    ...               [8, 7], [8, 8], [25, 80]])
    >>> clustering = DBSCAN(eps=3, min_samples=2).fit(X)
    >>> clustering.labels_
    array([ 0,  0,  0,  1,  1, -1])
    >>> clustering
    DBSCAN(eps=3, min_samples=2)

    For an example, see
    :ref:`sphx_glr_auto_examples_cluster_plot_dbscan.py`.

    For a comparison of DBSCAN with other clustering algorithms, see
    :ref:`sphx_glr_auto_examples_cluster_plot_cluster_comparison.py`
    """

    _parameter_constraints: dict = {
        "eps": [Interval(Real, 0.0, None, closed="neither")],
        "min_samples": [Interval(Integral, 1, None, closed="left")],
        "metric": [
            StrOptions(set(_VALID_METRICS) | {"precomputed"}),
            callable,
        ],
        "metric_params": [dict, None],
        "algorithm": [StrOptions({"auto", "ball_tree", "kd_tree", "brute"})],
        "leaf_size": [Interval(Integral, 1, None, closed="left")],
        "p": [Interval(Real, 0.0, None, closed="left"), None],
        "n_jobs": [Integral, None],
    }

    def __init__(
        self,
        eps=0.5,
        *,
        min_samples=5,
        metric="euclidean",
        metric_params=None,
        algorithm="auto",
        leaf_size=30,
        p=None,
        n_jobs=None,
    ):
        self.eps = eps
        self.min_samples = min_samples
        self.metric = metric
        self.metric_params = metric_params
        self.algorithm = algorithm
        self.leaf_size = leaf_size
        self.p = p
        self.n_jobs = n_jobs

    @_fit_context(
        # DBSCAN.metric is not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y=None, sample_weight=None, subsample=None, random_state=None):
        """Perform DBSCAN clustering from features, or distance matrix.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), or \
            (n_samples, n_samples)
            Training instances to cluster, or distances between instances if
            ``metric='precomputed'``. If a sparse matrix is provided, it will
            be converted into a sparse ``csr_matrix``.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            Weight of each sample, such that a sample with a weight of at least
            ``min_samples`` is by itself a core sample; a sample with a
            negative weight may inhibit its eps-neighbor from being core.
            Note that weights are absolute, and default to 1.

        subsample : float, default=None
            Should be between [0, 1]. By default, no sampling is done.
            Sampling probability, representing the proportion of the dataset
            that can be labeled a core sample. The lower the subsample, the
            less memory and computation is used.
            See: Jang, J. and Jiang, H. "DBSCAN++: Towards fast and scalable
            density clustering". Proceedings of the 36th International Conference
            on Machine Learning, 2019.

        random_state : int, RandomState instance or None, default=None
            Only relevant when ``subsample`` is set. Controls the randomness
            of the subsampling. Pass an int for reproducible output across
            multiple function calls. See :term:`Glossary <random_state>`.

        Returns
        -------
        self : object
            Returns a fitted instance of self.
        """
        X = validate_data(self, X, accept_sparse="csr")

        if sample_weight is not None:
            sample_weight = _check_sample_weight(sample_weight, X)
        if subsample is not None:
            if not isinstance(subsample, Number) or not 0 < subsample < 1:
                raise ValueError("Subsample needs to be float between 0 and 1.")

        # Calculate neighborhood for all samples. This leaves the original
        # point in, which needs to be considered later (i.e. point i is in the
        # neighborhood of point i. While True, its useless information)
        if self.metric == "precomputed" and sparse.issparse(X):
            # set the diagonal to explicit values, as a point is its own
            # neighbor
            X = X.copy()  # copy to avoid in-place modification
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", sparse.SparseEfficiencyWarning)
                X.setdiag(X.diagonal())

        neighbors_model = NearestNeighbors(
            radius=self.eps,
            algorithm=self.algorithm,
            leaf_size=self.leaf_size,
            metric=self.metric,
            metric_params=self.metric_params,
            p=self.p,
            n_jobs=self.n_jobs,
        )

        n = X.shape[0]
        neighbors_model.fit(X)

        if subsample:
            rng = check_random_state(random_state)
            mask = np.full(n, False)
            mask[: int(n * subsample)] = True
            rng.shuffle(mask)
            neighborhoods = np.full(n, None)
            neighborhoods[mask] = neighbors_model.radius_neighbors(
                X[mask], return_distance=False
            )
        else:
            # This has worst case O(n^2) memory complexity
            neighborhoods = neighbors_model.radius_neighbors(X, return_distance=False)

        if sample_weight is None:
            n_neighbors = np.array(
                [
                    0 if neighbors is None else len(neighbors)
                    for neighbors in neighborhoods
                ]
            )
        else:
            n_neighbors = np.array(
                [
                    0 if neighbors is None else np.sum(sample_weight[neighbors])
                    for neighbors in neighborhoods
                ]
            )

        # Initially, all samples are noise.
        labels = np.full(n, -1, dtype=np.intp)

        # A list of all core samples found.
        core_samples = np.asarray(n_neighbors >= self.min_samples, dtype=np.uint8)
        dbscan_inner(core_samples, neighborhoods, labels)

        self.core_sample_indices_ = np.where(core_samples)[0]
        self.labels_ = labels

        if len(self.core_sample_indices_):
            # fix for scipy sparse indexing issue
            self.components_ = X[self.core_sample_indices_].copy()
        else:
            # no core samples
            self.components_ = np.empty((0, X.shape[1]))
        return self

    def fit_predict(
        self, X, y=None, sample_weight=None, subsample=None, random_state=None
    ):
        """Compute clusters from a data or distance matrix and predict labels.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), or \
            (n_samples, n_samples)
            Training instances to cluster, or distances between instances if
            ``metric='precomputed'``. If a sparse matrix is provided, it will
            be converted into a sparse ``csr_matrix``.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            Weight of each sample, such that a sample with a weight of at least
            ``min_samples`` is by itself a core sample; a sample with a
            negative weight may inhibit its eps-neighbor from being core.
            Note that weights are absolute, and default to 1.

        subsample : float, default=None
            Should be between [0, 1]. By default, no sampling is done.
            Sampling probability, representing the proportion of the dataset
            that can be labeled a core sample. The lower the subsample, the
            less memory and computation is used.
            See: Jang, J. and Jiang, H. "DBSCAN++: Towards fast and scalable
            density clustering". Proceedings of the 36th International Conference
            on Machine Learning, 2019.

        random_state : int, RandomState instance or None, default=None
            Only relevant when ``subsample`` is set. Controls the randomness
            of the subsampling. Pass an int for reproducible output across
            multiple function calls. See :term:`Glossary <random_state>`.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Cluster labels. Noisy samples are given the label -1.
        """
        self.fit(
            X,
            sample_weight=sample_weight,
            subsample=subsample,
            random_state=random_state,
        )
        return self.labels_

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.pairwise = self.metric == "precomputed"
        tags.input_tags.sparse = True
        return tags
