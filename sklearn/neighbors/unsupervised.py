"""Unsupervised nearest neighbors learner"""
from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import RadiusNeighborsMixin
from .base import UnsupervisedMixin
from ..base import TransformerMixin
from ..utils.validation import check_is_fitted, check_array


class NearestNeighbors(NeighborsBase, KNeighborsMixin,
                       RadiusNeighborsMixin, UnsupervisedMixin):
    """Unsupervised learner for implementing neighbor searches.

    Read more in the :ref:`User Guide <unsupervised_neighbors>`.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`kneighbors` queries.

    radius : float, optional (default = 1.0)
        Range of parameter space to use by default for :meth:`radius_neighbors`
        queries.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    metric : string or callable, default 'minkowski'
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    p : integer, optional (default = 2)
        Parameter for the Minkowski metric from
        sklearn.metrics.pairwise.pairwise_distances. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.

    Examples
    --------
      >>> import numpy as np
      >>> from sklearn.neighbors import NearestNeighbors
      >>> samples = [[0, 0, 2], [1, 0, 0], [0, 0, 1]]

      >>> neigh = NearestNeighbors(2, 0.4)
      >>> neigh.fit(samples)  #doctest: +ELLIPSIS
      NearestNeighbors(...)

      >>> neigh.kneighbors([[0, 0, 1.3]], 2, return_distance=False)
      ... #doctest: +ELLIPSIS
      array([[2, 0]]...)

      >>> nbrs = neigh.radius_neighbors([[0, 0, 1.3]], 0.4, return_distance=False)
      >>> np.asarray(nbrs[0][0])
      array(2)

    See also
    --------
    KNeighborsClassifier
    RadiusNeighborsClassifier
    KNeighborsRegressor
    RadiusNeighborsRegressor
    BallTree

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    https://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto', leaf_size=30, metric='minkowski',
                 p=2, metric_params=None, n_jobs=1, **kwargs):
        super(NearestNeighbors, self).__init__(
              n_neighbors=n_neighbors,
              radius=radius,
              algorithm=algorithm,
              leaf_size=leaf_size, metric=metric, p=p,
              metric_params=metric_params, n_jobs=n_jobs, **kwargs)


class KNeighborsTransformer(NeighborsBase, KNeighborsMixin,
                            UnsupervisedMixin, TransformerMixin):
    """Transform X into a (weighted) graph of k nearest neighbors

    The transformed data is a sparse graph as return by kneighbors_graph.

    Parameters
    ----------
    mode : {'distance', 'connectivity'}, optional (default = 'connectivity')
        Type of returned matrix: 'connectivity' will return the connectivity
        matrix with ones and zeros, and 'distance' will return the distances
        between neighbors according to the given metric.

    include_self : bool, default=True.
        Whether or not to mark each sample as the first nearest neighbor to
        itself. If None, then True is used for mode='connectivity' and False
        for mode='distance'.

    n_neighbors : int, optional (default = 5)
        Number of neighbors to use for :meth:`kneighbors` queries.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    metric : string or callable, default 'minkowski'
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    p : integer, optional (default = 2)
        Parameter for the Minkowski metric from
        sklearn.metrics.pairwise.pairwise_distances. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.

    Examples
    --------
    >>> from sklearn.manifold import Isomap
    >>> from sklearn.neighbors import KNeighborsTransformer
    >>> from sklearn.pipeline import make_pipeline
    >>> estimator = make_pipeline(
    ...     KNeighborsTransformer(n_neighbors=5, mode='distance'),
    ...     Isomap(neighbors_algorithm='precomputed'))
    """
    def __init__(self, mode='connectivity', include_self=None,
                 n_neighbors=5, algorithm='auto', leaf_size=30,
                 metric='minkowski', p=2, metric_params=None, n_jobs=1):
        super(KNeighborsTransformer, self).__init__(
            n_neighbors=n_neighbors, radius=None, algorithm=algorithm,
            leaf_size=leaf_size, metric=metric, p=p,
            metric_params=metric_params, n_jobs=n_jobs)
        self.mode = mode
        self.include_self = include_self

    def transform(self, X):
        """Computes the (weighted) graph of Neighbors for points in X

        Parameters
        ----------
        X : array-like, shape = [n_samples_transform, n_features]
            Sample data

        Returns
        -------
        Xt : CSR sparse matrix, shape = [n_samples_fit, n_samples_transform]
            Xt[i, j] is assigned the weight of edge that connects i to j.
        """
        check_is_fitted(self, '_fit_X')
        check_array(X, accept_sparse='csr')

        if self.include_self is None:
            include_self = self.mode == 'connectivity'
        else:
            include_self = self.include_self

        # If we don't include each sample as its own neighbors
        if not include_self and X is self._fit_X:
            X = None

        return self.kneighbors_graph(X, self.n_neighbors, self.mode)


class RadiusNeighborsTransformer(NeighborsBase, RadiusNeighborsMixin,
                                 UnsupervisedMixin, TransformerMixin):
    """Transform X into a (weighted) graph of radius nearest neighbors

    The transformed data is a sparse graph as return by radius_neighbors_graph.

    Parameters
    ----------
    mode : {'distance', 'connectivity'}, optional (default = 'connectivity')
        Type of returned matrix: 'connectivity' will return the connectivity
        matrix with ones and zeros, and 'distance' will return the distances
        between neighbors according to the given metric.

    include_self : bool, default=True.
        Whether or not to mark each sample as the first nearest neighbor to
        itself. If None, then True is used for mode='connectivity' and False
        for mode='distance'.

    radius : float, optional (default = 1.)
        Range of parameter space to use for :meth:`radius_neighbors`

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    metric : string or callable, default 'minkowski'
        metric to use for distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    p : integer, optional (default = 2)
        Parameter for the Minkowski metric from
        sklearn.metrics.pairwise.pairwise_distances. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.

    Examples
    --------
    >>> from sklearn.cluster import DBSCAN
    >>> from sklearn.neighbors import RadiusNeighborsTransformer
    >>> from sklearn.pipeline import make_pipeline
    >>> estimator = make_pipeline(
    ...     RadiusNeighborsTransformer(radius=42.0, mode='distance'),
    ...     DBSCAN(min_samples=30, metric='precomputed'))
    """
    def __init__(self, mode='connectivity', include_self=None,
                 radius=1., algorithm='auto', leaf_size=30,
                 metric='minkowski', p=2, metric_params=None, n_jobs=1):
        super(RadiusNeighborsTransformer, self).__init__(
            n_neighbors=None, radius=radius, algorithm=algorithm,
            leaf_size=leaf_size, metric=metric, p=p,
            metric_params=metric_params, n_jobs=n_jobs)
        self.mode = mode
        self.include_self = include_self

    def transform(self, X):
        """Computes the (weighted) graph of Neighbors for points in X

        Parameters
        ----------
        X : array-like, shape = [n_samples_transform, n_features]
            Sample data

        Returns
        -------
        Xt : CSR sparse matrix, shape = [n_samples_fit, n_samples_transform]
            Xt[i, j] is assigned the weight of edge that connects i to j.
        """
        check_is_fitted(self, '_fit_X')
        check_array(X, accept_sparse='csr')

        if self.include_self is None:
            include_self = self.mode == 'connectivity'
        else:
            include_self = self.include_self

        # If we don't include each sample as its own neighbors
        if not include_self and X is self._fit_X:
            X = None

        return self.radius_neighbors_graph(X, self.radius, self.mode)
