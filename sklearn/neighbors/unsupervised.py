"""Unsupervised nearest neighbors learner"""

from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import RadiusNeighborsMixin
from .base import UnsupervisedMixin


class NearestNeighbors(NeighborsBase, KNeighborsMixin,
                       RadiusNeighborsMixin, UnsupervisedMixin):
    """Unsupervised learner for implementing neighbor searches.

    Read more in the :ref:`User Guide <unsupervised_neighbors>`.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`k_neighbors` queries.

    radius : float, optional (default = 1.0)
        Range of parameter space to use by default for :meth`radius_neighbors`
        queries.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDtree`
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

    p: integer, optional (default = 2)
        Parameter for the Minkowski metric from
        sklearn.metrics.pairwise.pairwise_distances. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

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
          'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
          'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
          'sqeuclidean', 'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`k_neighbors` and :meth:`kneighbors_graph` methods.

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
        self._init_params(n_neighbors=n_neighbors,
                          radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size, metric=metric, p=p,
                          metric_params=metric_params, n_jobs=n_jobs, **kwargs)
