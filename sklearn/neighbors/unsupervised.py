"""Unsupervised nearest neighbors learner"""
import numpy as np

from .base import \
    NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin, UnsupervisedMixin


class NearestNeighbors(NeighborsBase, KNeighborsMixin,
                       RadiusNeighborsMixin, UnsupervisedMixin):
    """Unsupervised learner for implementing neighbor searches.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`k_neighbors` queries.

    radius : float, optional (default = 1.0)
        Range of parameter space to use by default for :meth`radius_neighbors`
        queries.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors.
        ``'ball_tree'`` will use :class:`BallTree`,
        ``'kd_tree'`` will use :class:`scipy.spatial.cKDtree`, 
        and ``'brute'`` will use a brute-force search.
        ``'auto'`` will guess the most appropriate algorithm based on the
        values passed to :meth:`fit` method.
        Note: Fitting on sparse input will override the setting of this
        parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the construction and query, as well as the memory required to
        store the tree.  The optimal value depends on the nature of the
        problem.

    Examples
    --------
    >>> samples = [[0, 0, 2], [1, 0, 0], [0, 0, 1]]
    >>> from sklearn.neighbors import NearestNeighbors
    >>> neigh = NearestNeighbors(2, 0.4)
    >>> neigh.fit(samples)
    NearestNeighbors(algorithm='auto', leaf_size=30, n_neighbors=2, radius=0.4)
    >>> neigh.kneighbors([[0, 0, 1.3]], 2, return_distance=False)
    array([[2, 0]])
    >>> neigh.radius_neighbors([0, 0, 1.3], 0.4, return_distance=False)
    array([[2]], dtype=object)

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

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto', leaf_size=30):
        self._init_params(n_neighbors=n_neighbors,
                          radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
