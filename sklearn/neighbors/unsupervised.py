"""Unsupervised nearest neighbors learner"""
import numpy as np

from .base import \
    NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin, UnsupervisedMixin

class NearestNeighbors(NeighborsBase, KNeighborsMixin,
                       RadiusNeighborsMixin, UnsupervisedMixin):
    """Unsupervised learner for implementing neighbor searches.

    Parameters
    ----------
    n_neighbors : int, optional
        Number of neighbors to use by default for k_neighbors() queries.
        Default value is 5

    radius : float, optional
        Range of parameter space to use by default for radius_neighbors()
        queries.  Default value is 1.0

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree' will
        construct a BallTree, 'kd_tree' will use the KD-tree implementation
        in `scipy.spatial.ckdtree`, while 'brute' will perform brute-force
        search. 'auto' will guess the most appropriate based on current
        dataset.  See Discussion below for notes on choosing the optimal
        method. Fitting on sparse input will override the setting of this
        parameter.

    leaf_size : int, optional
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the nearest neighbors query.  The optimal value depends on the
        nature of the problem (see Discussion below).  Defaults to 20.

    Examples
    --------
    >>> samples = [[0, 0, 2], [1, 0, 0], [0, 0, 1]]
    >>> from sklearn.neighbors import NearestNeighbors
    >>> neigh = NearestNeighbors(2, 0.4)
    >>> neigh.fit(samples)
    NearestNeighbors(algorithm='auto', leaf_size=20, n_neighbors=2,
             radius=0.4)
    >>> neigh.kneighbors([[0, 0, 1.3]], 2, return_distance=False)
    array([[2, 0]])
    >>> neigh.radius_neighbors([0, 0, 1.3], 0.4, return_distance=False)
    array([[2]], dtype=object)

    See also
    --------
    NeighborsClassifier
    NeighborsRegressor
    BallTree

    Discussion
    ----------
    The optimal algorithm for a given dataset is a complicated choice, and
    depends on a number of factors:
    * number of samples N (n_samples).
      * 'brute' query time grows as O[N], while
      * 'ball_tree' and 'kd_tree' query time grows as O[log(N)]
    * dimensionality D (n_features)
      Intrinsic dimensionality refers to the dimension of a manifold which
      is linearly or nonlinearly embedded within the parameter space
      * 'brute' query time grows as O[D], and is unaffected by the value of d.
      * 'ball_tree' query time may grow faster or slower than this, depending
        on the structure of the data.
    * data structure: intrinsic dimensionality of the data and/or sparsity
      of the data. Intrinsic dimensionality refers to the dimension d<=D
      of a manifold on which the data lies, which can be linearly or
      nonlinearly embedded in the parameter space. Sparsity refers to the
      degree to which the data fills the parameter space (this is to be
      distinguished from the concept as used in "sparse" matrices.  The data
      matrix may have no zero entries, but the structure can still be
      "sparse" in this sense).
      * 'brute' query time is unchanged by data structure.
      * 'ball_tree' query time is greatly influenced by data structure.
        In general, the more sparse the data is, and the smaller the intrinsic
        dimension d, the faster the ball_tree algorithm will be compared to
        a brute-force search.  For data which densely fills the parameter
        space, brute-force is generally a better choice.
    * number of neighbors k requested for a query point.
      * 'brute' query time is unaffected by the value of k
      * 'ball_tree' query time is slower for k>1, mainly due to the internal
        queueing and sorting that takes place during the query.
    * leaf_size of the ball_tree
      The leaf_size parameter controls the point at which the ball_tree
      algorithm switches from a tree-based query to a brute-force search.
      As the number of points 'n' in a node becomes smaller, the difference
      between a brute-force query time O[n] and a ball-tree query time
      O[log(n)] becomes less significant.  At some point, depending on all
      the above-mentioned factors, brute-force becomes more efficient.
      The parameter leaf_size sets this threshold.
    * number of query points.
      The ball_tree algorithm requires a one-time building phase.  For very
      few queries, the time to build the tree may overwhelm any gain during
      the query time.

    Currently, the algorithm='auto' option chooses between 'brute' and
    'ball_tree' using an unsophisticated rubric.  In practice, it is
    suggested that the user experiment with different choices of
    `algorithm` and `leaf_size` to determine the optimal configuration.

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto', leaf_size=20):
        self._init_params(n_neighbors=n_neighbors,
                          radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
