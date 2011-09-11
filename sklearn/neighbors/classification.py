"""Nearest Neighbor Classification"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD, (C) INRIA, University of Amsterdam

import numpy as np
from scipy import stats

from .base import \
    NeighborsBase, KNeighborsMixin,\
    RadiusNeighborsMixin, SupervisedIntegerMixin
from ..base import ClassifierMixin
from ..utils import atleast2d_or_csr, deprecated


class KNeighborsClassifier(NeighborsBase, KNeighborsMixin,
                           SupervisedIntegerMixin, ClassifierMixin):
    """Classifier implementing the k-nearest neighbors vote.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`k_neighbors` queries.

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
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> neigh = KNeighborsClassifier(n_neighbors=2)
    >>> neigh.fit(X, y)
    KNeighborsClassifier(algorithm='auto', leaf_size=30, n_neighbors=2)
    >>> print neigh.predict([[1.5]])
    [0]

    See also
    --------
    RadiusNeighborsClassifier
    KNeighborsRegressor
    RadiusNeighborsRegressor
    NearestNeighbors

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """
    def __init__(self, n_neighbors=5, algorithm='auto', leaf_size=30):
        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size)

    def predict(self, X):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X: array
            A 2-D array representing the test points.

        Returns
        -------
        labels: array
            List of class labels (one for each data sample).
        """
        X = atleast2d_or_csr(X)

        neigh_ind = self.kneighbors(X, return_distance=False)
        pred_labels = self._y[neigh_ind]
        mode, _ = stats.mode(pred_labels, axis=1)
        return mode.flatten().astype(np.int)


class RadiusNeighborsClassifier(NeighborsBase, RadiusNeighborsMixin,
                                SupervisedIntegerMixin, ClassifierMixin):
    """Classifier implementing a vote among neighbors within a given radius

    Parameters
    ----------
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
        Note: fitting on sparse input will override the setting of this
        parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the construction and query, as well as the memory required to
        store the tree.  The optimal value depends on the nature of the
        problem.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import RadiusNeighborsClassifier
    >>> neigh = RadiusNeighborsClassifier(radius=1.0)
    >>> neigh.fit(X, y)
    RadiusNeighborsClassifier(algorithm='auto', leaf_size=30, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [0]

    See also
    --------
    KNeighborsClassifier
    RadiusNeighborsRegressor
    KNeighborsRegressor
    NearestNeighbors

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """
    def __init__(self, radius=5, algorithm='auto', leaf_size=30):
        self._init_params(radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)

    def predict(self, X):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X: array
            A 2-D array representing the test points.

        Returns
        -------
        labels: array
            List of class labels (one for each data sample).
        """
        X = atleast2d_or_csr(X)

        neigh_ind = self.radius_neighbors(X, return_distance=False)
        pred_labels = [self._y[ind] for ind in neigh_ind]
        mode = np.asarray([stats.mode(pl)[0] for pl in pred_labels],
                          dtype=np.int)
        return mode.flatten().astype(np.int)


@deprecated("deprecated in v0.9; will be removed in v0.11; "
            "use KNeighborsClassifier or RadiusNeighborsClassifier instead")
class NeighborsClassifier(NeighborsBase, KNeighborsMixin,
                          RadiusNeighborsMixin, SupervisedIntegerMixin,
                          ClassifierMixin):
    """Classifier implementing the nearest neighbors vote. (Deprecated)

    DEPRECATED IN VERSION 0.9; WILL BE REMOVED IN VERSION 0.11
    Please use :class:`KNeighborsClassifier` or
    :class:`RadiusNeighborsClassifier` instead.

    Samples participating in the vote are either the k-nearest neighbors
    (for some k)  or all neighbors within some fixed radius around the sample
    to classify.

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
        Note: fitting on sparse input will override the setting of this
        parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the construction and query, as well as the memory required to
        store the tree.  The optimal value depends on the nature of the
        problem.

    classification_type : {'knn_vote', 'radius_vote'}, optional
        Type of fit to use: 'knn_vote' specifies a k-NN classification.
        'radius_vote' specifies a r-NN classification.  Default is 'knn_vote'.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import NeighborsClassifier
    >>> neigh = NeighborsClassifier(n_neighbors=2)
    >>> neigh.fit(X, y)
    NeighborsClassifier(algorithm='auto', classification_type='knn_vote',
              leaf_size=30, n_neighbors=2, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [0]

    See also
    --------
    NearestNeighbors
    NeighborsRegressor

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    References
    ----------
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """
    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto', leaf_size=30,
                 classification_type='knn_vote'):
        if classification_type not in ('radius_vote', 'knn_vote'):
            raise ValueError("classification_type not recognized")
        self.classification_type = classification_type
        self._init_params(n_neighbors=n_neighbors,
                          radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)

    def predict(self, X):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X: array
            A 2-D array representing the test points.

        Returns
        -------
        labels: array
            List of class labels (one for each data sample).
        """
        X = atleast2d_or_csr(X)

        if self.classification_type == 'knn_vote':
            neigh_ind = self.kneighbors(X, return_distance=False)
            pred_labels = self._y[neigh_ind]
            mode, _ = stats.mode(pred_labels, axis=1)
            return mode.flatten().astype(np.int)
        else:
            neigh_ind = self.radius_neighbors(X, return_distance=False)
            pred_labels = [self._y[ind] for ind in neigh_ind]
            return np.asarray([stats.mode(pi) for pi in pred_labels],
                              dtype=np.int)
