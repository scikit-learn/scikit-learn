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
    NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin, SupervisedMixinInt
from ..base import ClassifierMixin
from ..utils import atleast2d_or_csr, deprecated


class KNeighborsClassifier(NeighborsBase, KNeighborsMixin,
                           SupervisedMixinInt, ClassifierMixin):
    """Classifier implementing the k-nearest neighbors vote.

    Parameters
    ----------
    n_neighbors : int, optional
        Number of neighbors to use by default for prediction. Default is 5.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree' will
        construct a BallTree, 'kd_tree' will use the KD-tree implementation
        in `scipy.spatial.ckdtree`, while 'brute' will perform brute-force
        search. 'auto' will guess the most appropriate based on current
        dataset.  See Discussion in `NearestNeighbors` for notes on choosing
        the optimal method. Fitting on sparse input will override the setting
        of this parameter.

    leaf_size : int, optional
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the nearest neighbors query.  The optimal value depends on the
        nature of the problem (see Discussion in NearestNeighbors).
        Defaults to 20.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> neigh = KNeighborsClassifier(n_neighbors=2)
    >>> neigh.fit(X, y)
    KNeighborsClassifier(algorithm='auto', leaf_size=20, n_neighbors=2)
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
    See Discussion in NearestNeighbors docstring regarding the optimal
    choice of algorithm and leaf_size.
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """
    def __init__(self, n_neighbors=5, algorithm='auto', leaf_size=20):
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
                                SupervisedMixinInt, ClassifierMixin):
    """Classifier implementing a vote among neighbors within a given radius

    Parameters
    ----------
    radius : float, optional
        Radius to use by default for prediction. Default is 1.0.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree' will
        construct a BallTree, 'kd_tree' will use the KD-tree implementation
        in `scipy.spatial.ckdtree`, while 'brute' will perform brute-force
        search. 'auto' will guess the most appropriate based on current
        dataset.  See Discussion in `NearestNeighbors` for notes on choosing
        the optimal method. Fitting on sparse input will override the setting
        of this parameter.

    leaf_size : int, optional
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the nearest neighbors query.  The optimal value depends on the
        nature of the problem (see Discussion in NearestNeighbors).
        Defaults to 20.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import RadiusNeighborsClassifier
    >>> neigh = RadiusNeighborsClassifier(radius=1.0)
    >>> neigh.fit(X, y)
    RadiusNeighborsClassifier(algorithm='auto', leaf_size=20, radius=1.0)
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
    See Discussion in NearestNeighbors docstring regarding the optimal
    choice of algorithm and leaf_size.
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """
    def __init__(self, radius=5, algorithm='auto', leaf_size=20):
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


#@deprecated("deprecated in v0.9; will be removed in v0.11; "
#            "use KNeighborsClassifier or RadiusNeighborsClassifier instead")
class NeighborsClassifier(NeighborsBase, KNeighborsMixin,
                          RadiusNeighborsMixin, SupervisedMixinInt,
                          ClassifierMixin):
    """Classifier implementing the nearest neighbors vote.

    Samples participating in the vote are either the k-nearest neighbors
    (for some k)  or all neighbors within some fixed radius around the sample
    to classify.

    Parameters
    ----------
    n_neighbors : int, optional
        Number of neighbors to use by default for k-NN prediction.
        Default is 5.

    radius : float, optional
        Radius to use by default for r-NN prediction. Default is 1.0.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors. 'ball_tree' will
        construct a BallTree, 'kd_tree' will use the KD-tree implementation
        in `scipy.spatial.ckdtree`, while 'brute' will perform brute-force
        search. 'auto' will guess the most appropriate based on current
        dataset.  See Discussion in `NearestNeighbors` for notes on choosing
        the optimal method. Fitting on sparse input will override the setting
        of this parameter.

    leaf_size : int, optional
        Leaf size passed to BallTree or cKDTree.  This can affect the speed
        of the nearest neighbors query.  The optimal value depends on the
        nature of the problem (see Discussion in NearestNeighbors).
        Defaults to 20.

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
              leaf_size=20, n_neighbors=2, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [0]

    See also
    --------
    NearestNeighbors
    NeighborsRegressor

    Notes
    -----
    See Discussion in NearestNeighbors docstring regarding the optimal
    choice of algorithm and leaf_size.
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """
    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto', leaf_size=20,
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
