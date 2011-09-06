"""Nearest Neighbor Classification"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD, (C) INRIA, University of Amsterdam

import numpy as np
from scipy import linalg

from .base import barycenter_weights, \
    NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin, SupervisedMixinFloat
from ..base import RegressorMixin
from ..utils import atleast2d_or_csr

class NeighborsRegressor(NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin,
                         SupervisedMixinFloat,
                         RegressorMixin):
    """Regression based on nearest neighbors.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.
    Samples used for the regression are either the k-nearest points, or
    all points within some fixed radius.  

    Different modes for estimating the result can be set via parameter
    mode. 'barycenter' will apply the weights that best reconstruct
    the point from its neighbors while 'mean' will apply constant
    weights to each point.

    Parameters
    ----------
    n_neighbors : int, optional
        Number of neighbors to use by default for k-NN prediction.
        Default is 5.

    radius : float, optional
        Radius to use by default for r-NN prediction. Default is 1.0.

    mode : {'mean', 'barycenter'}, optional
        Weights to apply to labels.

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
        nature of the problem (see discussion in NearestNeighbors).
        Defaults to 20.

    classification_type : {'knn_vote', 'radius_vote'}, optional
        Type of fit to use: 'knn_vote' specifies a k-NN classification.
        'radius_vote' specifies a r-NN classification.  Default is 'knn_vote'.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import NeighborsRegressor
    >>> neigh = NeighborsRegressor(n_neighbors=2)
    >>> neigh.fit(X, y)
    NeighborsRegressor(algorithm='auto', classification_type='knn_vote',
              leaf_size=20, mode='mean', n_neighbors=2, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    See also
    --------
    NearestNeighbors
    NeighborsClassifier

    Notes
    -----
    See Discussion in NearestNeighbors docstring regarding the optimal
    choice of algorithm and leaf_size.
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5, radius=1.0,
                 mode='mean', algorithm='auto',
                 leaf_size=20, classification_type='knn_vote'):
        if classification_type not in ('radius_vote', 'knn_vote'):
            raise ValueError("classification_type not recognized")
        self.classification_type = classification_type
        if mode not in ('mean', 'barycenter'):
            raise ValueError("mode must be one of 'mean' or 'barycenter'")
        self.mode = mode
        self._init_params(n_neighbors=n_neighbors,
                          radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)

    def predict(self, X):
        """Predict the target for the provided data

        Parameters
        ----------
        X : array
            A 2-D array representing the test data.

        Returns
        -------
        y: array
            List of target values (one for each data sample).
        """
        X = atleast2d_or_csr(X)

        if self.classification_type == 'knn_vote':
            neigh_ind = self.kneighbors(X, return_distance=False)

            # compute interpolation on y
            if self.mode == 'barycenter':
                W = barycenter_weights(X, self._fit_X[neigh_ind])
                return (W * self._y[neigh_ind]).sum(axis=1)

            elif self.mode == 'mean':
                return np.mean(self._y[neigh_ind], axis=1)

            else:
                raise ValueError(
                    'Unsupported mode, must be one of "barycenter" or '
                    '"mean" but got %s instead' % self.mode)
        else:
            # barycenter_weights() will have to be modified to allow
            # calculation of weights with a different `k` at each point.
            raise NotImplementedError("classification_type = 'radius_vote'")




class KNeighborsRegressor(NeighborsBase, KNeighborsMixin,
                          SupervisedMixinFloat,
                          RegressorMixin):
    """Regression based on k-nearest neighbors.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.

    Parameters
    ----------
    n_neighbors : int, optional
        Number of neighbors to use by default for k-NN prediction.
        Default is 5.

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
        nature of the problem (see discussion in NearestNeighbors).
        Defaults to 20.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsRegressor
    >>> neigh = KNeighborsRegressor(n_neighbors=2)
    >>> neigh.fit(X, y)
    KNeighborsRegressor(algorithm='auto', leaf_size=20, n_neighbors=2)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    See also
    --------
    RadiusNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier
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
        """Predict the target for the provided data

        Parameters
        ----------
        X : array
            A 2-D array representing the test data.

        Returns
        -------
        y: array
            List of target values (one for each data sample).
        """
        X = atleast2d_or_csr(X)

        neigh_ind = self.kneighbors(X, return_distance=False)

        return np.mean(self._y[neigh_ind], axis=1)


class RadiusNeighborsRegressor(NeighborsBase, RadiusNeighborsMixin,
                               SupervisedMixinFloat,
                               RegressorMixin):
    """Regression based on neighbors within a fixed radius.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.

    Parameters
    ----------
    radius : float, optional
        Fixed radius defining the neighborhood of each point.
        Default is 1.0.

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
        nature of the problem (see discussion in NearestNeighbors).
        Defaults to 20.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsRegressor
    >>> neigh = RadiusNeighborsRegressor(radius=1.0)
    >>> neigh.fit(X, y)
    RadiusNeighborsRegressor(algorithm='auto', leaf_size=20, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    See also
    --------
    KNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier
    NearestNeighbors

    Notes
    -----
    See Discussion in NearestNeighbors docstring regarding the optimal
    choice of algorithm and leaf_size.
    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, radius=1.0, algorithm='auto', leaf_size=20):
        self._init_params(radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)

    def predict(self, X):
        """Predict the target for the provided data

        Parameters
        ----------
        X : array
            A 2-D array representing the test data.

        Returns
        -------
        y: array
            List of target values (one for each data sample).
        """
        X = atleast2d_or_csr(X)

        neigh_ind = self.radius_neighbors(X, return_distance=False)

        return np.array([np.mean(self._y[ind])
                         for ind in neigh_ind])
