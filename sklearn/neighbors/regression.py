"""Nearest Neighbor Classification"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD, (C) INRIA, University of Amsterdam

import numpy as np
from scipy import linalg

from .base import \
    construct_docstring, _get_weights, _check_weights, \
    NeighborsBase, KNeighborsMixin, \
    RadiusNeighborsMixin, SupervisedFloatMixin
from ..base import RegressorMixin
from ..utils import atleast2d_or_csr, deprecated


class KNeighborsRegressor(NeighborsBase, KNeighborsMixin,
                          SupervisedFloatMixin,
                          RegressorMixin):
    """Regression based on k-nearest neighbors.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.

    Parameters
    ----------
    @INCLUDE n_neighbors

    @INCLUDE weights

    @INCLUDE algorithm

    @INCLUDE leaf_size

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsRegressor
    >>> neigh = KNeighborsRegressor(n_neighbors=2)
    >>> neigh.fit(X, y)
    KNeighborsRegressor(algorithm='auto', leaf_size=30, n_neighbors=2)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    See also
    --------
    NearestNeighbors
    RadiusNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

    def __init__(self, n_neighbors=5, weights='uniform',
                 algorithm='auto', leaf_size=30):
        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
        self.weights = _check_weights(weights)

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

        neigh_dist, neigh_ind = self.kneighbors(X)

        weights = _get_weights(neigh_dist, self.weights)

        if weights is None:
            return np.mean(self._y[neigh_ind], axis=1)
        else:
            num = np.sum(self._y[neigh_ind] * weights, axis=1)
            denom = np.sum(weights, axis=1)
            return num / denom


class RadiusNeighborsRegressor(NeighborsBase, RadiusNeighborsMixin,
                               SupervisedFloatMixin,
                               RegressorMixin):
    """Regression based on neighbors within a fixed radius.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.

    Parameters
    ----------
    @INCLUDE radius

    @INCLUDE weights

    @INCLUDE algorithm

    @INCLUDE leaf_size

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsRegressor
    >>> neigh = RadiusNeighborsRegressor(radius=1.0)
    >>> neigh.fit(X, y)
    RadiusNeighborsRegressor(algorithm='auto', leaf_size=30, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    See also
    --------
    NearestNeighbors
    KNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

    def __init__(self, radius=1.0, weights='uniform',
                 algorithm='auto', leaf_size=30):
        self._init_params(radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
        self.weights = _check_weights(weights)

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

        neigh_dist, neigh_ind = self.radius_neighbors(X)

        weights = _get_weights(neigh_dist)

        if weights is None:
            return np.array([np.mean(self._y[ind])
                             for ind in neigh_ind])
        else:
            return np.array([(np.sum(self._y[ind] * weights[i])
                              / np.sum(weights[i]))
                             for (i, ind) in enumerate(neigh_ind)])


@deprecated("deprecated in v0.9; will be removed in v0.11; "
            "use KNeighborsRegressor or RadiusNeighborsRegressor instead")
class NeighborsRegressor(NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin,
                         SupervisedFloatMixin,
                         RegressorMixin):
    """Regression based on nearest neighbors. (Deprecated)

    DEPRECATED IN VERSION 0.9; WILL BE REMOVED IN VERSION 0.11
    Please use :class:`KNeighborsRegressor` or
    :class:`RadiusNeighborsRegressor` instead.

    The target is predicted by local interpolation of the targets
    associated of the nearest neighbors in the training set.
    Samples used for the regression are either the k-nearest points, or
    all points within some fixed radius.

    Parameters
    ----------
    @INCLUDE n_neighbors

    @INCLUDE radius

    @INCLUDE algorithm

    @INCLUDE leaf_size

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
              leaf_size=30, n_neighbors=2, radius=1.0)
    >>> print neigh.predict([[1.5]])
    [ 0.5]

    See also
    --------
    NearestNeighbors
    KNeighborsRegressor
    RadiusNeighborsRegressor
    KNeighborsClassifier
    RadiusNeighborsClassifier

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto',
                 leaf_size=30, classification_type='knn_vote'):
        if classification_type not in ('radius_vote', 'knn_vote'):
            raise ValueError("classification_type not recognized")
        self.classification_type = classification_type
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
            return np.mean(self._y[neigh_ind], axis=1)

        else:
            neigh_ind = self.radius_neighbors(X, return_distance=False)

            # compute interpolation on y
            return np.array([np.mean(self._y[ind])
                             for ind in neigh_ind])
