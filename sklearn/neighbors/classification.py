"""Nearest Neighbor Classification"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#
# License: BSD, (C) INRIA, University of Amsterdam

import numpy as np
from scipy import stats
from ..utils.extmath import weighted_mode

from .base import \
    construct_docstring, _check_weights, _get_weights, \
    NeighborsBase, KNeighborsMixin,\
    RadiusNeighborsMixin, SupervisedIntegerMixin
from ..base import ClassifierMixin
from ..utils import atleast2d_or_csr, deprecated


class KNeighborsClassifier(NeighborsBase, KNeighborsMixin,
                           SupervisedIntegerMixin, ClassifierMixin):
    """Classifier implementing the k-nearest neighbors vote.

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
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> neigh = KNeighborsClassifier(n_neighbors=2)
    >>> neigh.fit(X, y) # doctest: +ELLIPSIS
    KNeighborsClassifier(...)
    >>> print neigh.predict([[1.5]])
    [0]

    See also
    --------
    RadiusNeighborsClassifier
    KNeighborsRegressor
    RadiusNeighborsRegressor
    NearestNeighbors

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

    def __init__(self, n_neighbors=5,
                 weights='uniform',
                 algorithm='auto', leaf_size=30):
        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
        self.weights = _check_weights(weights)

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

        neigh_dist, neigh_ind = self.kneighbors(X)
        pred_labels = self._y[neigh_ind]

        weights = _get_weights(neigh_dist, self.weights)

        if weights is None:
            mode, _ = stats.mode(pred_labels, axis=1)
        else:
            mode, _ = weighted_mode(pred_labels, weights, axis=1)

        return mode.flatten().astype(np.int)


class RadiusNeighborsClassifier(NeighborsBase, RadiusNeighborsMixin,
                                SupervisedIntegerMixin, ClassifierMixin):
    """Classifier implementing a vote among neighbors within a given radius

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
    >>> from sklearn.neighbors import RadiusNeighborsClassifier
    >>> neigh = RadiusNeighborsClassifier(radius=1.0)
    >>> neigh.fit(X, y) # doctest: +ELLIPSIS
    RadiusNeighborsClassifier(...)
    >>> print neigh.predict([[1.5]])
    [0]

    See also
    --------
    KNeighborsClassifier
    RadiusNeighborsRegressor
    KNeighborsRegressor
    NearestNeighbors

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

    def __init__(self, radius=5, weights='uniform',
                 algorithm='auto', leaf_size=30):
        self._init_params(radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
        self.weights = _check_weights(weights)

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

        neigh_dist, neigh_ind = self.radius_neighbors(X)
        pred_labels = [self._y[ind] for ind in neigh_ind]

        weights = _get_weights(neigh_dist, self.weights)

        if weights is None:
            mode = np.asarray([stats.mode(pl)[0] for pl in pred_labels],
                              dtype=np.int)
        else:
            mode = np.asarray([weighted_mode(pl, w)[0]
                               for (pl, w) in zip(pred_labels, weights)],
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

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

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
