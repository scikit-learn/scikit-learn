"""Unsupervised nearest neighbors learner"""
import numpy as np

from .base import construct_docstring, \
    NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin, UnsupervisedMixin


class NearestNeighbors(NeighborsBase, KNeighborsMixin,
                       RadiusNeighborsMixin, UnsupervisedMixin):
    """Unsupervised learner for implementing neighbor searches.

    Parameters
    ----------
    @INCLUDE n_neighbors

    @INCLUDE radius

    @INCLUDE algorithm

    @INCLUDE leaf_size

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

    @INCLUDE notes
    """
    __doc__ = construct_docstring(__doc__)

    def __init__(self, n_neighbors=5, radius=1.0,
                 algorithm='auto', leaf_size=30):
        self._init_params(n_neighbors=n_neighbors,
                          radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size)
