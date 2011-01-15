# -*- coding: utf-8 -*-

from ..base import BaseEstimator


class BaseEmbedding(BaseEstimator):
    """
    Basic embedding object

    Parameters
    ----------
    n_coords : int
      The dimension of the embedding space

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    ball_tree : BallTree
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments`
      will be passed to this class constructor.

    Attributes
    ----------
    embedding_ : array_like
        Embedding of the learning data

    X_ : array_like
        Original data that is embedded

    Notes
    -----

    """
    def __init__(self, n_coords, n_neighbors=None, ball_tree=None):
        self.n_coords = n_coords
        self.n_neighbors = n_neighbors if n_neighbors is not None else 9
        self.ball_tree = ball_tree
