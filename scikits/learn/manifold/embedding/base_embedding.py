# -*- coding: utf-8 -*-

from ...base import BaseEstimator

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

    mapping_kind : object
      The type of mapper to use. Can be:
          * None : no mapping built
          * "Barycenter" (default) : Barycenter mapping
          * a class object : a class that will be instantiated with the
              arguments of this function
          * an instance : an instance that will be fit() and then
              transform()ed

    Attributes
    ----------
    embedding_ : array_like
        Embedding of the learning data

    X_ : array_like
        Original data that is embedded

    See Also
    --------


    Notes
    -----

    """
    def __init__(self, n_coords, n_neighbors=None, ball_tree=None,
        mapping_kind="Barycenter"):
        self.n_coords = n_coords
        self.n_neighbors = n_neighbors if n_neighbors is not None else 9
        self.ball_tree = ball_tree
        self.mapping_kind = mapping_kind

    def transform(self, X):
        if self.mapping:
            return self.mapping.transform(X)
        else:
            raise RuntimeError("No mapping was built for this embedding")
