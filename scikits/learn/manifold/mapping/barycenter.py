# -*- coding: utf-8 -*-
"""
Barycenter Mapping

Use the barycenter in an original space to create a projection inside an embedded space
"""

import numpy as np

from ...base import BaseEstimator

from ..embedding.barycenters import barycenter
from ..embedding.tools import create_neighborer

class Barycenter(BaseEstimator):
    """
    Barycenter Mapping
    
    Parameters
    ----------
    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments` 
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor

    tol : float
      Tolerance for the Gram matrix construction for computing the barycenter
      weights
    Attributes
    ----------
    embedding_ : array_like
        Embedding of the learning data
    
    X_ : array_like
        Original data that is embedded

    """
    def __init__(self, n_neighbors = None, neigh = None,
        neigh_alternate_arguments = None, tol = 1e-3):
        self.n_neighbors = n_neighbors
        self.neigh = neigh
        self.neigh_alternate_arguments = neigh_alternate_arguments
        self.tol = tol

    def fit(self, embedding):
        """
        Train the mapping

        Parameters
        ----------
        embedding : 
            An embedding (for instance a trained Isomap instance)
        """
        self.__X = embedding.X_
        self.__Y = embedding.embedding_
        self.neigh = create_neighborer(self.__X, self.neigh, self.n_neighbors,
            self.neigh_alternate_arguments)
        self.neigh.fit(self.__X)
        return self

    def transform(self, X):
        """
        Parameters
        ----------
        X : arraylike
            The coordinates in an original space
        
        Returns
        -------
        Coordinates in the embedded space
        """
        X = np.atleast_2d(np.asanyarray(X))
        dist, X_neighbors = self.neigh.predict(X)
        return np.asanyarray([np.dot(barycenter(x, self.__X[neighbors], self.tol),
            self.__Y[neighbors]) for x, neighbors in zip(X, X_neighbors)])
