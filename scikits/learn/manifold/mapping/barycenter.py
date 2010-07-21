"""
Barycenter Mapping

Use the barycenter in an original space to create a projection inside an embedded space
"""

import numpy as np

from ...neighbors import Neighbors

from ..compression.barycenters import barycenter

class Barycenter(object):
    """
    Barycenter Mapping
    
    Parameters
    ----------
    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this function will be passed to its constructor.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not given.

    """
    def __init__(self, **kwargs):
        neigh = kwargs.get('neigh', None)
        if neigh is None:
            self.neigh = Neighbors(k=kwargs.get('n_neighbors', 9))
        else:
            self.neigh = neigh(**kwargs)
        self.kwargs = kwargs

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
        self.neigh.fit(self.__X)
      
    def predict(self, X):
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
        dist, X_neighbors = self.neigh.kneighbors(X)
        return np.asanyarray([np.dot(barycenter(x, self.__X[neighbors], **self.kwargs), self.__Y[neighbors]) for x, neighbors in zip(X, X_neighbors)])
