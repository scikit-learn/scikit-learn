
from .embedding.geodesic_mds import isomap

class Isomap(object):
    """
    Isomap embedding object

    Parameters
    ----------
    temp_file : string
      name of a file for caching the distance matrix

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this function will be passed to its constructor.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not given.

    Attributes
    ----------
    embedding_ : array_like
       Embedding of the learning data
      
    See Also
    --------

   
    Notes
    -----
    
    .. [1] Tenenbaum, J. B., de Silva, V. and Langford, J. C.,
           "A Global Geometric Framework for Nonlinear Dimensionality 
           Reduction,"
           Science, 290(5500), pp. 2319-2323, 2000
    
    Examples
    --------  
    
    """
    def __init__(self, **embedded_opts):
        self.__embedded_opts = embedded_opts

    def transform(self, X):
        """
        Parameters
        ----------
        X : array_like
        The learning dataset
        
        Returns
        -------
        Self
        """
        self.embedding_, reduced_parameter_set = isomap(X,
                                **self.__embedded_opts)
        return self
