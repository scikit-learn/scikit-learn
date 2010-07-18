
from .mapping import Barycenter
from .embedding.geodesic_mds import isomap

class Isomap(object):
    """
    Isomap embedding object

    Parameters
    ----------
    embedded_opts : dict
       Dictionary of embedding options
   
    mapping : dict
       mapping technique that will be used ("barycenter" by default)

    mapping_opts : dict
       Dictionary of mapping options

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
    def __init__(self, embedded_opts, mapping = None, mapping_opts = None):
        self.__embedded_opts = embedded_opts
        self.__mapping = mapping if mapping is not None else Barycenter
        self.__mapping_opts = mapping_opts if mapping_opts is not None else {}

    def fit(self, X):
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

      self.mapping = self.__mapping(**self.__mapping_opts)
      self.mapping.fit(X, self.embedding_)
    
      return self
      
    def predict(self, X):
      """
      Parameters
      ----------
      X : array_like
        A new sample
      
      Returns
      -------
      The embedding of the new dataset
      """
      return self.mapping.predict(X)
