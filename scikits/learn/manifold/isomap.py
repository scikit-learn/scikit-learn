
from .mapping import Barycenter

class Isomap(object):
  """
  Isomap reduction object

  Parameters
  ----------
  reduction_opts : dict
    Dictionary of reduction options
 
  projection : dict
    Projection technique that will be used ("barycenter" by default)

  projection_opts : dict
    Dictionary of projection options

  Attributes
  ----------
  embedding_ : array_like
    Embedding of the learning data
    
  See Also
  --------

 
  Notes
  -----
  
  .. [1] Tenenbaum, J. B., de Silva, V. and Langford, J. C.,
         "A Global Geometric Framework for Nonlinear Dimensionality Reduction,"
         Science, 290(5500), pp. 2319-2323, 2000
  
  Examples
  --------  
  
  """
  def __init__(self, reduction_opts, projection = None, projection_opts = None):
    self.__reduction_opts = reduction_opts
    self.__projection = projection if projection is not None else Barycenter
    self.__projection_opts = projection_opts if projection_opts is not None else {}

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
    from compression.geodesic_mds import isomap
    
    self.embedding_, reduced_parameter_set = isomap(X, **self.__reduction_opts)
    self.mapping = self.__projection(**self.__projection_opts)
	
    return self
    
  def predict(self, Y):
    """
    Parameters
    ----------
    Y : array_like
      The learning dataset
    
    Returns
    -------
    The embedding of the new dataset
    """
    pass
	