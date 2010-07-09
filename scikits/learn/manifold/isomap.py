
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
    Dictionary of reduction options

  Attributes
  ----------
  embedding : array_like
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
    pass

  def fit(self, X):
    """
    Parameters
    ----------
    X : array_like
      The learning dataset
      
    Returns
    -------
    The embedding of the training dataset
    """
    pass
    
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
	