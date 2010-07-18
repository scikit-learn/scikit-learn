"""
Barycenter Mapping

Use the barycenter in an original space to create a projection inside an embedded space
"""

class Barycenter(object):
  """
  Barycenter Mapping
  
  Parameters
  ----------
  
  
  """
  def __init__(self, **kwargs):
    pass

  def fit(self, X, Y):
    """
    Train the mapping
    
    Parameters
    ----------
    X : arraylike
      The coordinates in an original space
    Y : arraylike
      The coordinates in an embedded space
    """
    self.__X = X
    self.__Y = Y
    
  
  def predict(self, X):
    """
    """
    pass
    