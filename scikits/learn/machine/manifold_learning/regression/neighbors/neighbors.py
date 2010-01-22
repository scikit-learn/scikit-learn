"""
k-Nearest Neighbor Algorithm
"""

from scipy.stats import mode
from scipy.spatial.ckdtree import cKDTree as KDTree
import numpy as np

class Neighbors:
  """
  Classifier implementing k-Nearest Neighbor Algorithm.

  Parameters
  ----------
  data : array-like, shape (n, k)
      The data points to be indexed. This array is not copied, and so
      modifying this data will result in bogus results.
  labels : array
      An array representing labels for the data (only arrays of
      integers are supported).
  k : int
      default number of neighbors.
  window_size : float
      the default window size.

  Examples
  --------
  >>> samples = [[0.,0.,1.], [1.,0.,0.], [2.,2.,2.], [2.,5.,4.]]
  >>> labels = [0,0,1,1]
  >>> neigh = Neighbors(samples, labels=labels)
  >>> print neigh.predict([[0,0,0]])
  [0]
  """

  def __init__(self, data, labels, k = 1, window_size = 1.):
    """
    Internally uses scipy.spatial.KDTree for most of its algorithms.
    """
    self.kdtree = KDTree(data, leafsize=20)
    self._k = k
    self.window_size = window_size
    self.points = np.ascontiguousarray(data) # needed for saving the state
    self.labels = np.asarray(labels)
    self.label_range = [self.labels.min(), self.labels.max()]

  def __getinitargs__(self):
    """
    Returns the state of the neighboorhood
    """
    return (self.points, self._k, self.window_size)

  def __setstate__(self, state):
    pass

  def __getstate__(self):
    return {}

  def kneighbors(self, data, k=None):
    """
    Finds the K-neighbors of a point.

    Parameters
    ----------
    point : array-like
        The new point.
    k : int
        Number of neighbors to get (default is the value
        passed to the constructor).

    Returns
    -------
    dist : array
        Array representing the lenghts to point.
    ind : array
        Array representing the indices of the nearest points in the
        population matrix.

    Examples
    --------
    In the following example, we construnct a Neighbors class from an
    array representing our data set and ask who's the closest point to
    [1,1,1]

    >>> import numpy as np
    >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
    >>> labels = [0, 0, 1]
    >>> neigh = Neighbors(samples, labels=labels)
    >>> print neigh.kneighbors([1., 1., 1.])
    (0.5, 2)

    As you can see, it returns [0.5], and [2], which means that the
    element is at distance 0.5 and is the third element of samples
    (indexes start at 0). You can also query for multiple points:

    >>> print neigh.kneighbors([[0., 1., 0.], [1., 0., 1.]])
    (array([ 0.5       ,  1.11803399]), array([1, 2]))

    """
    if k is None: k = self._k
    return self.kdtree.query(data, k=k)


  def parzen(self, point, window_size=None):
    """
    Finds the neighbors of a point in a Parzen window
    Parameters :
      - point is a new point
      - window_size is the size of the window (default is the value passed to the constructor)
    """
    if window_size is None: window_size = self.window_size
    return self.kdtree.query_ball_point(data, p=1.)

  def predict(self, data):
    """
    Predict the class labels for the provided data.

    Parameters
    ----------
    data: matrix
        An array representing the test point.

    Returns
    -------
    labels: array
        List of class labels (one for each data sample).

    Examples
    --------
    >>> import numpy as np
    >>> labels = [0,0,1]
    >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
    >>> neigh = Neighbors(samples, labels=labels)
    >>> print neigh.predict([.2, .1, .2])
    0
    >>> print neigh.predict([[0., -1., 0.], [3., 2., 0.]])
    [0 1]
    """
    dist, ind = self.kneighbors(data)
    labels = self.labels[ind]
    if self._k == 1: return labels
    # search most common values along axis 1 of labels
    # this is much faster than scipy.stats.mode
    return np.apply_along_axis(
      lambda x: np.bincount(x).argmax(),
      axis=1,
      arr=labels)

class Kneighbors(Neighbors):
  """
  Wrapper for K-neighbors only
  """
  __call__ = Neighbors.kneighbors


class Parzen(Neighbors):
  """
  Wrapper for Parzen Window only
  """
  __call__ = Neighbors.parzen
