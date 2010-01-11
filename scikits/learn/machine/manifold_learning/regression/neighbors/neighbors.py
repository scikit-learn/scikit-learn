#!/usr/bin/env python

# Matthieu Brucher
# Last Change : 2008-04-15 10:55

from numpy.ctypeslib import ndpointer, load_library
import math
import numpy
import sys
import ctypes

# Load the library
#if sys.platform == 'win32':
#  _neighbors = load_library('neighbors', "\\".join(__file__.split("\\")[:-1]) + "\\release")
#else:
_neighbors = load_library('_neighbors', __file__)

_neighbors.allocate_neighborer.restype = ctypes.c_void_p
_neighbors.allocate_neighborer.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_ulong, ctypes.c_ulong, ctypes.c_ulong]

_neighbors.delete_neighborer.restype = ctypes.c_int
_neighbors.delete_neighborer.argtypes = [ctypes.c_void_p]

CALLBACK = ctypes.CFUNCTYPE(ctypes.c_ulong, ctypes.c_double, ctypes.c_ulong)

_neighbors.find_kneighbors.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_ulong, ctypes.c_ulong, CALLBACK]
_neighbors.find_parzen.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_ulong, ctypes.c_double, CALLBACK]

class Neighbors:
  """
  Wrapper with ctypes around the neighbors tree
  """
  def __init__(self, points, neighbors = 0, window_size = 1., **kwargs):
    """
    Creates the tree, with a number of level depending on the log of the number of elements and the number of coordinates
    Parameters :
      - points is the matrix with the points populating the space
      - neighbors is the default number of neighbors
      - window_size is the default window size
    """
    self.neighbors = neighbors
    self.window_size = window_size
    self.points = numpy.ascontiguousarray(points) # needed for saving the state
    self._neigh = _neighbors.allocate_neighborer(self.points.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), self.points.shape[0], self.points.shape[1], int(math.log(self.points.shape[0]) / (self.points.shape[1] * math.log(2))))

  def __getinitargs__(self):
    """
    Returns the state of the neighboorhood
    """
    return (self.points, self.neighbors, self.window_size)

  def __setstate__(self, state):
    pass

  def __getstate__(self):
    return {}

  def __del__(self, close_func = _neighbors.delete_neighborer):
    """
    Deletes the cost function
    """
    if not (self._neigh == 0):
      close_func(self._neigh)
      self._neigh = 0

  def kneighbors(self, point, *args, **kwargs):
    """
    Finds the K-neighbors of a point

    Parameters :
      - point is a new point
      - neighbors is the number of neighbors to get (default is the
        value passed to the constructor)

    Outputs a list of tuples in which each tuple has two components,
    the first one indicates the length to point, whereas the second
    one is the index of that point in the population matrix.

    In the following example, we construnct a Neighbors class from an
    array representing our data set and ask who's the closest point to
    [1,1,1]

    >>> import numpy as np
    >>> samples = np.array([ \
      [0., 0., 0.], \
      [0., .5, 0.], \
      [1., 1., .5]])
    >>> n = Neighbors(samples)
    >>> print n.kneighbors(np.array([1., 1., 1.]), 1)
    [(0.5, 2L)]

    As you can see, it returns [(0.5, 2L)], which means that the
    element is at distance 0.5 and is the third element of samples
    (indexes start at 0)
    """
    if len(args) > 0:
      neighbors = args[0]
    elif 'neighbors' in kwargs:
      neighbors = kwargs['neighbors']
    else:
      neighbors = self.neighbors

    self.results = []
    _neighbors.find_kneighbors(self._neigh, point.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(point), neighbors, CALLBACK(self.callback))
    return self.results


  def parzen(self, point, *args, **kwargs):
    """
    Finds the neighbors of a point in a Parzen window
    Parameters :
      - point is a new point
      - window_size is the size of the window (default is the value passed to the constructor)
    """
    if len(args) > 0:
      window_size = args[0]
    elif 'window_size' in kwargs:
      window_size = kwargs['window_size']
    else:
      window_size = self.window_size

    self.results = []
    _neighbors.find_parzen(self._neigh, point.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(point), window_size, CALLBACK(self.callback))
    return self.results

  def callback(self, distance, indice):
    """
    Callback for the searchs, populates self.results
    Parameters :
      - distance is the distance of the point to another point
      - indice is the indice of the other point
    """
    self.results.append((distance, indice))
    return len(self.results)

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

if __name__ == "__main__":
  points = numpy.array([[0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [1., 1., 0.],
                [1., 0., 1.],
                [0., 1., 1.],
                [1., 1., 1.]], numpy.float64)

  n = Neighbors(points, 4, 1.2)
  print n.kneighbors(points[0])
  print n.parzen(points[0])

  kn = Kneighbors(points, 4, 1.2)
  print kn(points[0])
  print kn.parzen(points[0])

  pn = Parzen(points, 4, 1.2)
  print pn.kneighbors(points[0])
  print pn(points[0])
