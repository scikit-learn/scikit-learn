#!/usr/bin/env python

# Matthieu Brucher
# Last Change : 2007-10-25 14:01

from numpy.ctypeslib import ndpointer, load_library
import math
import numpy
import sys
import ctypes

# Load the library
#if sys.platform == 'win32':
#  _neighboors = load_library('neighboors', "\\".join(__file__.split("\\")[:-1]) + "\\release")
#else:
_neighboors = load_library('_neighboors', __file__)

_neighboors.allocate_neighbooorer.restype = ctypes.c_void_p
_neighboors.allocate_neighbooorer.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_ulong, ctypes.c_ulong, ctypes.c_ulong]

_neighboors.delete_neighboorer.restype = ctypes.c_int
_neighboors.delete_neighboorer.argtypes = [ctypes.c_void_p]

CALLBACK = ctypes.CFUNCTYPE(ctypes.c_ulong, ctypes.c_double, ctypes.c_ulong)

_neighboors.find_kneighboors.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_ulong, ctypes.c_ulong, CALLBACK]
_neighboors.find_parzen.argtypes = [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_ulong, ctypes.c_double, CALLBACK]

class Neighboors:
  """
  Wrapper with ctypes around the neighboors tree
  """
  def __init__(self, points, neighboors = 0, window_size = 1., **kwargs):
    """
    Creates the tree, with a number of level depending on the log of the number of elements and the number of coordinates
    Parameters :
      - points is the matrix with the points populating the space
      - neighboors is the default number of neighboors
      - window_size is the default window size
    """
    self.neighboors = neighboors
    self.window_size = window_size
    self.points = numpy.ascontiguousarray(points) # needed for saving the state
    self._neigh = _neighboors.allocate_neighbooorer(self.points.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), self.points.shape[0], self.points.shape[1], int(math.log(self.points.shape[0]) / (self.points.shape[1] * math.log(2))))

  def __getinitargs__(self):
    """
    Returns the state of the neighboorhood
    """
    return (self.points, self.neighboors, self.window_size)

  def __setstate__(self, state):
    pass

  def __getstate__(self):
    return {}

  def __del__(self, close_func = _neighboors.delete_neighboorer):
    """
    Deletes the cost function
    """
    if not (self._neigh == 0):
      close_func(self._neigh)
      self._neigh = 0

  def kneighboors(self, point, *args, **kwargs):
    """
    Finds the K-Neighboors of a point
    Parameters :
      - point is a new point
      - neighboors is the number of neighboors to get (default is the value passed to the constructor)
    """
    if len(args) > 0:
      neighboors = args[0]
    elif 'neighboors' in kwargs:
      neighboors = kwargs['neighboors']
    else:
      neighboors = self.neighboors

    self.results = []
    _neighboors.find_kneighboors(self._neigh, point.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(point), neighboors, CALLBACK(self.callback))
    return self.results


  def parzen(self, point, *args, **kwargs):
    """
    Finds the neighboors of a point in a Parzen window
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
    _neighboors.find_parzen(self._neigh, point.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), len(point), window_size, CALLBACK(self.callback))
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

class KNeighboors(Neighboors):
  """
  Wrapper for K-Neighboors only
  """
  __call__ = Neighboors.kneighboors


class Parzen(Neighboors):
  """
  Wrapper for Parzen Window only
  """
  __call__ = Neighboors.parzen

if __name__ == "__main__":
  points = numpy.array([[0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., 0., 1.],
                [1., 1., 0.],
                [1., 0., 1.],
                [0., 1., 1.],
                [1., 1., 1.]], numpy.float64)

  n = Neighboors(points, 4, 1.2)
  print n.kneighboors(points[0])
  print n.parzen(points[0])

  kn = KNeighboors(points, 4, 1.2)
  print kn(points[0])
  print kn.parzen(points[0])

  pn = Parzen(points, 4, 1.2)
  print pn.kneighboors(points[0])
  print pn(points[0])
