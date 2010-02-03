
"""
Tools for computation
"""

# Matthieu Brucher
# Last Change : 2008-02-28 09:33

__all__ = ['create_graph', 'create_sym_graph', 'centered_normalized', 'dist2hd']

import numpy as np

def create_graph(samples, **kwargs):
  """
  Creates a list of list containing the nearest neighboors for each point in the dataset

  Parameters
  ----------
  samples : matrix
    The points to consider.

  neigh : Neighbors
    A neighboorer (optional).

  k : int
    The number of K-neighboors to use (optional, default 9) if neigh is not given.

  Examples
  --------
  The following example creates a graph from samples and outputs the
  first item, that is a tuple representing the distance from that
  element to all other elements in sample:

  >>> import numpy as np
  >>> samples = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [1., 1., 0.5]])
  >>> graph = create_graph(samples, k=3)
  >>> print graph[0]
  [array([ 0.,  1.,  1.]), array([0, 2, 1])]

  That is, it's at distance 0 from itself, at distance 1.0 from the
  second element and equally distance 1.0 from the third element.
  """
  from ..regression.neighbors import Neighbors

  n = len(samples)
  labels, graph = np.zeros(n), [None]*n

  neigh = kwargs.get('neigh', None)
  if neigh is None:
    neigh = Neighbors(samples, labels=labels, k=kwargs.get('k', 9))

  for i in range(0, len(samples)):
    graph[i] = [neighboor for neighboor in neigh.kneighbors(samples[i])]

  return graph

def create_sym_graph(samples, **kwargs):
  """
  Creates a list of list containing the nearest neighboors for each point in the dataset. The list of lists is symmetric
  Parameters :
    - samples is the points to consider
    - neigh is a neighboorer (optional)
    - neighboors is the number of K-neighboors to use (optional, default 9) if neigh is not given
  """
  import toolbox.neighboors
  if 'neigh' in kwargs:
    neighboorer = kwargs['neigh'](samples, **kwargs)
  else:
    neighboorer = toolbox.neighboors.distances.kneigh(samples, kwargs.get('neighboors', 9))

  graph = [set() for i in range(len(samples))]

  for point in range(0, len(samples)):
    for vertex in neighboorer[point][1:]:
      graph[point].add(vertex)
      graph[vertex].add(point)

  return [list(el) for el in graph]

def centered_normalized(samples):
  """
  Returns a set of samles that are centered and of variance 1
  """
  centered = samples - np.mean(samples, axis=0)
  centered /= np.std(centered, axis=0)
  return centered

def dist2hd(x,y):
   """
   Generate a 'coordinate' of the solution at a time
   """
   d = np.zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
   for i in xrange(x.shape[1]):
       diff2 = x[:,i,None] - y[:,i]
       diff2 **= 2
       d += diff2
   np.sqrt(d,d)
   return d
