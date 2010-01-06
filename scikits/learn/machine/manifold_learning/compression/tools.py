
"""
Tools for computation
"""

# Matthieu Brucher
# Last Change : 2008-02-28 09:33

__all__ = ['create_graph', 'create_sym_graph', 'centered_normalized', 'dist2hd']

import numpy

def create_graph(samples, **kwargs):
  """
  Creates a list of list containing the nearest neighboors for each point in the dataset
  Parameters :
  - samples is the points to consider
  - neigh is a neighboorer (optional)
  - neighboors is the number of K-neighboors to use (optional, default 9) if neigh is not given
  """
  from ..regression.neighbors import Neighbors
  if 'neigh' in kwargs:
    neighboorer = kwargs['neigh'](samples, **kwargs)
  else:
    neighboorer = Neighbors(samples, kwargs.get('neighbors', 9), 1.2)

  graph = [None] * len(samples)

  for i in range(0, len(samples)):
    graph[i] = [neighboor for neighboor in neighboorer.kneighbors(samples[i])]

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
  centered = samples - numpy.mean(samples, axis=0)
  centered /= numpy.std(centered, axis=0)
  return centered

def dist2hd(x,y):
   """
   Generate a 'coordinate' of the solution at a time
   """
   d = numpy.zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
   for i in xrange(x.shape[1]):
       diff2 = x[:,i,None] - y[:,i]
       diff2 **= 2
       d += diff2
   numpy.sqrt(d,d)
   return d
