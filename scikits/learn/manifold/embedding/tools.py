
"""
Tools for computation
"""

__all__ = ['create_graph', 'create_sym_graph', 'centered_normalized', 'dist2hd']

import numpy as np

from ...neighbors import Neighbors

def create_graph(samples, **kwargs):
    """
    Creates a list of list containing the nearest neighboors for each point in the dataset

    Parameters
    ----------
    samples : matrix
      The points to consider.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this function will be passed to its constructor.

    neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not given.

    Examples
    --------
    The following example creates a graph from samples and outputs the
    first item, that is a tuple representing the distance from that
    element to all other elements in sample:
    
    >>> from scikits.learn.manifold.compression.tools import create_graph
    >>> graph = create_graph(data)
    >>> print graph[0]
    """
    n = len(samples)
    labels, graph = np.zeros(n), [None]*n

    neigh = kwargs.get('neigh', None)
    if neigh is None:
        neigh = Neighbors(k=kwargs.get('neighbors', 9))
        neigh.fit(samples, labels)
        neigh = neigh.kneighbors
    else:
        neigh = neigh(**kwargs)
        neigh.fit(X)

    for i in range(0, len(samples)):
        graph[i] = [neighboor for neighboor in neigh(samples[i])]

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
