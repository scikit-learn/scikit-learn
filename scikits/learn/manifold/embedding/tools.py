# -*- coding: utf-8 -*-

"""
Tools for computation
"""

__all__ = \
    ['create_graph', 'create_sym_graph', 'centered_normalized', 'dist2hd']

import numpy as np

from ...neighbors import Neighbors

def create_neighborer(samples, neigh = None, n_neighbors = None,
    neigh_alternate_arguments = None):
    """
    Computes the barycenters of samples given as parameters and returns them.
    
    Parameters
    ----------
    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor class . `neigh_alternate_arguments` 
      will be passed to this class constructor.

    neigh_alternate_arguments : dictionary
      Dictionary of arguments that will be passed to the `neigh` constructor

    Returns
    -------
    A neighbor instance
    """
    if neigh is None:
        neigh = Neighbors(k = n_neighbors if n_neighbors is not None else 9)
        neigh.fit(samples)
        neigh.predict = neigh.kneighbors
    else:
        neigh = neigh(**neigh_alternate_arguments)
        neigh.fit(samples)
    return neigh

def create_graph(samples, neigh = None, n_neighbors = None, neigh_alternate_arguments = None):
    """
    Creates a list of list containing the nearest neighboors for each point in
    the dataset

    Parameters
    ----------
    samples : matrix
      The points to consider.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this
      function will be passed to its constructor.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    Examples
    --------
    The following example creates a graph from samples and outputs the
    first item, that is a tuple representing the distance from that
    element to all other elements in sample:

    >>> import numpy
    >>> from scikits.learn.manifold.embedding.tools import create_graph
    >>> samples = numpy.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> graph = create_graph(samples, n_neighbors = 3)
    >>> print graph[0]
    [array([ 0. ,  0.5,  0.5]), array([0, 5, 4])]
    """
    n = len(samples)
    labels, graph = np.zeros(n), [None]*n

    neigh = create_neighborer(samples, neigh = neigh,
        n_neighbors = n_neighbors,
        neigh_alternate_arguments = neigh_alternate_arguments)

    for i in range(0, len(samples)):
        graph[i] = [neighboor for neighboor in neigh.predict(samples[i])]

    return graph

def create_sym_graph(samples, neigh = None, n_neighbors = None, neigh_alternate_arguments = None):
    """
    Creates a list of lists containing the nearest neighboors for each point in
    the dataset. The list of lists is symmetric

    Parameters
    ----------
    samples : matrix
      The points to consider.

    neigh : Neighbors
      A neighboorer (optional). By default, a K-Neighbor research is done.
      If provided, neigh must be a functor. All parameters passed to this
      function will be passed to its constructor.

    n_neighbors : int
      The number of K-neighboors to use (optional, default 9) if neigh is not
      given.

    Examples
    --------
    The following example creates a graph from samples and outputs the
    first item, that is a tuple representing the distance from that
    element to all other elements in sample:

    >>> import numpy
    >>> from scikits.learn.manifold.embedding.tools import create_sym_graph
    >>> samples = numpy.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> graph = create_sym_graph(samples, n_neighbors = 3)
    >>> print graph[0]
    [0, 1, 4, 5]
    """
    neigh = create_neighborer(samples, neigh = neigh,
        n_neighbors = n_neighbors,
        neigh_alternate_arguments = neigh_alternate_arguments)

    graph = [set() for i in range(len(samples))]

    for point in range(0, len(samples)):
        for vertex in neigh.predict(samples[point])[1]:
            graph[point].add(vertex)
            graph[vertex].add(point)

    return [list(el) for el in graph]

def centered_normalized(samples):
    """
    Returns a set of samples that are centered and of variance 1

    >>> import numpy
    >>> from  scikits.learn.manifold.embedding.tools import centered_normalized
    >>> samples = numpy.array((0., 0., 0., \
      1., 0., 0., \
      0., 1., 0., \
      1., 1., 0., \
      0., .5, 0., \
      .5, 0., 0., \
      1., 1., 0.5, \
      )).reshape((-1,3))
    >>> centered_normalized(samples)
    array([[-1.08012345, -1.08012345, -0.40824829],
           [ 1.08012345, -1.08012345, -0.40824829],
           [-1.08012345,  1.08012345, -0.40824829],
           [ 1.08012345,  1.08012345, -0.40824829],
           [-1.08012345,  0.        , -0.40824829],
           [ 0.        , -1.08012345, -0.40824829],
           [ 1.08012345,  1.08012345,  2.44948974]])
    """
    centered = samples - np.mean(samples, axis=0)
    centered /= np.std(centered, axis=0)
    return centered

def dist2hd(x,y):
    """
    Generates a distance matrix
    """
    d = np.zeros((x.shape[0],y.shape[0]),dtype=x.dtype)
    for i in xrange(x.shape[1]):
        diff2 = x[:,i,None] - y[:,i]
        diff2 **= 2
        d += diff2
    np.sqrt(d,d)
    return d
