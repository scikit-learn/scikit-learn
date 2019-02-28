"""
Graph utilities and algorithms

Graphs are represented with their adjacency matrices, preferably using
sparse matrices.
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD 3 clause

from scipy import sparse

from .graph_shortest_path import graph_shortest_path  # noqa


###############################################################################
# Path and connected component analysis.
# Code adapted from networkx

def single_source_shortest_path_length(graph, source, cutoff=None):
    """Return the shortest path length from source to all reachable nodes.

    Returns a dictionary of shortest path lengths keyed by target.

    Parameters
    ----------
    graph : sparse matrix or 2D array (preferably LIL matrix)
        Adjacency matrix of the graph
    source : integer
       Starting node for path
    cutoff : integer, optional
        Depth to stop the search - only
        paths of length <= cutoff are returned.

    Examples
    --------
    >>> from sklearn.utils.graph import single_source_shortest_path_length
    >>> import numpy as np
    >>> graph = np.array([[ 0, 1, 0, 0],
    ...                   [ 1, 0, 1, 0],
    ...                   [ 0, 1, 0, 1],
    ...                   [ 0, 0, 1, 0]])
    >>> list(sorted(single_source_shortest_path_length(graph, 0).items()))
    [(0, 0), (1, 1), (2, 2), (3, 3)]
    >>> graph = np.ones((6, 6))
    >>> list(sorted(single_source_shortest_path_length(graph, 2).items()))
    [(0, 1), (1, 1), (2, 0), (3, 1), (4, 1), (5, 1)]
    """
    if sparse.isspmatrix(graph):
        graph = graph.tolil()
    else:
        graph = sparse.lil_matrix(graph)
    seen = {}                   # level (number of hops) when seen in BFS
    level = 0                   # the current level
    next_level = [source]       # dict of nodes to check at next level
    while next_level:
        this_level = next_level     # advance to next level
        next_level = set()          # and start a new list (fringe)
        for v in this_level:
            if v not in seen:
                seen[v] = level     # set the level of vertex v
                next_level.update(graph.rows[v])
        if cutoff is not None and cutoff <= level:
            break
        level += 1
    return seen  # return all path lengths as dictionary
