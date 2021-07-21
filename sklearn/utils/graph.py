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

from .deprecation import deprecated


###############################################################################
# Path and connected component analysis.
# Code adapted from networkx
def single_source_shortest_path_length(graph, source, *, cutoff=None):
    """Return the shortest path length from source to all reachable nodes.

    Returns a dictionary of shortest path lengths keyed by target.

    Parameters
    ----------
    graph : {sparse matrix, ndarray} of shape (n, n)
        Adjacency matrix of the graph. Sparse matrix of format LIL is
        preferred.

    source : int
       Starting node for path.

    cutoff : int, default=None
        Depth to stop the search - only paths of length <= cutoff are returned.

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
    seen = {}  # level (number of hops) when seen in BFS
    level = 0  # the current level
    next_level = [source]  # dict of nodes to check at next level
    while next_level:
        this_level = next_level  # advance to next level
        next_level = set()  # and start a new list (fringe)
        for v in this_level:
            if v not in seen:
                seen[v] = level  # set the level of vertex v
                next_level.update(graph.rows[v])
        if cutoff is not None and cutoff <= level:
            break
        level += 1
    return seen  # return all path lengths as dictionary


@deprecated(
    "`graph_shortest_path` is deprecated in 1.0 (renaming of 0.25) and will "
    "be removed in 1.2. Use `scipy.sparse.csgraph.shortest_path` instead."
)
def graph_shortest_path(dist_matrix, directed=True, method="auto"):
    """Shortest-path graph search on a positive directed or undirected graph.

    Parameters
    ----------
    dist_matrix : arraylike or sparse matrix, shape = (N,N)
        Array of positive distances.
        If vertex i is connected to vertex j, then dist_matrix[i,j] gives
        the distance between the vertices.
        If vertex i is not connected to vertex j, then dist_matrix[i,j] = 0

    directed : boolean
        if True, then find the shortest path on a directed graph: only
        progress from a point to its neighbors, not the other way around.
        if False, then find the shortest path on an undirected graph: the
        algorithm can progress from a point to its neighbors and vice versa.

    method : string ['auto'|'FW'|'D']
        method to use.  Options are
        'auto' : attempt to choose the best method for the current problem
        'FW' : Floyd-Warshall algorithm.  O[N^3]
        'D' : Dijkstra's algorithm with Fibonacci stacks.  O[(k+log(N))N^2]

    Returns
    -------
    G : np.ndarray, float, shape = [N,N]
        G[i,j] gives the shortest distance from point i to point j
        along the graph.

    Notes
    -----
    As currently implemented, Dijkstra's algorithm does not work for
    graphs with direction-dependent distances when directed == False.
    i.e., if dist_matrix[i,j] and dist_matrix[j,i] are not equal and
    both are nonzero, method='D' will not necessarily yield the correct
    result.
    Also, these routines have not been tested for graphs with negative
    distances.  Negative distances can lead to infinite cycles that must
    be handled by specialized algorithms.
    """
    return sparse.csgraph.shortest_path(dist_matrix, method=method, directed=directed)
