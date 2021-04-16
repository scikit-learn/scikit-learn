"""
Graph utilities and algorithms

Graphs are represented with their adjacency matrices, preferably using
sparse matrices.
"""

# Authors: Aric Hagberg <hagberg@lanl.gov>
#          Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Jake Vanderplas <vanderplas@astro.washington.edu>
# License: BSD 3 clause

import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components

from .graph_shortest_path import graph_shortest_path  # noqa
from .validation import _deprecate_positional_args


###############################################################################
# Path and connected component analysis.
# Code adapted from networkx
@_deprecate_positional_args
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

def graph_connected_component(graph, node_id):
    """Find the largest graph connected components that contains one
    given node.

    Parameters
    ----------
    graph : array-like of shape (n_samples, n_samples)
        Adjacency matrix of the graph, non-zero weight means an edge
        between the nodes.

    node_id : int
        The index of the query node of the graph.

    Returns
    -------
    connected_components_matrix : array-like of shape (n_samples,)
        An array of bool value indicating the indexes of the nodes
        belonging to the largest connected components of the given query
        node.
    """
    n_node = graph.shape[0]
    if sparse.issparse(graph):
        # speed up row-wise access to boolean connection mask
        graph = graph.tocsr()
    connected_nodes = np.zeros(n_node, dtype=bool)
    nodes_to_explore = np.zeros(n_node, dtype=bool)
    nodes_to_explore[node_id] = True
    for _ in range(n_node):
        last_num_component = connected_nodes.sum()
        np.logical_or(connected_nodes, nodes_to_explore, out=connected_nodes)
        if last_num_component >= connected_nodes.sum():
            break
        indices = np.where(nodes_to_explore)[0]
        nodes_to_explore.fill(False)
        for i in indices:
            if sparse.issparse(graph):
                neighbors = graph[i].toarray().ravel()
            else:
                neighbors = graph[i]
            np.logical_or(nodes_to_explore, neighbors, out=nodes_to_explore)
    return connected_nodes

def graph_is_connected(graph):
    """ Return whether the graph is connected (True) or Not (False).

    Parameters
    ----------
    graph : {array-like, sparse matrix} of shape (n_samples, n_samples)
        Adjacency matrix of the graph, non-zero weight means an edge
        between the nodes.

    Returns
    -------
    is_connected : bool
        True means the graph is fully connected and False means not.
    """
    if sparse.isspmatrix(graph):
        # sparse graph, find all the connected components
        n_connected_components, _ = connected_components(graph)
        return n_connected_components == 1
    else:
        # dense graph, find all connected components start from node 0
        return graph_connected_component(graph, 0).sum() == graph.shape[0]
