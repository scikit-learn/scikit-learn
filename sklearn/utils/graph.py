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

from .graph_shortest_path import graph_shortest_path


###############################################################################
# Path and connected component analysis.
# Code adapted from networkx

def single_source_shortest_path_length(graph, source, cutoff=None):
    """Return the shortest path length from source to all reachable nodes.

    Returns a dictionary of shortest path lengths keyed by target.

    Parameters
    ----------
    graph: sparse matrix or 2D array (preferably LIL matrix)
        Adjacency matrix of the graph
    source : node label
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
    >>> single_source_shortest_path_length(graph, 0)
    {0: 0, 1: 1, 2: 2, 3: 3}
    >>> single_source_shortest_path_length(np.ones((6, 6)), 2)
    {0: 1, 1: 1, 2: 0, 3: 1, 4: 1, 5: 1}
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


if hasattr(sparse, 'connected_components'):
    connected_components = sparse.connected_components
else:
    from .sparsetools import connected_components


###############################################################################
# Graph laplacian
def graph_laplacian(csgraph, normed=False, return_diag=False):
    """ Return the Laplacian matrix of a directed graph.

    For non-symmetric graphs the out-degree is used in the computation.

    Parameters
    ----------
    csgraph : array_like or sparse matrix, 2 dimensions
        compressed-sparse graph, with shape (N, N).
    normed : bool, optional
        If True, then compute normalized Laplacian.
    return_diag : bool, optional
        If True, then return diagonal as well as laplacian.

    Returns
    -------
    lap : ndarray
        The N x N laplacian matrix of graph.
    diag : ndarray
        The length-N diagonal of the laplacian matrix.
        diag is returned only if return_diag is True.

    Notes
    -----
    The Laplacian matrix of a graph is sometimes referred to as the
    "Kirchoff matrix" or the "admittance matrix", and is useful in many
    parts of spectral graph theory.  In particular, the eigen-decomposition
    of the laplacian matrix can give insight into many properties of the graph.

    For non-symmetric directed graphs, the laplacian is computed using the
    out-degree of each node.
    """
    if csgraph.ndim != 2 or csgraph.shape[0] != csgraph.shape[1]:
        raise ValueError('csgraph must be a square matrix or array')

    if normed and (np.issubdtype(csgraph.dtype, np.int)
                   or np.issubdtype(csgraph.dtype, np.uint)):
        csgraph = csgraph.astype(np.float)

    if sparse.isspmatrix(csgraph):
        return _laplacian_sparse(csgraph, normed=normed,
                                 return_diag=return_diag)
    else:
        return _laplacian_dense(csgraph, normed=normed,
                                return_diag=return_diag)


def _laplacian_sparse(graph, normed=False, return_diag=False):
    n_nodes = graph.shape[0]
    if not graph.format == 'coo':
        lap = (-graph).tocoo()
    else:
        lap = -graph.copy()
    diag_mask = (lap.row == lap.col)
    if not diag_mask.sum() == n_nodes:
        # The sparsity pattern of the matrix has holes on the diagonal,
        # we need to fix that
        diag_idx = lap.row[diag_mask]
        diagonal_holes = list(set(range(n_nodes)).difference(diag_idx))
        new_data = np.concatenate([lap.data, np.ones(len(diagonal_holes))])
        new_row = np.concatenate([lap.row, diagonal_holes])
        new_col = np.concatenate([lap.col, diagonal_holes])
        lap = sparse.coo_matrix((new_data, (new_row, new_col)),
                                shape=lap.shape)
        diag_mask = (lap.row == lap.col)

    lap.data[diag_mask] = 0
    w = -np.asarray(lap.sum(axis=1)).squeeze()
    if normed:
        w = np.sqrt(w)
        w_zeros = (w == 0)
        w[w_zeros] = 1
        lap.data /= w[lap.row]
        lap.data /= w[lap.col]
        lap.data[diag_mask] = (1 - w_zeros[lap.row[diag_mask]]).astype(
            lap.data.dtype)
    else:
        lap.data[diag_mask] = w[lap.row[diag_mask]]

    if return_diag:
        return lap, w
    return lap


def _laplacian_dense(graph, normed=False, return_diag=False):
    n_nodes = graph.shape[0]
    lap = -np.asarray(graph)  # minus sign leads to a copy

    # set diagonal to zero
    lap.flat[::n_nodes + 1] = 0
    w = -lap.sum(axis=0)
    if normed:
        w = np.sqrt(w)
        w_zeros = (w == 0)
        w[w_zeros] = 1
        lap /= w
        lap /= w[:, np.newaxis]
        lap.flat[::n_nodes + 1] = (1 - w_zeros).astype(lap.dtype)
    else:
        lap.flat[::n_nodes + 1] = w.astype(lap.dtype)

    if return_diag:
        return lap, w
    return lap
