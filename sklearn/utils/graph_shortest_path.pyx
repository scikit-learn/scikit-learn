"""
Routines for performing shortest-path graph searches

The main interface is in the function :func:`shortest_path`.  This
calls cython routines that compute the shortest path using
the Floyd-Warshall algorithm, Dijkstra's algorithm with Fibonacci Heaps,
the Bellman-Ford algorithm, or Johnson's Algorithm.

This module is a backdoor for scipy.sparse.csgraph.
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2011
import warnings

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix, isspmatrix_csr, isspmatrix_csc
from .sparsetools._graph_validation import validate_graph

cimport cython

from libc.stdlib cimport malloc, free

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# Fused type for int32 and int64
ctypedef fused int32_or_int64:
    np.int32_t
    np.int64_t

# EPS is the precision of DTYPE
cdef DTYPE_t DTYPE_EPS = 1E-15

# NULL_IDX is the index used in predecessor matrices to store a non-path
cdef ITYPE_t NULL_IDX = -9999


class NegativeCycleError(Exception):
    def __init__(self, message=''):
        Exception.__init__(self, message)


def shortest_path(csgraph, method='auto',
                  directed=True,
                  indices=None,
                  return_predecessors=False,
                  unweighted=False,
                  overwrite=False):
    """
    shortest_path(csgraph, method='auto', directed=True, return_predecessors=False,
                  unweighted=False, overwrite=False)

    Perform a shortest-path graph search on a positive directed or
    undirected graph.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    method : string ['auto'|'FW'|'D'], optional
        Algorithm to use for shortest paths.  Options are:

           'auto' -- (default) select the best among 'FW', 'D', 'BF', or 'J'
                     based on the input data.

           'FW'   -- Floyd-Warshall algorithm.  Computational cost is
                     approximately ``O[N^3]``.  The input csgraph will be
                     converted to a dense representation.

           'D'    -- Dijkstra's algorithm with Fibonacci heaps.  Computational
                     cost is approximately ``O[N(N*k + N*log(N))]``, where
		     ``k`` is the average number of connected edges per node.
		     The input csgraph will be converted to a csr
		     representation.

           'BF'   -- Bellman-Ford algorithm.  This algorithm can be used when
                     weights are negative.  If a negative cycle is encountered,
                     an error will be raised.  Computational cost is
                     approximately ``O[N(N^2 k)]``, where ``k`` is the average
                     number of connected edges per node. The input csgraph will
                     be converted to a csr representation.

           'J'    -- Johnson's algorithm.  Like the Bellman-Ford algorithm,
                     Johnson's algorithm is designed for use when the weights
                     are negative.  It combines the Bellman-Ford algorithm
                     with Dijkstra's algorithm for faster computation.

    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths for the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecesor matrix
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.
    overwrite : bool, optional
        If True, overwrite csgraph with the result.  This applies only if
        method == 'FW' and csgraph is a dense, c-ordered array with
        dtype=float64.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray
        Returned only if return_predecessors == True.
        The N x N matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    Notes
    -----
    As currently implemented, Dijkstra's algorithm and Johnson's algorithm
    do not work for graphs with direction-dependent distances when
    directed == False.  i.e., if csgraph[i,j] and csgraph[j,i] are non-equal
    edges, method='D' may yield an incorrect result.
    """
    # validate here to catch errors early but don't store the result;
    # we'll validate again later
    validate_graph(csgraph, directed, DTYPE,
                   copy_if_dense=(not overwrite),
                   copy_if_sparse=(not overwrite))

    if method == 'auto':
        # guess fastest method based on number of nodes and edges
        N = csgraph.shape[0]
        if isspmatrix(csgraph):
            Nk = csgraph.nnz
        else:
            Nk = np.sum(csgraph > 0)

        if indices is not None:
            # Indices were passed, so only indices*k
            # edges should be analyzed with fl.
            if isinstance(indices, int):
                n_indices = indices
            else:
                n_indices = len(indices)

            Nk = n_indices * Nk / N

        if Nk < N * N / 4:
            if ((isspmatrix(csgraph) and np.any(csgraph.data < 0))
                      or (not isspmatrix(csgraph) and np.any(csgraph < 0))):
                method = 'J'
            else:
                method = 'D'
        else:
            method = 'FW'

    if method == 'FW':
        return floyd_warshall(csgraph, directed,
                              indices=indices,
                              return_predecessors=return_predecessors,
                              unweighted=unweighted,
                              overwrite=overwrite)

    elif method == 'D':
        return dijkstra(csgraph, directed,
                        indices=indices,
                        return_predecessors=return_predecessors,
                        unweighted=unweighted)

    elif method == 'BF':
        return bellman_ford(csgraph, directed,
                            indices=indices,
                            return_predecessors=return_predecessors,
                            unweighted=unweighted)

    elif method == 'J':
        return johnson(csgraph, directed,
                       indices=indices,
                       return_predecessors=return_predecessors,
                       unweighted=unweighted)

    else:
        raise ValueError("unrecognized method '%s'" % method)


def floyd_warshall(csgraph, directed=True,
                   indices=None,
                   return_predecessors=False,
                   unweighted=False,
                   overwrite=False):
    """
    floyd_warshall(csgraph, directed=True, indices=None,
                   return_predecessors=False, unweighted=False,
                   overwrite=False)

    Compute the shortest path lengths using the Floyd-Warshall algorithm

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths for the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecesor matrix
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.
    overwrite : bool, optional
        If True, overwrite csgraph with the result.  This applies only if
        csgraph is a dense, c-ordered array with dtype=float64.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray
        Returned only if return_predecessors == True.
        The N x N matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph
    """
    dist_matrix = validate_graph(csgraph, directed, DTYPE,
                                 csr_output=False,
                                 copy_if_dense=not overwrite)

    N = csgraph.shape[0]

    #------------------------------
    # intitialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
        return_shape = indices.shape + (N,)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE, copy=True)
        return_shape = indices.shape + (N,)
        indices = np.atleast_1d(indices).reshape(-1)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")

    if unweighted:
        dist_matrix[~np.isinf(dist_matrix)] = 1

    if return_predecessors:
        predecessor_matrix = np.empty(dist_matrix.shape,
                                      dtype=ITYPE, order='C')
    else:
        predecessor_matrix = np.empty((0, 0), dtype=ITYPE)

    _floyd_warshall(dist_matrix,
                    predecessor_matrix,
                    int(directed))

    if np.any(dist_matrix.diagonal() < 0):
        raise NegativeCycleError("Negative cycle in nodes %s"
                                 % np.where(dist_matrix.diagonal() < 0)[0])

    dist_matrix = dist_matrix[indices]

    if return_predecessors:
        return dist_matrix, predecessor_matrix
    else:
        return dist_matrix


@cython.boundscheck(False)
cdef void _floyd_warshall(
               np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
               np.ndarray[ITYPE_t, ndim=2, mode='c'] predecessor_matrix,
               int directed=0):
    # dist_matrix : in/out
    #    on input, the graph
    #    on output, the matrix of shortest paths
    # dist_matrix should be a [N,N] matrix, such that dist_matrix[i, j]
    # is the distance from point i to point j.  Zero-distances imply that
    # the points are not connected.
    global NULL_IDX
    cdef int N = dist_matrix.shape[0]
    assert dist_matrix.shape[1] == N

    cdef unsigned int i, j, k

    cdef DTYPE_t infinity = np.inf
    cdef DTYPE_t d_ijk

    #----------------------------------------------------------------------
    #  Initialize distance matrix
    #   - set non-edges to infinity
    #   - set diagonal to zero
    #   - symmetrize matrix if non-directed graph is desired
    dist_matrix[dist_matrix == 0] = infinity
    dist_matrix.flat[::N + 1] = 0
    if not directed:
        for i from 0 <= i < N:
            for j from i + 1 <= j < N:
                if dist_matrix[j, i] <= dist_matrix[i, j]:
                    dist_matrix[i, j] = dist_matrix[j, i]
                else:
                    dist_matrix[j, i] = dist_matrix[i, j]

    #----------------------------------------------------------------------
    #  Initialize predecessor matrix
    #   - check matrix size
    #   - initialize diagonal and all non-edges to NULL
    #   - initialize all edges to the row index
    cdef int store_predecessors = False

    if predecessor_matrix.size > 0:
        store_predecessors = True
        assert predecessor_matrix.shape[0] == N
        assert predecessor_matrix.shape[1] == N
        predecessor_matrix.fill(NULL_IDX)
        i_edge = np.where(~np.isinf(dist_matrix))
        predecessor_matrix[i_edge] = i_edge[0]
        predecessor_matrix.flat[::N + 1] = NULL_IDX

    # Now perform the Floyd-Warshall algorithm.
    # In each loop, this finds the shortest path from point i
    #  to point j using intermediate nodes 0 ... k
    if store_predecessors:
        for k from 0 <= k < N:
            for i from 0 <= i < N:
                if dist_matrix[i, k] == infinity:
                    continue
                for j from 0 <= j < N:
                    d_ijk = dist_matrix[i, k] + dist_matrix[k, j]
                    if d_ijk < dist_matrix[i, j]:
                        dist_matrix[i, j] = d_ijk
                        predecessor_matrix[i, j] = predecessor_matrix[k, j]
    else:
        for k from 0 <= k < N:
            for i from 0 <= i < N:
                if dist_matrix[i, k] == infinity:
                    continue
                for j from 0 <= j < N:
                    d_ijk = dist_matrix[i, k] + dist_matrix[k, j]
                    if d_ijk < dist_matrix[i, j]:
                        dist_matrix[i, j] = d_ijk


def dijkstra(csgraph, directed=True, indices=None,
             return_predecessors=False,
             unweighted=False, limit=np.inf):
    """
    dijkstra(csgraph, directed=True, indices=None, return_predecessors=False,
             unweighted=False)

    Dijkstra algorithm using Fibonacci Heaps

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, 2 dimensions
        The N x N array of non-negative distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths for the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecesor matrix
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.
    limit : float, optional
        The maximum distance to calculate, must be >= 0. Using a smaller limit
        will decrease computation time by aborting calculations between pairs
        that are separated by a distance > limit. For such pairs, the distance
        will be equal to np.inf (i.e., not connected).
        .. versionadded:: 0.14.0

    Returns
    -------
    dist_matrix : ndarray
        The matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray
        Returned only if return_predecessors == True.
        The matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Notes
    -----
    As currently implemented, Dijkstra's algorithm does not work for
    graphs with direction-dependent distances when directed == False.
    i.e., if csgraph[i,j] and csgraph[j,i] are not equal and
    both are nonzero, setting directed=False will not yield the correct
    result.

    Also, this routine does not work for graphs with negative
    distances.  Negative distances can lead to infinite cycles that must
    be handled by specialized algorithms such as Bellman-Ford's algorithm
    or Johnson's algorithm.
    """
    global NULL_IDX

    #------------------------------
    # validate csgraph and convert to csr matrix
    csgraph = validate_graph(csgraph, directed, DTYPE,
                             dense_output=False)

    if np.any(csgraph.data < 0):
        warnings.warn("Graph has negative weights: dijkstra will give "
                      "inaccurate results if the graph contains negative "
                      "cycles. Consider johnson or bellman_ford.")

    N = csgraph.shape[0]

    #------------------------------
    # intitialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
        return_shape = indices.shape + (N,)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE, copy=True)
        return_shape = indices.shape + (N,)
        indices = np.atleast_1d(indices).reshape(-1)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")

    if not np.isscalar(limit):
        raise TypeError('limit must be numeric (float)')
    limit = float(limit)
    if limit < 0:
        raise ValueError('limit must be >= 0')

    #------------------------------
    # initialize dist_matrix for output
    dist_matrix = np.zeros((len(indices), N), dtype=DTYPE)
    dist_matrix.fill(np.inf)
    dist_matrix[np.arange(len(indices)), indices] = 0

    #------------------------------
    # initialize predecessors for output
    if return_predecessors:
        predecessor_matrix = np.empty((len(indices), N), dtype=ITYPE)
        predecessor_matrix.fill(NULL_IDX)
    else:
        predecessor_matrix = np.empty((0, N), dtype=ITYPE)

    if unweighted:
        csr_data = np.ones(csgraph.data.shape)
    else:
        csr_data = csgraph.data

    if directed:
        _dijkstra_directed(indices,
                           csr_data, csgraph.indices, csgraph.indptr,
                           dist_matrix, predecessor_matrix, limit)
    else:
        csgraphT = csgraph.T.tocsr()
        if unweighted:
            csrT_data = csr_data
        else:
            csrT_data = csgraphT.data
        _dijkstra_undirected(indices,
                             csr_data, csgraph.indices, csgraph.indptr,
                             csrT_data, csgraphT.indices, csgraphT.indptr,
                             dist_matrix, predecessor_matrix, limit)

    if return_predecessors:
        return (dist_matrix.reshape(return_shape),
                predecessor_matrix.reshape(return_shape))
    else:
        return dist_matrix.reshape(return_shape)


cdef _dijkstra_directed(
            np.ndarray[ITYPE_t, ndim=1, mode='c'] source_indices,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
            np.ndarray[ITYPE_t, ndim=2, mode='c'] pred,
            DTYPE_t limit):
    cdef unsigned int Nind = dist_matrix.shape[0]
    cdef unsigned int N = dist_matrix.shape[1]
    cdef unsigned int i, k, j_source, j_current
    cdef ITYPE_t j

    cdef DTYPE_t next_val

    cdef int return_pred = (pred.size > 0)

    cdef FibonacciHeap heap
    cdef FibonacciNode *v
    cdef FibonacciNode *current_node
    cdef FibonacciNode* nodes = <FibonacciNode*> malloc(N *
                                                        sizeof(FibonacciNode))

    for i from 0 <= i < Nind:
        j_source = source_indices[i]

        for k from 0 <= k < N:
            initialize_node(&nodes[k], k)

        dist_matrix[i, j_source] = 0
        heap.min_node = NULL
        insert_node(&heap, &nodes[j_source])

        while heap.min_node:
            v = remove_min(&heap)
            v.state = SCANNED

            for j from csr_indptr[v.index] <= j < csr_indptr[v.index + 1]:
                j_current = csr_indices[j]
                current_node = &nodes[j_current]
                if current_node.state != SCANNED:
                    next_val = v.val + csr_weights[j]
                    if next_val <= limit:
                        if current_node.state == NOT_IN_HEAP:
                            current_node.state = IN_HEAP
                            current_node.val = next_val
                            insert_node(&heap, current_node)
                            if return_pred:
                                pred[i, j_current] = v.index
                        elif current_node.val > next_val:
                            decrease_val(&heap, current_node,
                                         next_val)
                            if return_pred:
                                pred[i, j_current] = v.index

            #v has now been scanned: add the distance to the results
            dist_matrix[i, v.index] = v.val

    free(nodes)


cdef _dijkstra_undirected(
            np.ndarray[ITYPE_t, ndim=1, mode='c'] source_indices,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csrT_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csrT_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csrT_indptr,
            np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
            np.ndarray[ITYPE_t, ndim=2, mode='c'] pred,
            DTYPE_t limit):
    cdef unsigned int Nind = dist_matrix.shape[0]
    cdef unsigned int N = dist_matrix.shape[1]
    cdef unsigned int i, k, j_source, j_current
    cdef ITYPE_t j

    cdef DTYPE_t next_val

    cdef int return_pred = (pred.size > 0)

    cdef FibonacciHeap heap
    cdef FibonacciNode *v
    cdef FibonacciNode *current_node
    cdef FibonacciNode* nodes = <FibonacciNode*> malloc(N *
                                                        sizeof(FibonacciNode))

    for i from 0 <= i < Nind:
        j_source = source_indices[i]

        for k from 0 <= k < N:
            initialize_node(&nodes[k], k)

        dist_matrix[i, j_source] = 0
        heap.min_node = NULL
        insert_node(&heap, &nodes[j_source])

        while heap.min_node:
            v = remove_min(&heap)
            v.state = SCANNED

            for j from csr_indptr[v.index] <= j < csr_indptr[v.index + 1]:
                j_current = csr_indices[j]
                current_node = &nodes[j_current]
                if current_node.state != SCANNED:
                    next_val = v.val + csr_weights[j]
                    if next_val <= limit:
                        if current_node.state == NOT_IN_HEAP:
                            current_node.state = IN_HEAP
                            current_node.val = next_val
                            insert_node(&heap, current_node)
                            if return_pred:
                                pred[i, j_current] = v.index
                        elif current_node.val > next_val:
                            decrease_val(&heap, current_node,
                                         next_val)
                            if return_pred:
                                pred[i, j_current] = v.index

            for j from csrT_indptr[v.index] <= j < csrT_indptr[v.index + 1]:
                j_current = csrT_indices[j]
                current_node = &nodes[j_current]
                if current_node.state != SCANNED:
                    next_val = v.val + csrT_weights[j]
                    if next_val <= limit:
                        if current_node.state == NOT_IN_HEAP:
                            current_node.state = IN_HEAP
                            current_node.val = next_val
                            insert_node(&heap, current_node)
                            if return_pred:
                                pred[i, j_current] = v.index
                        elif current_node.val > next_val:
                            decrease_val(&heap, current_node, next_val)
                            if return_pred:
                                pred[i, j_current] = v.index

            #v has now been scanned: add the distance to the results
            dist_matrix[i, v.index] = v.val

    free(nodes)


def bellman_ford(csgraph, directed=True, indices=None,
                 return_predecessors=False,
                 unweighted=False):
    """
    bellman_ford(csgraph, directed=True, indices=None, return_predecessors=False,
                 unweighted=False)

    Compute the shortest path lengths using the Bellman-Ford algorithm.

    The Bellman-ford algorithm can robustly deal with graphs with negative
    weights.  If a negative cycle is detected, an error is raised.  For
    graphs without negative edge weights, dijkstra's algorithm may be faster.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths for the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecesor matrix
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray
        Returned only if return_predecessors == True.
        The N x N matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    Notes
    -----
    This routine is specially designed for graphs with negative edge weights.
    If all edge weights are positive, then Dijkstra's algorithm is a better
    choice.
    """
    global NULL_IDX

    #------------------------------
    # validate csgraph and convert to csr matrix
    csgraph = validate_graph(csgraph, directed, DTYPE,
                             dense_output=False)
    N = csgraph.shape[0]

    #------------------------------
    # intitialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")
    return_shape = indices.shape + (N,)
    indices = np.atleast_1d(indices).reshape(-1)

    #------------------------------
    # initialize dist_matrix for output
    dist_matrix = np.empty((len(indices), N), dtype=DTYPE)
    dist_matrix.fill(np.inf)
    dist_matrix[np.arange(len(indices)), indices] = 0

    #------------------------------
    # initialize predecessors for output
    if return_predecessors:
        predecessor_matrix = np.empty((len(indices), N), dtype=ITYPE)
        predecessor_matrix.fill(NULL_IDX)
    else:
        predecessor_matrix = np.empty((0, N), dtype=ITYPE)

    if unweighted:
        csr_data = np.ones(csgraph.data.shape)
    else:
        csr_data = csgraph.data

    if directed:
        ret = _bellman_ford_directed(indices,
                                     csr_data, csgraph.indices,
                                     csgraph.indptr,
                                     dist_matrix, predecessor_matrix)
    else:
        ret = _bellman_ford_undirected(indices,
                                       csr_data, csgraph.indices,
                                       csgraph.indptr,
                                       dist_matrix, predecessor_matrix)

    if ret >= 0:
        raise NegativeCycleError("Negative cycle detected on node %i" % ret)

    if return_predecessors:
        return (dist_matrix.reshape(return_shape),
                predecessor_matrix.reshape(return_shape))
    else:
        return dist_matrix.reshape(return_shape)


cdef int _bellman_ford_directed(
            np.ndarray[ITYPE_t, ndim=1, mode='c'] source_indices,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
            np.ndarray[ITYPE_t, ndim=2, mode='c'] pred):
    global DTYPE_EPS
    cdef unsigned int Nind = dist_matrix.shape[0]
    cdef unsigned int N = dist_matrix.shape[1]
    cdef unsigned int i, j, k, j_source, count

    cdef DTYPE_t d1, d2, w12

    cdef int return_pred = (pred.size > 0)

    for i from 0 <= i < Nind:
        j_source = source_indices[i]

        # relax all edges N-1 times
        for count from 0 <= count < N - 1:
            for j from 0 <= j < N:
                d1 = dist_matrix[i, j]
                for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
                    w12 = csr_weights[k]
                    d2 = dist_matrix[i, csr_indices[k]]
                    if d1 + w12 < d2:
                        dist_matrix[i, csr_indices[k]] = d1 + w12
                        if return_pred:
                            pred[i, csr_indices[k]] = j

        # check for negative-weight cycles
        for j from 0 <= j < N:
            d1 = dist_matrix[i, j]
            for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
                w12 = csr_weights[k]
                d2 = dist_matrix[i, csr_indices[k]]
                if d1 + w12 + DTYPE_EPS < d2:
                    return j_source

    return -1


cdef int _bellman_ford_undirected(
            np.ndarray[ITYPE_t, ndim=1, mode='c'] source_indices,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=2, mode='c'] dist_matrix,
            np.ndarray[ITYPE_t, ndim=2, mode='c'] pred):
    global DTYPE_EPS
    cdef unsigned int Nind = dist_matrix.shape[0]
    cdef unsigned int N = dist_matrix.shape[1]
    cdef unsigned int i, j, k, j_source, ind_k, count

    cdef DTYPE_t d1, d2, w12

    cdef int return_pred = (pred.size > 0)

    for i from 0 <= i < Nind:
        j_source = source_indices[i]

        # relax all edges N-1 times
        for count from 0 <= count < N - 1:
            for j from 0 <= j < N:
                d1 = dist_matrix[i, j]
                for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
                    w12 = csr_weights[k]
                    ind_k = csr_indices[k]
                    d2 = dist_matrix[i, ind_k]
                    if d1 + w12 < d2:
                        dist_matrix[i, ind_k] = d2 = d1 + w12
                        if return_pred:
                            pred[i, ind_k] = j
                    if d2 + w12 < d1:
                        dist_matrix[i, j] = d1 = d2 + w12
                        if return_pred:
                            pred[i, j] = ind_k

        # check for negative-weight cycles
        for j from 0 <= j < N:
            d1 = dist_matrix[i, j]
            for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
                w12 = csr_weights[k]
                d2 = dist_matrix[i, csr_indices[k]]
                if abs(d2 - d1) > w12 + DTYPE_EPS:
                    return j_source

    return -1


def johnson(csgraph, directed=True, indices=None,
            return_predecessors=False,
            unweighted=False):
    """
    johnson(csgraph, directed=True, indices=None, return_predecessors=False,
            unweighted=False)

    Compute the shortest path lengths using Johnson's algorithm.

    Johnson's algorithm combines the Bellman-Ford algorithm and Dijkstra's
    algorithm to quickly find shortest paths in a way that is robust to
    the presence of negative cycles.  If a negative cycle is detected,
    an error is raised.  For graphs without negative edge weights,
    dijkstra() may be faster.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array, matrix, or sparse matrix, 2 dimensions
        The N x N array of distances representing the input graph.
    directed : bool, optional
        If True (default), then find the shortest path on a directed graph:
        only move from point i to point j along paths csgraph[i, j].
        If False, then find the shortest path on an undirected graph: the
        algorithm can progress from point i to j along csgraph[i, j] or
        csgraph[j, i]
    indices : array_like or int, optional
        if specified, only compute the paths for the points at the given
        indices.
    return_predecessors : bool, optional
        If True, return the size (N, N) predecesor matrix
    unweighted : bool, optional
        If True, then find unweighted distances.  That is, rather than finding
        the path between each point such that the sum of weights is minimized,
        find the path such that the number of edges is minimized.

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between graph nodes. dist_matrix[i,j]
        gives the shortest distance from point i to point j along the graph.

    predecessors : ndarray
        Returned only if return_predecessors == True.
        The N x N matrix of predecessors, which can be used to reconstruct
        the shortest paths.  Row i of the predecessor matrix contains
        information on the shortest paths from point i: each entry
        predecessors[i, j] gives the index of the previous node in the
        path from point i to point j.  If no path exists between point
        i and j, then predecessors[i, j] = -9999

    Raises
    ------
    NegativeCycleError:
        if there are negative cycles in the graph

    Notes
    -----
    This routine is specially designed for graphs with negative edge weights.
    If all edge weights are positive, then Dijkstra's algorithm is a better
    choice.
    """
    #------------------------------
    # if unweighted, there are no negative weights: we just use dijkstra
    if unweighted:
        return dijkstra(csgraph, directed, indices,
                        return_predecessors, unweighted)

    #------------------------------
    # validate csgraph and convert to csr matrix
    csgraph = validate_graph(csgraph, directed, DTYPE,
                             dense_output=False)
    N = csgraph.shape[0]

    #------------------------------
    # initialize/validate indices
    if indices is None:
        indices = np.arange(N, dtype=ITYPE)
        return_shape = indices.shape + (N,)
    else:
        indices = np.array(indices, order='C', dtype=ITYPE)
        return_shape = indices.shape + (N,)
        indices = np.atleast_1d(indices).reshape(-1)
        indices[indices < 0] += N
        if np.any(indices < 0) or np.any(indices >= N):
            raise ValueError("indices out of range 0...N")

    #------------------------------
    # initialize dist_matrix for output
    dist_matrix = np.empty((len(indices), N), dtype=DTYPE)
    dist_matrix.fill(np.inf)
    dist_matrix[np.arange(len(indices)), indices] = 0

    #------------------------------
    # initialize predecessors for output
    if return_predecessors:
        predecessor_matrix = np.empty((len(indices), N), dtype=ITYPE)
        predecessor_matrix.fill(NULL_IDX)
    else:
        predecessor_matrix = np.empty((0, N), dtype=ITYPE)

    #------------------------------
    # initialize distance array
    dist_array = np.zeros(N, dtype=DTYPE)

    csr_data = csgraph.data.copy()

    #------------------------------
    # here we first add a single node to the graph, connected by a
    # directed edge of weight zero to each node, and perform bellman-ford
    if directed:
        ret = _johnson_directed(csr_data, csgraph.indices,
                                csgraph.indptr, dist_array)
    else:
        ret = _johnson_undirected(csr_data, csgraph.indices,
                                  csgraph.indptr, dist_array)

    if ret >= 0:
        raise NegativeCycleError("Negative cycle detected on node %i" % ret)

    #------------------------------
    # add the bellman-ford weights to the data
    _johnson_add_weights(csr_data, csgraph.indices,
                         csgraph.indptr, dist_array)

    if directed:
        _dijkstra_directed(indices,
                           csr_data, csgraph.indices, csgraph.indptr,
                           dist_matrix, predecessor_matrix, np.inf)
    else:
        csgraphT = csr_matrix((csr_data, csgraph.indices, csgraph.indptr),
                          csgraph.shape).T.tocsr()
        _johnson_add_weights(csgraphT.data, csgraphT.indices,
                             csgraphT.indptr, dist_array)
        _dijkstra_undirected(indices,
                             csr_data, csgraph.indices, csgraph.indptr,
                             csgraphT.data, csgraphT.indices, csgraphT.indptr,
                             dist_matrix, predecessor_matrix, np.inf)

    #------------------------------
    # correct the distance matrix for the bellman-ford weights
    dist_matrix += dist_array
    dist_matrix -= dist_array[:, None][indices]

    if return_predecessors:
        return (dist_matrix.reshape(return_shape),
                predecessor_matrix.reshape(return_shape))
    else:
        return dist_matrix.reshape(return_shape)


cdef void _johnson_add_weights(
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] dist_array):
    # let w(u, v) = w(u, v) + h(u) - h(v)
    cdef unsigned int j, k, N = dist_array.shape[0]

    for j from 0 <= j < N:
        for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
            csr_weights[k] += dist_array[j]
            csr_weights[k] -= dist_array[csr_indices[k]]


cdef int _johnson_directed(
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] dist_array):
    global DTYPE_EPS
    cdef unsigned int N = dist_array.shape[0]
    cdef unsigned int j, k, j_source, count

    cdef DTYPE_t d1, d2, w12

    # relax all edges (N+1) - 1 times
    for count from 0 <= count < N:
        for k from 0 <= k < N:
            if dist_array[k] < 0:
                dist_array[k] = 0

        for j from 0 <= j < N:
            d1 = dist_array[j]
            for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
                w12 = csr_weights[k]
                d2 = dist_array[csr_indices[k]]
                if d1 + w12 < d2:
                    dist_array[csr_indices[k]] = d1 + w12

    # check for negative-weight cycles
    for j from 0 <= j < N:
        d1 = dist_array[j]
        for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
            w12 = csr_weights[k]
            d2 = dist_array[csr_indices[k]]
            if d1 + w12 + DTYPE_EPS < d2:
                return j

    return -1


cdef int _johnson_undirected(
            np.ndarray[DTYPE_t, ndim=1, mode='c'] csr_weights,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indices,
            np.ndarray[ITYPE_t, ndim=1, mode='c'] csr_indptr,
            np.ndarray[DTYPE_t, ndim=1, mode='c'] dist_array):
    global DTYPE_EPS
    cdef unsigned int N = dist_array.shape[0]
    cdef unsigned int j, k, j_source, count

    cdef DTYPE_t d1, d2, w12

    # relax all edges (N+1) - 1 times
    for count from 0 <= count < N:
        for k from 0 <= k < N:
            if dist_array[k] < 0:
                dist_array[k] = 0

        for j from 0 <= j < N:
            d1 = dist_array[j]
            for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
                w12 = csr_weights[k]
                ind_k = csr_indices[k]
                d2 = dist_array[ind_k]
                if d1 + w12 < d2:
                    dist_array[ind_k] = d1 + w12
                if d2 + w12 < d1:
                    dist_array[j] = d1 = d2 + w12

    # check for negative-weight cycles
    for j from 0 <= j < N:
        d1 = dist_array[j]
        for k from csr_indptr[j] <= k < csr_indptr[j + 1]:
            w12 = csr_weights[k]
            d2 = dist_array[csr_indices[k]]
            if abs(d2 - d1) > w12 + DTYPE_EPS:
                return j

    return -1


######################################################################
# FibonacciNode structure
#  This structure and the operations on it are the nodes of the
#  Fibonacci heap.
#
cdef enum FibonacciState:
    SCANNED
    NOT_IN_HEAP
    IN_HEAP


cdef struct FibonacciNode:
    unsigned int index
    unsigned int rank
    FibonacciState state
    DTYPE_t val
    FibonacciNode* parent
    FibonacciNode* left_sibling
    FibonacciNode* right_sibling
    FibonacciNode* children


cdef void initialize_node(FibonacciNode* node,
                          unsigned int index,
                          DTYPE_t val=0):
    # Assumptions: - node is a valid pointer
    #              - node is not currently part of a heap
    node.index = index
    node.val = val
    node.rank = 0
    node.state = NOT_IN_HEAP

    node.parent = NULL
    node.left_sibling = NULL
    node.right_sibling = NULL
    node.children = NULL


cdef FibonacciNode* rightmost_sibling(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    cdef FibonacciNode* temp = node
    while(temp.right_sibling):
        temp = temp.right_sibling
    return temp


cdef FibonacciNode* leftmost_sibling(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    cdef FibonacciNode* temp = node
    while(temp.left_sibling):
        temp = temp.left_sibling
    return temp


cdef void add_child(FibonacciNode* node, FibonacciNode* new_child):
    # Assumptions: - node is a valid pointer
    #              - new_child is a valid pointer
    #              - new_child is not the sibling or child of another node
    new_child.parent = node

    if node.children:
        add_sibling(node.children, new_child)
    else:
        node.children = new_child
        new_child.right_sibling = NULL
        new_child.left_sibling = NULL
        node.rank = 1


cdef void add_sibling(FibonacciNode* node, FibonacciNode* new_sibling):
    # Assumptions: - node is a valid pointer
    #              - new_sibling is a valid pointer
    #              - new_sibling is not the child or sibling of another node
    cdef FibonacciNode* temp = rightmost_sibling(node)
    temp.right_sibling = new_sibling
    new_sibling.left_sibling = temp
    new_sibling.right_sibling = NULL
    new_sibling.parent = node.parent
    if new_sibling.parent:
        new_sibling.parent.rank += 1


cdef void remove(FibonacciNode* node):
    # Assumptions: - node is a valid pointer
    if node.parent:
        node.parent.rank -= 1
        if node.left_sibling:
            node.parent.children = node.left_sibling
        elif node.right_sibling:
            node.parent.children = node.right_sibling
        else:
            node.parent.children = NULL

    if node.left_sibling:
        node.left_sibling.right_sibling = node.right_sibling
    if node.right_sibling:
        node.right_sibling.left_sibling = node.left_sibling

    node.left_sibling = NULL
    node.right_sibling = NULL
    node.parent = NULL


######################################################################
# FibonacciHeap structure
#  This structure and operations on it use the FibonacciNode
#  routines to implement a Fibonacci heap

ctypedef FibonacciNode* pFibonacciNode


cdef struct FibonacciHeap:
    FibonacciNode* min_node
    pFibonacciNode[100] roots_by_rank  # maximum number of nodes is ~2^100.


cdef void insert_node(FibonacciHeap* heap,
                      FibonacciNode* node):
    # Assumptions: - heap is a valid pointer
    #              - node is a valid pointer
    #              - node is not the child or sibling of another node
    if heap.min_node:
        add_sibling(heap.min_node, node)
        if node.val < heap.min_node.val:
            heap.min_node = node
    else:
        heap.min_node = node


cdef void decrease_val(FibonacciHeap* heap,
                       FibonacciNode* node,
                       DTYPE_t newval):
    # Assumptions: - heap is a valid pointer
    #              - newval <= node.val
    #              - node is a valid pointer
    #              - node is not the child or sibling of another node
    #              - node is in the heap
    node.val = newval
    if node.parent and (node.parent.val >= newval):
        remove(node)
        insert_node(heap, node)
    elif heap.min_node.val > node.val:
        heap.min_node = node


cdef void link(FibonacciHeap* heap, FibonacciNode* node):
    # Assumptions: - heap is a valid pointer
    #              - node is a valid pointer
    #              - node is already within heap

    cdef FibonacciNode *linknode
    cdef FibonacciNode *parent
    cdef FibonacciNode *child

    if heap.roots_by_rank[node.rank] == NULL:
        heap.roots_by_rank[node.rank] = node
    else:
        linknode = heap.roots_by_rank[node.rank]
        heap.roots_by_rank[node.rank] = NULL

        if node.val < linknode.val or node == heap.min_node:
            remove(linknode)
            add_child(node, linknode)
            link(heap, node)
        else:
            remove(node)
            add_child(linknode, node)
            link(heap, linknode)


cdef FibonacciNode* remove_min(FibonacciHeap* heap):
    # Assumptions: - heap is a valid pointer
    #              - heap.min_node is a valid pointer
    cdef FibonacciNode *temp
    cdef FibonacciNode *temp_right
    cdef FibonacciNode *out
    cdef unsigned int i

    # make all min_node children into root nodes
    if heap.min_node.children:
        temp = leftmost_sibling(heap.min_node.children)
        temp_right = NULL

        while temp:
            temp_right = temp.right_sibling
            remove(temp)
            add_sibling(heap.min_node, temp)
            temp = temp_right

        heap.min_node.children = NULL

    # choose a root node other than min_node
    temp = leftmost_sibling(heap.min_node)
    if temp == heap.min_node:
        if heap.min_node.right_sibling:
            temp = heap.min_node.right_sibling
        else:
            out = heap.min_node
            heap.min_node = NULL
            return out

    # remove min_node, and point heap to the new min
    out = heap.min_node
    remove(heap.min_node)
    heap.min_node = temp

    # re-link the heap
    for i from 0 <= i < 100:
        heap.roots_by_rank[i] = NULL

    while temp:
        if temp.val < heap.min_node.val:
            heap.min_node = temp
        temp_right = temp.right_sibling
        link(heap, temp)
        temp = temp_right

    return out


######################################################################
# Debugging: Functions for printing the Fibonacci heap
#
#cdef void print_node(FibonacciNode* node, int level=0):
#    print '%s(%i,%i) %i' % (level*'   ', node.index, node.val, node.rank)
#    if node.children:
#        print_node(leftmost_sibling(node.children), level+1)
#    if node.right_sibling:
#        print_node(node.right_sibling, level)
#
#
#cdef void print_heap(FibonacciHeap* heap):
#    print "---------------------------------"
#    print "min node: (%i, %i)" % (heap.min_node.index, heap.min_node.val)
#    if heap.min_node:
#        print_node(leftmost_sibling(heap.min_node))
#    else:
#        print "[empty heap]"


@cython.boundscheck(False)
cdef void dijkstra_directed_one_row(unsigned int i_node,
                    unsigned int pos_in_graph,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors,
                    np.ndarray[DTYPE_t, ndim=1, mode='c'] distances,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                    np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes):
    """
    Calculate distances from a single point to all targets using a
    directed graph.

    Parameters
    ----------
    i_node : index of source point
    pos_in_graph : int
        index of source point in graph matrix.
        This might differ from i_node value if some
        shortest-path doesn't need to be computed.
    neighbors : array, shape = [N,]
        indices of neighbors for each point
    distances : array, shape = [N,]
        lengths of edges to each neighbor
    indptr : array, shape = (N+1,)
        the neighbors of point i are given by
        neighbors[indptr[i]:indptr[i+1]]
    graph : array, shape = (N,N)
        on return, graph[i_node] contains the path lengths from
        i_node to each target
    heap : the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int N = graph.shape[1]
    cdef unsigned int i
    cdef FibonacciNode *v
    cdef FibonacciNode *current_neighbor
    cdef DTYPE_t dist

    # initialize nodes
    for i from 0 <= i < N:
        initialize_node(&nodes[i], i)

    heap.min_node = NULL
    insert_node(heap, &nodes[i_node])

    while heap.min_node:
        v = remove_min(heap)
        v.state = 2  # 2 -> SCANNED

        for i from indptr[v.index] <= i < indptr[v.index + 1]:
            current_neighbor = &nodes[neighbors[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        #v has now been scanned: add the distance to the results
        graph[pos_in_graph, v.index] = v.val


@cython.boundscheck(False)
cdef void dijkstra_one_row(unsigned int i_node,
                    unsigned int pos_in_graph,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors1,
                    np.ndarray[DTYPE_t, ndim=1, mode='c'] distances1,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr1,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] neighbors2,
                    np.ndarray[DTYPE_t, ndim=1, mode='c'] distances2,
                    np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr2,
                    np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                    FibonacciHeap* heap,
                    FibonacciNode* nodes):
    """
    Calculate distances from a single point to all targets using an
    undirected graph.

    Parameters
    ----------
    i_node : index of source point
    pos_in_graph : int
        index of source point in graph matrix.
        This might differ from i_node value if some
        shortest-path doesn't need to be computed.
    neighbors[1,2] : array, shape = [N,]
        indices of neighbors for each point
    distances[1,2] : array, shape = [N,]
        lengths of edges to each neighbor
    indptr[1,2] : array, shape = (N+1,)
        the neighbors of point i are given by
        neighbors1[indptr1[i]:indptr1[i+1]] and
        neighbors2[indptr2[i]:indptr2[i+1]]
    graph : array, shape = (N,)
        on return, graph[i_node_in_graph] contains the path lengths from
        i_node to each target
    heap : the Fibonacci heap object to use
    nodes : the array of nodes to use
    """
    cdef unsigned int N = graph.shape[1]
    cdef unsigned int i
    cdef FibonacciNode *v
    cdef FibonacciNode *current_neighbor
    cdef DTYPE_t dist

    # re-initialize nodes
    # children, parent, left_sibling, right_sibling should already be NULL
    # rank should already be 0, index will already be set
    # we just need to re-set state and val
    for i from 0 <= i < N:
        nodes[i].state = 0  # 0 -> NOT_IN_HEAP
        nodes[i].val = 0

    insert_node(heap, &nodes[i_node])

    while heap.min_node:
        v = remove_min(heap)
        v.state = 2  # 2 -> SCANNED

        for i from indptr1[v.index] <= i < indptr1[v.index + 1]:
            current_neighbor = &nodes[neighbors1[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances1[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        for i from indptr2[v.index] <= i < indptr2[v.index + 1]:
            current_neighbor = &nodes[neighbors2[i]]
            if current_neighbor.state != 2:      # 2 -> SCANNED
                dist = distances2[i]
                if current_neighbor.state == 0:  # 0 -> NOT_IN_HEAP
                    current_neighbor.state = 1   # 1 -> IN_HEAP
                    current_neighbor.val = v.val + dist
                    insert_node(heap, current_neighbor)
                elif current_neighbor.val > v.val + dist:
                    decrease_val(heap, current_neighbor,
                                 v.val + dist)

        #v has now been scanned: add the distance to the results
        graph[pos_in_graph, v.index] = v.val
