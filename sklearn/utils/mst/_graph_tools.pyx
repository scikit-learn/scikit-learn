"""
Tools and utilities for working with compressed sparse graphs
"""

# Author: Jake Vanderplas  -- <vanderplas@astro.washington.edu>
# License: BSD, (C) 2012

import numpy as np
cimport numpy as np

from scipy.sparse import csr_matrix, isspmatrix,\
    isspmatrix_csr, isspmatrix_csc, isspmatrix_lil

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# EPS is the precision of DTYPE
cdef DTYPE_t DTYPE_EPS = 1E-15

# NULL_IDX is the index used in predecessor matrices to store a non-path
cdef ITYPE_t NULL_IDX = -9999

def csgraph_from_masked(graph):
    """
    csgraph_from_masked(graph)

    Construct a CSR-format graph from a masked array.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    graph : MaskedArray
        Input graph.  Shape should be (n_nodes, n_nodes).

    Returns
    -------
    csgraph : csr_matrix
        Compressed sparse representation of graph, 
    """
    # check that graph is a square matrix
    graph = np.ma.asarray(graph)

    if graph.ndim != 2:
        raise ValueError("graph should have two dimensions")
    N = graph.shape[0]
    if graph.shape[1] != N:
        raise ValueError("graph should be a square array")

    # construct the csr matrix using graph and mask
    if np.ma.is_masked(graph):
        data = graph.compressed()
        mask = ~graph.mask
    else:
        data = graph.data
        mask = np.ones(graph.shape, dtype='bool')

    data = np.asarray(data, dtype=DTYPE, order='c')

    idx_grid = np.empty((N, N), dtype=ITYPE)
    idx_grid[:] = np.arange(N, dtype=ITYPE)
    indices = np.asarray(idx_grid[mask], dtype=ITYPE, order='c')

    indptr = np.zeros(N + 1, dtype=ITYPE)
    indptr[1:] = mask.sum(1).cumsum()

    return csr_matrix((data, indices, indptr), (N, N))


def csgraph_masked_from_dense(graph,
                              null_value=0,
                              nan_null=True,
                              infinity_null=True,
                              copy=True):
    """
    csgraph_masked_from_dense(graph, null_value=0, nan_null=True,
                              infinity_null=True, copy=True)

    Construct a masked array graph representation from a dense matrix.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    graph : array_like
        Input graph.  Shape should be (n_nodes, n_nodes).
    null_value : float or None (optional)
        Value that denotes non-edges in the graph.  Default is zero.
    infinity_null : bool
        If True (default), then infinite entries (both positive and negative)
        are treated as null edges.
    nan_null : bool
        If True (default), then NaN entries are treated as non-edges

    Returns
    -------
    csgraph : MaskedArray
        masked array representation of graph
    """
    graph = np.array(graph, copy=copy)

    # check that graph is a square matrix
    if graph.ndim != 2:
        raise ValueError("graph should have two dimensions")
    N = graph.shape[0]
    if graph.shape[1] != N:
        raise ValueError("graph should be a square array")

    # check whether null_value is infinity or NaN
    if null_value is not None:
        null_value = DTYPE(null_value)        
        if np.isnan(null_value):
            nan_null = True
            null_value = None
        elif np.isinf(null_value):
            infinity_null = True
            null_value = None
    
    # flag all the null edges
    if null_value is None:
        mask = np.zeros(graph.shape, dtype='bool')
        graph = np.ma.masked_array(graph, mask, copy=False)
    else:
        graph = np.ma.masked_values(graph, null_value, copy=False)

    if infinity_null:
        graph.mask |= np.isinf(graph)

    if nan_null:
        graph.mask |= np.isnan(graph)

    return graph


def csgraph_from_dense(graph,
                       null_value=0,
                       nan_null=True,
                       infinity_null=True):
    """
    csgraph_from_dense(graph, null_value=0, nan_null=True, infinity_null=True)

    Construct a CSR-format sparse graph from a dense matrix.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    graph : array_like
        Input graph.  Shape should be (n_nodes, n_nodes).
    null_value : float or None (optional)
        Value that denotes non-edges in the graph.  Default is zero.
    infinity_null : bool
        If True (default), then infinite entries (both positive and negative)
        are treated as null edges.
    nan_null : bool
        If True (default), then NaN entries are treated as non-edges

    Returns
    -------
    csgraph : csr_matrix
        Compressed sparse representation of graph, 
    """
    return csgraph_from_masked(csgraph_masked_from_dense(graph,
                                                         null_value,
                                                         nan_null,
                                                         infinity_null))


def csgraph_to_dense(csgraph, null_value=0):
    """
    csgraph_to_dense(csgraph, null_value=0)

    Convert a sparse graph representation to a dense representation

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : csr_matrix, csc_matrix, or lil_matrix
        Sparse representation of a graph.
    null_value : float, optional
        The value used to indicate null edges in the dense representation.
        Default is 0.

    Returns
    -------
    graph : ndarray
        The dense representation of the sparse graph.

    Notes
    -----
    For normal sparse graph representations, calling csgraph_to_dense with
    null_value=0 produces an equivalent result to using dense format
    conversions in the main sparse package.  When the sparse representations
    have repeated values, however, the results will differ.  The tools in
    scipy.sparse will add repeating values to obtain a final value.  This
    function will select the minimum among repeating values to obtain a
    final value.  For example, here we'll create a two-node directed sparse
    graph with multiple edges from node 0 to node 1, of weights 2 and 3.
    This illustrates the difference in behavior:

    >>> from scipy.sparse import csr_matrix
    >>> data = np.array([2, 3])
    >>> indices = np.array([1, 1])
    >>> indptr = np.array([0, 2, 2])
    >>> M = csr_matrix((data, indices, indptr), shape=(2, 2))
    >>> M.toarray()
    array([[0, 5],
           [0, 0]])
    >>> csgraph_to_dense(M)
    array([[0, 2],
           [0, 0]])

    The reason for this difference is to allow a compressed sparse graph to
    represent multiple edges between any two nodes.  As most sparse graph
    algorithms are concerned with the single lowest-cost edge between any
    two nodes, the default scipy.sparse behavior of summming multiple weights
    does not make sense in this context.

    The other reason for using this routine is to allow for graphs with
    zero-weight edges.  Let's look at the example of a two-node directed
    graph, connected by an edge of weight zero:

    >>> from scipy.sparse import csr_matrix
    >>> data = np.array([0.0])
    >>> indices = np.array([1])
    >>> indptr = np.array([0, 2, 2])
    >>> M = csr_matrix((data, indices, indptr), shape=(2, 2))
    >>> M.toarray()
    array([[0, 0],
           [0, 0]])
    >>> csgraph_to_dense(M, np.inf)
    array([[ Inf,   0.],
           [ Inf,  Inf]])

    In the first case, the zero-weight edge gets lost in the dense
    representation.  In the second case, we can choose a different null value
    and see the true form of the graph.
    """
    # Allow only csr, lil and csc matrices: other formats when converted to csr
    # combine duplicated edges: we don't want this to happen in the background.
    if isspmatrix_csc(csgraph) or isspmatrix_lil(csgraph):
        csgraph = csgraph.tocsr()
    elif not isspmatrix_csr(csgraph):
        raise ValueError("csgraph must be lil, csr, or csc format")

    N = csgraph.shape[0]
    if csgraph.shape[1] != N:
        raise ValueError('csgraph should be a square matrix')

    # get attribute arrays
    data = np.asarray(csgraph.data, dtype=DTYPE, order='C')
    indices = np.asarray(csgraph.indices, dtype=ITYPE, order='C')
    indptr = np.asarray(csgraph.indptr, dtype=ITYPE, order='C')

    # create the output array
    graph = np.empty(csgraph.shape, dtype=DTYPE)
    graph.fill(np.inf)
    _populate_graph(data, indices, indptr, graph, null_value)
    return graph


def csgraph_to_masked(csgraph):
    """
    csgraph_to_masked(csgraph)

    Convert a sparse graph representation to a masked array representation

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : csr_matrix, csc_matrix, or lil_matrix
        Sparse representation of a graph.

    Returns
    -------
    graph : MaskedArray
        The masked dense representation of the sparse graph.
    """
    return np.ma.masked_invalid(csgraph_to_dense(csgraph, np.nan))


cdef void _populate_graph(np.ndarray[DTYPE_t, ndim=1, mode='c'] data,
                          np.ndarray[ITYPE_t, ndim=1, mode='c'] indices,
                          np.ndarray[ITYPE_t, ndim=1, mode='c'] indptr,
                          np.ndarray[DTYPE_t, ndim=2, mode='c'] graph,
                          DTYPE_t null_value):
    # data, indices, indptr are the csr attributes of the sparse input.
    # on input, graph should be filled with infinities, and should be
    # of size [N, N], which is also the size of the sparse matrix
    cdef unsigned int N = graph.shape[0]
    cdef np.ndarray null_flag = np.ones((N, N), dtype=bool, order='C')
    cdef np.npy_bool* null_ptr = <np.npy_bool*> null_flag.data
    cdef unsigned int row, col, i

    for row from 0 <= row < N:
        for i from indptr[row] <= i < indptr[row + 1]:
            col = indices[i]
            null_ptr[col] = 0
            # in case of multiple edges, we'll choose the smallest
            if data[i] < graph[row, col]:
                graph[row, col] = data[i]
        null_ptr += N

    graph[null_flag] = null_value


def reconstruct_path(csgraph, predecessors, directed=True):
    """
    reconstruct_path(csgraph, predecessors, directed=True)

    Construct a tree from a graph and a predecessor list.

    .. versionadded:: 0.11.0

    Parameters
    ----------
    csgraph : array_like or sparse matrix
        The N x N matrix representing the directed or undirected graph
        from which the predecessors are drawn.
    predecessors : array_like, one dimension
        The length-N array of indices of predecessors for the tree.  The
        index of the parent of node i is given by predecessors[i].
    directed : bool, optional
        If True (default), then operate on a directed graph: only move from
        point i to point j along paths csgraph[i, j].
        If False, then operate on an undirected graph: the algorithm can
        progress from point i to j along csgraph[i, j] or csgraph[j, i].

    Returns
    -------
    cstree : csr matrix
        The N x N directed compressed-sparse representation of the tree drawn
        from csgraph which is encoded by the predecessor list.
    """
    from _validation import validate_graph
    csgraph = validate_graph(csgraph, directed, dense_output=False)

    N = csgraph.shape[0]

    nnull = (predecessors < 0).sum()

    indices = np.argsort(predecessors)[nnull:].astype(ITYPE)
    pind = predecessors[indices]
    indptr = pind.searchsorted(np.arange(N + 1)).astype(ITYPE)

    if directed == True:
        data = csgraph[pind, indices]
    else:
        data1 = csgraph[pind, indices]
        data2 = csgraph[indices, pind]
        data1[data1 == 0] = np.inf
        data2[data2 == 0] = np.inf
        data = np.minimum(data1, data2)

    data = np.asarray(data).ravel()

    return csr_matrix((data, indices, indptr), shape=(N, N))


def construct_dist_matrix(graph,
                          predecessors,
                          directed=True,
                          null_value=np.inf):
    """
    construct_dist_matrix(graph, predecessors, directed=True, null_value=np.inf)

    Construct distance matrix from a predecessor matrix

    .. versionadded:: 0.11.0

    Parameters
    ----------
    graph : array_like or sparse
        The N x N matrix representation of a directed or undirected graph.
        If dense, then non-edges are indicated by zeros or infinities.
    predecessors : array_like
        The N x N matrix of predecessors of each node (see Notes below).
    directed : bool, optional
        If True (default), then operate on a directed graph: only move from
        point i to point j along paths csgraph[i, j].
        If False, then operate on an undirected graph: the algorithm can
        progress from point i to j along csgraph[i, j] or csgraph[j, i].
    null_value : bool, optional
        value to use for distances between unconnected nodes.  Default is
        np.inf

    Returns
    -------
    dist_matrix : ndarray
        The N x N matrix of distances between nodes along the path specified
        by the predecessor matrix.  If no path exists, the distance is zero.

    Notes
    -----
    The predecessor matrix is of the form returned by
    :func:`graph_shortest_path`.  Row i of the predecessor matrix contains
    information on the shortest paths from point i: each entry
    predecessors[i, j] gives the index of the previous node in the path from
    point i to point j.  If no path exists between point i and j, then
    predecessors[i, j] = -9999
    """
    from _validation import validate_graph
    graph = validate_graph(graph, directed, dtype=DTYPE,
                           csr_output=False,
                           copy_if_dense=not directed)
    predecessors = np.asarray(predecessors)

    if predecessors.shape != graph.shape:
        raise ValueError("graph and predecessors must have the same shape")

    dist_matrix = np.zeros(graph.shape, dtype=DTYPE)
    _construct_dist_matrix(graph, predecessors, dist_matrix,
                           directed, null_value)
    
    return dist_matrix


cdef void _construct_dist_matrix(np.ndarray[DTYPE_t, ndim=2] graph,
                                 np.ndarray[ITYPE_t, ndim=2] pred,
                                 np.ndarray[DTYPE_t, ndim=2] dist,
                                 int directed,
                                 DTYPE_t null_value):
    # All matrices should be size N x N
    # note that graph will be modified if directed == False
    # dist should be all zero on entry
    global NULL_IDX

    cdef int i, j, k1, k2, N, null_path
    N = graph.shape[0]

    #------------------------------------------
    # symmetrize matrix if necessary
    if not directed:
        graph[graph == 0] = np.inf
        for i from 0 <= i < N:
            for j from i + 1 <= j < N:
                if graph[j, i] <= graph[i, j]:
                    graph[i, j] = graph[j, i]
                else:
                    graph[j, i] = graph[i, j]
    #------------------------------------------

    for i from 0 <= i < N:
        for j from 0 <= j < N:
            null_path = True
            k2 = j
            while k2 != i:
                k1 = pred[i, k2]
                if k1 == NULL_IDX:
                    break
                dist[i, j] += graph[k1, k2]
                null_path = False
                k2 = k1
            if null_path and i != j:
                dist[i, j] = null_value
