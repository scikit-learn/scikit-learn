# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE
ctypedef np.npy_intp INTP
ctypedef np.int8_t INT8

np.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_ward_dist(np.ndarray[DOUBLE, ndim=1, mode='c'] m_1,
                      np.ndarray[DOUBLE, ndim=2, mode='c'] m_2,
                      np.ndarray[INTP, ndim=1, mode='c'] coord_row,
                      np.ndarray[INTP, ndim=1, mode='c'] coord_col,
                      np.ndarray[DOUBLE, ndim=1, mode='c'] res):
    cdef INTP size_max = coord_row.shape[0]
    cdef INTP n_features = m_2.shape[1]
    cdef INTP i, j, row, col
    cdef DOUBLE pa, n

    for i in range(size_max):
        row = coord_row[i]
        col = coord_col[i]
        n = (m_1[row] * m_1[col]) / (m_1[row] + m_1[col])
        pa = 0.
        for j in range(n_features):
            pa += (m_2[row, j] / m_1[row] - m_2[col, j] / m_1[col]) ** 2
        res[i] = pa * n
    return res


def _hc_get_descendent(INTP node, children, INTP n_leaves):
    """
    Function returning all the descendent leaves of a set of nodes in the tree.

    Parameters
    ----------
    node : integer
        The node for which we want the descendents.

    children : list of pairs, length n_nodes
        The children of each non-leaf node. Values less than `n_samples` refer
        to leaves of the tree. A greater value `i` indicates a node with
        children `children[i - n_samples]`.

    n_leaves : integer
        Number of leaves.

    Returns
    -------
    descendent : list of int
    """
    ind = [node]
    if node < n_leaves:
        return ind
    descendent = []

    # It is actually faster to do the accounting of the number of
    # elements is the list ourselves: len is a lengthy operation on a
    # chained list
    cdef INTP i, n_indices = 1

    while n_indices:
        i = ind.pop()
        if i < n_leaves:
            descendent.append(i)
            n_indices -= 1
        else:
            ind.extend(children[i - n_leaves])
            n_indices += 1
    return descendent


@cython.boundscheck(False)
@cython.wraparound(False)
def hc_get_heads(np.ndarray[INTP, ndim=1] parents, copy=True):
    """Returns the heads of the forest, as defined by parents.

    Parameters
    ----------
    parents : array of integers
        The parent structure defining the forest (ensemble of trees)
    copy : boolean
        If copy is False, the input 'parents' array is modified inplace

    Returns
    -------
    heads : array of integers of same shape as parents
        The indices in the 'parents' of the tree heads

    """
    cdef INTP parent, node0, node, size
    if copy:
        parents = np.copy(parents)
    size = parents.size

    # Start from the top of the tree and go down
    for node0 in range(size - 1, -1, -1):
        node = node0
        parent = parents[node]
        while parent != node:
            parents[node0] = parent
            node = parent
            parent = parents[node]
    return parents


@cython.boundscheck(False)
@cython.wraparound(False)
def _get_parents(nodes, heads, np.ndarray[INTP, ndim=1] parents,
                 np.ndarray[INT8, ndim=1, mode='c'] not_visited):
    """Returns the heads of the given nodes, as defined by parents.

    Modifies 'heads' and 'not_visited' in-place.

    Parameters
    ----------
    nodes : list of integers
        The nodes to start from
    heads : list of integers
        A list to hold the results (modified inplace)
    parents : array of integers
        The parent structure defining the tree
    not_visited
        The tree nodes to consider (modified inplace)

    """
    cdef INTP parent, node

    for node in nodes:
        parent = parents[node]
        while parent != node:
            node = parent
            parent = parents[node]
        if not_visited[node]:
            not_visited[node] = 0
            heads.append(node)
    return heads
