import numpy as np
cimport numpy as np
cimport cython
ctypedef np.float64_t DOUBLE
ctypedef np.int_t INT


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_ward_dist(np.ndarray[DOUBLE, ndim=1] m_1,\
                    np.ndarray[DOUBLE, ndim=2] m_2,\
                    np.ndarray[INT, ndim=1] coord_row,
                    np.ndarray[INT, ndim=1] coord_col,\
                    np.ndarray[DOUBLE, ndim=1] res):
    cdef int size_max = coord_row.shape[0]
    cdef int n_features = m_2.shape[1]
    cdef int i, j, row, col
    cdef DOUBLE pa, n
    for i in range(size_max):
        row = coord_row[i]
        col = coord_col[i]
        n = (m_1[row] * m_1[col]) / (m_1[row] + m_1[col])
        pa = 0.
        for j in range(n_features):
            pa += (m_2[row, j] / m_1[row] - m_2[col, j] / m_1[col])**2
        res[i] = pa * n
    return res


def _hc_get_descendent(int node, children, int n_leaves):
    """
    Function returning all the descendent leaves of a set of nodes in the tree.

    Parameters
    ----------
    node : int
        The node for which we want the descendents.

    children : list of pairs. Length of n_nodes
        List of the children of each nodes.
        This is not defined for leaves.

    n_leaves : int
        Number of leaves.

    Return
    ------
    descendent : list of int
    """
    ind = [node]
    if node < n_leaves:
        return ind
    descendent = []
    # It is actually faster to do the accounting of the number of
    # elements is the list ourselves: len is a lengthy operation on a
    # chained list
    cdef int n_indices = 1
    cdef int i
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
def _get_parent(int node, np.ndarray[INT, ndim=1] parents):
    cdef int parent
    parent = parents[node]
    while parent != node:
        node = parent
        parent = parents[node]
    return node
