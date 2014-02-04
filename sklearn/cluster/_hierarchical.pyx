# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>

import numpy as np
cimport numpy as np
cimport cython

ctypedef np.float64_t DOUBLE
ctypedef np.npy_intp INTP
ctypedef np.int8_t INT8

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults

np.import_array()

from sklearn.utils.fast_dict cimport IntFloatDict

# C++
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map

DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.intp
ctypedef np.intp_t ITYPE_t

# Reimplementation for MSVC support
cdef inline double fmax(double a, double b):
    return max(a, b)

###############################################################################
# Utilities for computing the ward momentum

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


###############################################################################
# Utilities for cutting and exploring a hierarchical tree

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


###############################################################################
# merge strategies implemented on IntFloatDicts

# These are used in the hierarchical clustering code, to implement
# merging between two clusters, defined as a dict containing node number
# as keys and edge weights as values.


@cython.boundscheck(False)
@cython.wraparound(False)
def max_merge(IntFloatDict a, IntFloatDict b,
              np.ndarray[ITYPE_t, ndim=1] mask,
              ITYPE_t n_a, ITYPE_t n_b):
    """Merge two IntFloatDicts with the max strategy: when the same key is
    present in the two dicts, the max of the two values is used.

    Parameters
    ==========
    a, b : IntFloatDict object
        The IntFloatDicts to merge
    mask : ndarray array of dtype integer and of dimension 1
        a mask for keys to ignore: if not mask[key] the corresponding key
        is skipped in the output dictionnary
    n_a, n_b : float
        n_a and n_b are weights for a and b for the merge strategy.
        They are not used in the case of a max merge.

    Returns
    =======
    out : IntFloatDict object
        The IntFloatDict resulting from the merge
    """
    cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_it = a.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_end = a.my_map.end()
    cdef ITYPE_t key
    cdef DTYPE_t value
    # First copy a into out
    while a_it != a_end:
        key = deref(a_it).first
        if mask[key]:
            out_obj.my_map[key] = deref(a_it).second
        inc(a_it)

    # Then merge b into out
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_it = out_obj.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_end = out_obj.my_map.end()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_it = b.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_end = b.my_map.end()
    while b_it != b_end:
        key = deref(b_it).first
        value = deref(b_it).second
        if mask[key]:
            out_it = out_obj.my_map.find(key)
            if out_it == out_end:
                # Key not found
                out_obj.my_map[key] = value
            else:
                deref(out_it).second = fmax(deref(out_it).second, value)
        inc(b_it)
    return out_obj


@cython.boundscheck(False)
@cython.wraparound(False)
def average_merge(IntFloatDict a, IntFloatDict b,
              np.ndarray[ITYPE_t, ndim=1] mask,
              ITYPE_t n_a, ITYPE_t n_b):
    """Merge two IntFloatDicts with the average strategy: when the 
    same key is present in the two dicts, the weighted average of the two 
    values is used.

    Parameters
    ==========
    a, b : IntFloatDict object
        The IntFloatDicts to merge
    mask : ndarray array of dtype integer and of dimension 1
        a mask for keys to ignore: if not mask[key] the corresponding key
        is skipped in the output dictionnary
    n_a, n_b : float
        n_a and n_b are weights for a and b for the merge strategy.
        They are used for a weighted mean.

    Returns
    =======
    out : IntFloatDict object
        The IntFloatDict resulting from the merge
    """
    cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_it = a.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator a_end = a.my_map.end()
    cdef ITYPE_t key
    cdef DTYPE_t value
    cdef DTYPE_t n_out = <DTYPE_t> (n_a + n_b)
    # First copy a into out
    while a_it != a_end:
        key = deref(a_it).first
        if mask[key]:
            out_obj.my_map[key] = deref(a_it).second
        inc(a_it)

    # Then merge b into out
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_it = out_obj.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator out_end = out_obj.my_map.end()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_it = b.my_map.begin()
    cdef cpp_map[ITYPE_t, DTYPE_t].iterator b_end = b.my_map.end()
    while b_it != b_end:
        key = deref(b_it).first
        value = deref(b_it).second
        if mask[key]:
            out_it = out_obj.my_map.find(key)
            if out_it == out_end:
                # Key not found
                out_obj.my_map[key] = value
            else:
                deref(out_it).second = (n_a * deref(out_it).second
                                        + n_b * value) / n_out
        inc(b_it)
    return out_obj


###############################################################################
# An edge object for fast comparisons 

cdef class WeightedEdge:
    cdef public ITYPE_t a
    cdef public ITYPE_t b
    cdef public DTYPE_t weight
    
    def __init__(self, DTYPE_t weight, ITYPE_t a, ITYPE_t b):
        self.weight = weight
        self.a = a
        self.b = b

    @cython.nonecheck(False)
    def __richcmp__(self, WeightedEdge other, int op):
        """Cython-specific comparison method.

        op is the comparison code::
            <   0
            ==  2
            >   4
            <=  1
            !=  3
            >=  5
        """
        if op == 0:
            return self.weight < other.weight
        elif op == 1:
            return self.weight <= other.weight
        elif op == 2:
            return self.weight == other.weight
        elif op == 3:
            return self.weight != other.weight
        elif op == 4:
            return self.weight > other.weight
        elif op == 5:
            return self.weight >= other.weight
        
    def __repr__(self):
        return "%s(weight=%f, a=%i, b=%i)" % (self.__class__.__name__,
                                              self.weight,
                                              self.a, self.b)

