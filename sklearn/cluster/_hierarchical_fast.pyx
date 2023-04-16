# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>

import numpy as np
cimport cython

from ..metrics._dist_metrics cimport DistanceMetric
from ..utils._fast_dict cimport IntFloatDict
from ..utils._typedefs cimport float64_t, intp_t, uint8_t

# C++
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp.map cimport map as cpp_map
from libc.math cimport fmax, INFINITY


###############################################################################
# Utilities for computing the ward momentum

def compute_ward_dist(
    const float64_t[::1] m_1,
    const float64_t[:, ::1] m_2,
    const intp_t[::1] coord_row,
    const intp_t[::1] coord_col,
    float64_t[::1] res
):
    cdef intp_t size_max = coord_row.shape[0]
    cdef intp_t n_features = m_2.shape[1]
    cdef intp_t i, j, row, col
    cdef float64_t pa, n

    for i in range(size_max):
        row = coord_row[i]
        col = coord_col[i]
        n = (m_1[row] * m_1[col]) / (m_1[row] + m_1[col])
        pa = 0.
        for j in range(n_features):
            pa += (m_2[row, j] / m_1[row] - m_2[col, j] / m_1[col]) ** 2
        res[i] = pa * n


###############################################################################
# Utilities for cutting and exploring a hierarchical tree

def _hc_get_descendent(intp_t node, children, intp_t n_leaves):
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
    cdef intp_t i, n_indices = 1

    while n_indices:
        i = ind.pop()
        if i < n_leaves:
            descendent.append(i)
            n_indices -= 1
        else:
            ind.extend(children[i - n_leaves])
            n_indices += 1
    return descendent


def hc_get_heads(intp_t[:] parents, copy=True):
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
    cdef intp_t parent, node0, node, size
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


def _get_parents(
    nodes,
    heads,
    const intp_t[:] parents,
    uint8_t[::1] not_visited
):
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
    cdef intp_t parent, node

    for node in nodes:
        parent = parents[node]
        while parent != node:
            node = parent
            parent = parents[node]
        if not_visited[node]:
            not_visited[node] = 0
            heads.append(node)


###############################################################################
# merge strategies implemented on IntFloatDicts

# These are used in the hierarchical clustering code, to implement
# merging between two clusters, defined as a dict containing node number
# as keys and edge weights as values.


def max_merge(
    IntFloatDict a,
    IntFloatDict b,
    const intp_t[:] mask,
    intp_t n_a,
    intp_t n_b
):
    """Merge two IntFloatDicts with the max strategy: when the same key is
    present in the two dicts, the max of the two values is used.

    Parameters
    ==========
    a, b : IntFloatDict object
        The IntFloatDicts to merge
    mask : ndarray array of dtype integer and of dimension 1
        a mask for keys to ignore: if not mask[key] the corresponding key
        is skipped in the output dictionary
    n_a, n_b : float
        n_a and n_b are weights for a and b for the merge strategy.
        They are not used in the case of a max merge.

    Returns
    =======
    out : IntFloatDict object
        The IntFloatDict resulting from the merge
    """
    cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
    cdef cpp_map[intp_t, float64_t].iterator a_it = a.my_map.begin()
    cdef cpp_map[intp_t, float64_t].iterator a_end = a.my_map.end()
    cdef intp_t key
    cdef float64_t value
    # First copy a into out
    while a_it != a_end:
        key = deref(a_it).first
        if mask[key]:
            out_obj.my_map[key] = deref(a_it).second
        inc(a_it)

    # Then merge b into out
    cdef cpp_map[intp_t, float64_t].iterator out_it = out_obj.my_map.begin()
    cdef cpp_map[intp_t, float64_t].iterator out_end = out_obj.my_map.end()
    cdef cpp_map[intp_t, float64_t].iterator b_it = b.my_map.begin()
    cdef cpp_map[intp_t, float64_t].iterator b_end = b.my_map.end()
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


def average_merge(
    IntFloatDict a,
    IntFloatDict b,
    const intp_t[:] mask,
    intp_t n_a,
    intp_t n_b
):
    """Merge two IntFloatDicts with the average strategy: when the
    same key is present in the two dicts, the weighted average of the two
    values is used.

    Parameters
    ==========
    a, b : IntFloatDict object
        The IntFloatDicts to merge
    mask : ndarray array of dtype integer and of dimension 1
        a mask for keys to ignore: if not mask[key] the corresponding key
        is skipped in the output dictionary
    n_a, n_b : float
        n_a and n_b are weights for a and b for the merge strategy.
        They are used for a weighted mean.

    Returns
    =======
    out : IntFloatDict object
        The IntFloatDict resulting from the merge
    """
    cdef IntFloatDict out_obj = IntFloatDict.__new__(IntFloatDict)
    cdef cpp_map[intp_t, float64_t].iterator a_it = a.my_map.begin()
    cdef cpp_map[intp_t, float64_t].iterator a_end = a.my_map.end()
    cdef intp_t key
    cdef float64_t value
    cdef float64_t n_out = <float64_t> (n_a + n_b)
    # First copy a into out
    while a_it != a_end:
        key = deref(a_it).first
        if mask[key]:
            out_obj.my_map[key] = deref(a_it).second
        inc(a_it)

    # Then merge b into out
    cdef cpp_map[intp_t, float64_t].iterator out_it = out_obj.my_map.begin()
    cdef cpp_map[intp_t, float64_t].iterator out_end = out_obj.my_map.end()
    cdef cpp_map[intp_t, float64_t].iterator b_it = b.my_map.begin()
    cdef cpp_map[intp_t, float64_t].iterator b_end = b.my_map.end()
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
    cdef public intp_t a
    cdef public intp_t b
    cdef public float64_t weight

    def __init__(self, float64_t weight, intp_t a, intp_t b):
        self.weight = weight
        self.a = a
        self.b = b

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


################################################################################
# Efficient labelling/conversion of MSTs to single linkage hierarchies

cdef class UnionFind(object):

    cdef intp_t next_label
    cdef intp_t[:] parent
    cdef intp_t[:] size

    def __init__(self, N):
        self.parent = np.full(2 * N - 1, -1., dtype=np.intp, order='C')
        self.next_label = N
        self.size = np.hstack((np.ones(N, dtype=np.intp),
                               np.zeros(N - 1, dtype=np.intp)))

    cdef void union(self, intp_t m, intp_t n):
        self.parent[m] = self.next_label
        self.parent[n] = self.next_label
        self.size[self.next_label] = self.size[m] + self.size[n]
        self.next_label += 1

        return

    @cython.wraparound(True)
    cdef intp_t fast_find(self, intp_t n):
        cdef intp_t p
        p = n
        # find the highest node in the linkage graph so far
        while self.parent[n] != -1:
            n = self.parent[n]
        # provide a shortcut up to the highest node
        while self.parent[p] != n:
            p, self.parent[p] = self.parent[p], n
        return n


def _single_linkage_label(const float64_t[:, :] L):
    """
    Convert an linkage array or MST to a tree by labelling clusters at merges.
    This is done by using a Union find structure to keep track of merges
    efficiently. This is the private version of the function that assumes that
    ``L`` has been properly validated. See ``single_linkage_label`` for the
    user facing version of this function.

    Parameters
    ----------
    L: array of shape (n_samples - 1, 3)
        The linkage array or MST where each row specifies two samples
        to be merged and a distance or weight at which the merge occurs. This
         array is assumed to be sorted by the distance/weight.

    Returns
    -------
    A tree in the format used by scipy.cluster.hierarchy.
    """

    cdef float64_t[:, ::1] result_arr

    cdef intp_t left, left_cluster, right, right_cluster, index
    cdef float64_t delta

    result_arr = np.zeros((L.shape[0], 4), dtype=np.float64)
    U = UnionFind(L.shape[0] + 1)

    for index in range(L.shape[0]):

        left = <intp_t> L[index, 0]
        right = <intp_t> L[index, 1]
        delta = L[index, 2]

        left_cluster = U.fast_find(left)
        right_cluster = U.fast_find(right)

        result_arr[index][0] = left_cluster
        result_arr[index][1] = right_cluster
        result_arr[index][2] = delta
        result_arr[index][3] = U.size[left_cluster] + U.size[right_cluster]

        U.union(left_cluster, right_cluster)

    return np.asarray(result_arr)


@cython.wraparound(True)
def single_linkage_label(L):
    """
    Convert an linkage array or MST to a tree by labelling clusters at merges.
    This is done by using a Union find structure to keep track of merges
    efficiently.

    Parameters
    ----------
    L: array of shape (n_samples - 1, 3)
        The linkage array or MST where each row specifies two samples
        to be merged and a distance or weight at which the merge occurs. This
         array is assumed to be sorted by the distance/weight.

    Returns
    -------
    A tree in the format used by scipy.cluster.hierarchy.
    """
    # Validate L
    if L[:, :2].min() < 0 or L[:, :2].max() >= 2 * L.shape[0] + 1:
        raise ValueError("Input MST array is not a validly formatted MST array")

    is_sorted = lambda x: np.all(x[:-1] <= x[1:])
    if not is_sorted(L[:, 2]):
        raise ValueError("Input MST array must be sorted by weight")

    return _single_linkage_label(L)


# Implements MST-LINKAGE-CORE from https://arxiv.org/abs/1109.2378
def mst_linkage_core(
        const float64_t [:, ::1] raw_data,
        DistanceMetric dist_metric):
    """
    Compute the necessary elements of a minimum spanning
    tree for computation of single linkage clustering. This
    represents the MST-LINKAGE-CORE algorithm (Figure 6) from
    :arxiv:`Daniel Mullner, "Modern hierarchical, agglomerative clustering
    algorithms" <1109.2378>`.

    In contrast to the scipy implementation is never computes
    a full distance matrix, generating distances only as they
    are needed and releasing them when no longer needed.

    Parameters
    ----------
    raw_data: array of shape (n_samples, n_features)
        The array of feature data to be clustered. Must be C-aligned

    dist_metric: DistanceMetric
        A DistanceMetric object conforming to the API from
        ``sklearn.metrics._dist_metrics.pxd`` that will be
        used to compute distances.

    Returns
    -------
    mst_core_data: array of shape (n_samples, 3)
        An array providing information from which one
        can either compute an MST, or the linkage hierarchy
        very efficiently. See :arxiv:`Daniel Mullner, "Modern hierarchical,
        agglomerative clustering algorithms" <1109.2378>` algorithm
        MST-LINKAGE-CORE for more details.
    """
    cdef:
        intp_t n_samples = raw_data.shape[0]
        uint8_t[:] in_tree = np.zeros(n_samples, dtype=bool)
        float64_t[:, ::1] result = np.zeros((n_samples - 1, 3))

        intp_t current_node = 0
        intp_t new_node
        intp_t i
        intp_t j
        intp_t num_features = raw_data.shape[1]

        float64_t right_value
        float64_t left_value
        float64_t new_distance

        float64_t[:] current_distances = np.full(n_samples, INFINITY)

    for i in range(n_samples - 1):

        in_tree[current_node] = 1

        new_distance = INFINITY
        new_node = 0

        for j in range(n_samples):
            if in_tree[j]:
                continue

            right_value = current_distances[j]
            left_value = dist_metric.dist(&raw_data[current_node, 0],
                                          &raw_data[j, 0],
                                          num_features)

            if left_value < right_value:
                current_distances[j] = left_value

            if current_distances[j] < new_distance:
                new_distance = current_distances[j]
                new_node = j

        result[i, 0] = current_node
        result[i, 1] = new_node
        result[i, 2] = new_distance
        current_node = new_node

    return np.array(result)
