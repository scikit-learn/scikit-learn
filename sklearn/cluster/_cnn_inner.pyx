# Fast inner loop for CNN.
# Author: Jan-Oliver Joswig
# License: 3-clause BSD
#
# cython: boundscheck = False
# cython: wraparound = False

cimport cython
from libcpp.queue cimport queue as cppqueue
cimport numpy as np


ctypedef np.intp_t ARRAYINDEX_DTYPE_t


cdef inline bint check_similarity(
        ARRAYINDEX_DTYPE_t[::1] a, ARRAYINDEX_DTYPE_t sa,
        ARRAYINDEX_DTYPE_t[::1] b,
        ARRAYINDEX_DTYPE_t c):
    """Check if the CNN criterion is fullfilled

    Check if `a` and `b` have at least `c` common elements.  Faster than
    computing the intersection (say with `numpy.intersect1d`) and
    comparing its length to `c`.

    Args:
        a: 1D Array of integers (neighbor point indices) that
           supports the buffer protocol.
        b: 1D Array of integers (neighbor point indices) that
           supports the buffer protocol.
        sa: Length of `a`. Only received here because already computed.
        c: Check if `a` and `b` share this many elements.

    Returns:
        True (1) or False (0)
    """

    cdef ARRAYINDEX_DTYPE_t i, j, sb = b.shape[0]  # Control variables
    cdef ARRAYINDEX_DTYPE_t ai, bj                 # Checked elements
    cdef ARRAYINDEX_DTYPE_t common = 0             # Common neighbors count

    if c == 0:
        return 1

    for i in range(sa):
        ai = a[i]
        for j in range(sb):
            bj = b[j]
            if ai == bj:
                # Check similarity and return/move on early
                common += 1
                if common == c:
                    return 1
                break
    return 0


def cnn_inner(
        object[::1] neighborhoods,
        ARRAYINDEX_DTYPE_t[::1] labels,
        np.uint8_t[::1] core_candidates,
        ARRAYINDEX_DTYPE_t min_samples):

    cdef ARRAYINDEX_DTYPE_t init_point, point, member, member_i
    cdef ARRAYINDEX_DTYPE_t m, n = neighborhoods.shape[0]
    cdef ARRAYINDEX_DTYPE_t[::1] neighbors, neighbor_neighbors
    cdef ARRAYINDEX_DTYPE_t current = 0  # Cluster (start at 0; noise = -1)
    cdef unsigned long membercount       # Current cluster size
    cdef cppqueue[ARRAYINDEX_DTYPE_t] q  # FIFO queue

    # BFS find connected components
    for init_point in range(n):
        # Loop over points and find source node
        if core_candidates[init_point] == 0:
            # Point already assigned or not enough neighbors
            continue
        core_candidates[init_point] = 0  # Mark point as included

        neighbors = neighborhoods[init_point]
        m = neighbors.shape[0]

        labels[init_point] = current     # Assign cluster label
        membercount = 1

        while True:
            for member_i in range(m):
                # Loop over connected neighboring points
                member = neighbors[member_i]
                if core_candidates[member] == 0:
                    # Point already assigned or not enough neighbors
                    continue

                neighbor_neighbors = neighborhoods[member]
                if check_similarity(             # Conditional growth
                        neighbors, m, neighbor_neighbors, min_samples):
                    core_candidates[member] = 0  # Point included
                    labels[member] = current     # Assign cluster label
                    membercount += 1             # Cluster grows
                    q.push(member)               # Add point to queue

            if q.empty():
                # No points left to check
                if membercount == 1:
                    # Cluster to small -> effectively noise
                    labels[init_point] = -1  # Revert assignment
                    current -= 1             # Revert cluster number
                break

            point = q.front()  # Get the next point from the queue
            q.pop()

            neighbors = neighborhoods[point]
            m = neighbors.shape[0]

        current += 1  # Increase cluster number
