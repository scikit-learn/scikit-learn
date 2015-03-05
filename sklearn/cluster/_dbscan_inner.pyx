# Fast inner loop for DBSCAN.
# Author: Lars Buitinck
# License: 3-clause BSD

cimport cython
from libcpp.vector cimport vector
cimport numpy as np
import numpy as np


# Work around Cython bug: C++ exceptions are not caught unless thrown within
# a cdef function with an "except +" declaration.
cdef inline void push(vector[np.npy_intp] &stack, np.npy_intp u) except +:
    stack.push_back(u)


@cython.boundscheck(False)
@cython.wraparound(False)
def dbscan_inner(np.ndarray[np.uint8_t, ndim=1, mode='c'] is_core,
                 np.ndarray[object, ndim=1] neighborhoods,
                 np.ndarray[np.npy_intp, ndim=1, mode='c'] labels,
                 np.ndarray[np.npy_intp, ndim=1, mode='c'] index_order):
    cdef np.npy_intp i, label_num = 0, u, v
    cdef np.ndarray[np.npy_intp, ndim=1] neighb
    cdef vector[np.npy_intp] stack

    for i in range(index_order.shape[0]):
        u = index_order[i]
        if labels[u] != -1 or not is_core[u]:
            continue

        # Depth-first search starting from u, ending at the non-core points.
        # This is very similar to the classic algorithm for computing connected
        # components, the difference being that we label non-core points as
        # part of a cluster (component), but don't expand their neighborhoods.
        while True:
            if labels[u] == -1:
                labels[u] = label_num
                if is_core[u]:
                    neighb = neighborhoods[u]
                    for i in range(neighb.shape[0]):
                        v = neighb[i]
                        if labels[v] == -1:
                            push(stack, v)

            if stack.size() == 0:
                break
            u = stack.back()
            stack.pop_back()

        label_num += 1
