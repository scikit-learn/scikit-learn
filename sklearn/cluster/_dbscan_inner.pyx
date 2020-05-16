# Fast inner loop for DBSCAN.
# Author: Lars Buitinck
# License: 3-clause BSD
#
# cython: boundscheck=False, wraparound=False

cimport cython
from libcpp.vector cimport vector
cimport numpy as np
import numpy as np

np.import_array()


# Work around Cython bug: C++ exceptions are not caught unless thrown within
# a cdef function with an "except +" declaration.
cdef inline void push(vector[np.npy_intp] &stack, np.npy_intp i) except +:
    stack.push_back(i)


def dbscan_inner(np.ndarray[np.uint8_t, ndim=1, mode='c'] is_core,
                 np.ndarray[object, ndim=1] neighborhoods,
                 np.ndarray[np.npy_intp, ndim=1, mode='c'] labels):
    cdef np.npy_intp i, label_num = 0, v
    cdef np.ndarray[np.npy_intp, ndim=1] neighb
    cdef vector[np.npy_intp] stack

    for i in range(labels.shape[0]):
        if labels[i] == -1 and not is_core[i]:

            # Depth-first search starting from i, ending at the non-core points.
            # This is very similar to the classic algorithm for computing connected
            # components, the difference being that we label non-core points as
            # part of a cluster (component), but don't expand their neighborhoods.
            while True:
                labels[i] = label_num
                for j in range(neighborhoods[i].shape[0]):
                    if labels[neighborhoods[i][j]] == -1:
                        push(stack, neighborhoods[i][j])

                if stack.size() == 0:
                    break
                i = stack.back()
                stack.pop_back()

            label_num += 1
