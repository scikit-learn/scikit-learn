# Fast inner loop for DBSCAN.
# Author: Lars Buitinck
# License: 3-clause BSD

cimport cython
from libcpp.vector cimport vector
cimport numpy as np
import numpy as np
from scipy import sparse


# Work around Cython bug: C++ exceptions are not caught unless thrown within
# a cdef function with an "except +" declaration.
cdef inline void push(vector[np.npy_intp] &stack, np.npy_intp i) except +:
    stack.push_back(i)

def get_neighborhood(X, neighbors_model, eps, sample_weight, min_samples):
    def query_nn(i):
        if sparse.issparse(X):
            pt = X[i]
        else:
            pt = X[i].reshape(1, -1)

        pt_neighborhoods = neighbors_model.radius_neighbors(pt, eps, return_distance=False)
        if sample_weight is None:
            n_neighbors = len(pt_neighborhoods[0])
        else:
            n_neighbors = np.sum(sample_weight[pt_neighborhoods[0]])

        # Is current point a core sample?
        if n_neighbors >= min_samples:
            is_core = 1
        else:
            is_core = 0

        return (pt_neighborhoods[0], is_core)

    return query_nn


@cython.boundscheck(False)
@cython.wraparound(False)
def dbscan_inner(np.ndarray[np.uint8_t, ndim=1, mode='c'] is_core,
                 np.ndarray[object, ndim=1] neighborhoods,
                 np.ndarray[np.npy_intp, ndim=1, mode='c'] labels,
                 save_memory,
                 get_neighborhood):
    cdef np.npy_intp i, label_num = 0, v
    cdef np.ndarray[np.npy_intp, ndim=1] neighb
    cdef vector[np.npy_intp] stack

    for i in range(labels.shape[0]):
        if save_memory:
            if labels[i] != -1 or is_core[i] == 0:
                continue

            # We don't know if the point is a core sample or not. Let's query it.
            if is_core[i] == 2:
                # Query nearest neighbors of the point.
                (pt_neighborhoods, is_core[i]) = get_neighborhood(i)

        if labels[i] != -1 or not is_core[i]:
            continue
        # Depth-first search starting from i, ending at the non-core points.
        # This is very similar to the classic algorithm for computing connected
        # components, the difference being that we label non-core points as
        # part of a cluster (component), but don't expand their neighborhoods.
        while True:
            if labels[i] == -1:
                labels[i] = label_num

                if is_core[i] == 1:
                    if save_memory:
                        neighb = pt_neighborhoods
                    else:
                        neighb = neighborhoods[i]

                    for i in range(neighb.shape[0]):
                        v = neighb[i]
                        if labels[v] == -1:
                            push(stack, v)

            if stack.size() == 0:
                break
            i = stack.back()

            if save_memory and labels[i] == -1 and is_core[i] == 2:
                (pt_neighborhoods, is_core[i]) = get_neighborhood(i)

            stack.pop_back()

        label_num += 1
