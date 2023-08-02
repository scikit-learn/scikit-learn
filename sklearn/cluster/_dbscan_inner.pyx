# Fast inner loop for DBSCAN.
# Author: Lars Buitinck
# License: 3-clause BSD

from libcpp.vector cimport vector
cimport numpy as cnp

cnp.import_array()


def dbscan_inner(const cnp.uint8_t[::1] is_core,
                 object[:] neighborhoods,
                 cnp.npy_intp[::1] labels):
    cdef cnp.npy_intp i, label_num = 0, v
    cdef cnp.npy_intp[:] neighb
    cdef vector[cnp.npy_intp] stack

    for i in range(labels.shape[0]):
        if labels[i] != -1 or not is_core[i]:
            continue

        # Depth-first search starting from i, ending at the non-core points.
        # This is very similar to the classic algorithm for computing connected
        # components, the difference being that we label non-core points as
        # part of a cluster (component), but don't expand their neighborhoods.
        while True:
            if labels[i] == -1:
                labels[i] = label_num
                if is_core[i]:
                    neighb = neighborhoods[i]
                    for i in range(neighb.shape[0]):
                        v = neighb[i]
                        if labels[v] == -1:
                            stack.push_back(v)

            if stack.size() == 0:
                break
            i = stack.back()
            stack.pop_back()

        label_num += 1
