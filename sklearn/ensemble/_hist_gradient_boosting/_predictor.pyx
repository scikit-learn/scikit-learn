# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3

# Author: Nicolas Hug

cimport cython
from cython.parallel import prange
import numpy as np
cimport numpy as np

from .types cimport X_DTYPE_C
from .types cimport Y_DTYPE_C
from .types cimport X_BINNED_DTYPE_C


cdef packed struct node_struct:
    # Equivalent struct to PREDICTOR_RECORD_DTYPE to use in memory views. It
    # needs to be packed since by default numpy dtypes aren't aligned
    Y_DTYPE_C value
    unsigned int count
    unsigned int feature_idx
    X_DTYPE_C threshold
    unsigned int left
    unsigned int right
    Y_DTYPE_C gain
    unsigned int depth
    unsigned char is_leaf
    X_BINNED_DTYPE_C bin_threshold


def _predict_from_numeric_data(nodes, numeric_data, out):
    _predict_from_numeric_data_parallel(nodes, numeric_data, out)


def _predict_from_binned_data(nodes, binned_data, out):
    _predict_from_binned_data_parallel(nodes, binned_data, out)


cdef void _predict_from_numeric_data_parallel(
        node_struct [:] nodes,
        const X_DTYPE_C [:, :] numeric_data,
        Y_DTYPE_C [:] out):

    cdef:
        int i

    for i in prange(numeric_data.shape[0], schedule='static', nogil=True):
        out[i] = _predict_one_from_numeric_data(nodes, numeric_data, i)


cdef inline Y_DTYPE_C _predict_one_from_numeric_data(
        node_struct [:] nodes,
        const X_DTYPE_C [:, :] numeric_data,
        const int row) nogil:
    # Need to pass the whole array and the row index, else prange won't work.
    # See issue Cython #2798

    cdef:
        node_struct node = nodes[0]

    while True:
        if node.is_leaf:
            return node.value
        if numeric_data[row, node.feature_idx] <= node.threshold:
            node = nodes[node.left]
        else:
            node = nodes[node.right]


cdef void _predict_from_binned_data_parallel(
        node_struct [:] nodes,
        const X_BINNED_DTYPE_C [:, :] binned_data,
        Y_DTYPE_C [:] out):

    cdef:
        int i

    for i in prange(binned_data.shape[0], schedule='static', nogil=True):
        out[i] = _predict_one_from_binned_data(nodes, binned_data, i)


cdef inline Y_DTYPE_C _predict_one_from_binned_data(
        node_struct [:] nodes,
        const X_BINNED_DTYPE_C [:, :] binned_data,
        const int row) nogil:
    # Need to pass the whole array and the row index, else prange won't work.
    # See issue Cython #2798

    cdef:
        node_struct node = nodes[0]

    while True:
        if node.is_leaf:
            return node.value
        if binned_data[row, node.feature_idx] <= node.bin_threshold:
            node = nodes[node.left]
        else:
            node = nodes[node.right]

def _compute_partial_dependence(
    node_struct [:] nodes,
    const X_DTYPE_C [:, ::1] X,
    int [:] target_feature,
    Y_DTYPE_C [:] out):

    cdef:

        unsigned int current_node_idx
        unsigned int [:] node_idx_stack = np.zeros(shape=nodes.shape[0],
                                                   dtype=np.uint32)
        double[::1] weight_stack = np.zeros(shape=nodes.shape[0],
                                            dtype=np.double)
        node_struct * current_node  # pointer to avoid copying attributes

        unsigned int sample_idx
        unsigned feature_idx
        unsigned stack_size
        double left_sample_frac
        double current_weight
        double total_weight  # used for sanity check only
        bint is_target_feature

    for sample_idx in range(X.shape[0]):
        # init stacks for current sample
        stack_size = 1
        node_idx_stack[0] = 0  # root node
        weight_stack[0] = 1  # all the samples are in the root node
        total_weight = 0

        while stack_size > 0:

            # pop the stack
            stack_size -= 1
            current_node_idx = node_idx_stack[stack_size]
            current_node = &nodes[current_node_idx]

            if current_node.is_leaf:
                out[sample_idx] += (weight_stack[stack_size] *
                                    current_node.value)
                total_weight += weight_stack[stack_size]
            else:
                # determine if the split feature is a target feature
                is_target_feature = False
                for feature_idx in range(target_feature.shape[0]):
                    if target_feature[feature_idx] == current_node.feature_idx:
                        is_target_feature = True
                        break

                if is_target_feature:
                    # In this case, we push left or right child on stack
                    if X[sample_idx, feature_idx] <= current_node.threshold:
                        node_idx_stack[stack_size] = current_node.left
                    else:
                        node_idx_stack[stack_size] = current_node.right
                    stack_size += 1
                else:
                    # In this case, we push both children onto the stack,
                    # and give a weight proportional to the number of
                    # samples going through each branch.

                    # push left child
                    node_idx_stack[stack_size] = current_node.left
                    left_sample_frac = (
                        nodes[current_node.left].count / current_node.count)
                    current_weight = weight_stack[stack_size]
                    weight_stack[stack_size] = current_weight * left_sample_frac
                    stack_size += 1

                    # push right child
                    node_idx_stack[stack_size] = current_node.right
                    weight_stack[stack_size] = (
                        current_weight * (1 - left_sample_frac))
                    stack_size += 1

        # Sanity check. Should never happen.
        if not (0.999 < total_weight < 1.001):
            raise ValueError("Total weight should be 1.0 but was %.9f" %
                                total_weight)
