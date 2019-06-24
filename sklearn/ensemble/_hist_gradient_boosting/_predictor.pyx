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
from .types cimport node_struct


def _predict_from_numeric_data(
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


def _predict_from_binned_data(
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
