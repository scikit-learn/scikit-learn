from ...utils._typedefs cimport intp_t, float64_t, uint8_t

# This corresponds to the scipy.cluster.hierarchy format
ctypedef packed struct HIERARCHY_t:
    intp_t left_node
    intp_t right_node
    float64_t value
    intp_t cluster_size

# Effectively an edgelist encoding a parent/child pair, along with a value and
# the corresponding cluster_size in each row providing a tree structure.
ctypedef packed struct CONDENSED_t:
    intp_t parent
    intp_t child
    float64_t value
    intp_t cluster_size

cdef extern from "numpy/arrayobject.h":
    ctypedef struct PyArrayObject
    intp_t * PyArray_SHAPE(PyArrayObject *)
