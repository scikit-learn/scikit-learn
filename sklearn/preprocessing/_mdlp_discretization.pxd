# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
# cython: language_level=3


# Author: Juan Carlos Alfaro Jiménez <JuanCarlos.Alfaro@uclm.es>
# License: BSD


# These statements are required to use the C-API
# and to avoid segmentation faults, respectively
cimport numpy as np
np.import_array()


# =============================================================================
# Types
# =============================================================================

# Type for the data
ctypedef np.npy_float64 DTYPE_t

# Type for the labels
ctypedef np.npy_int32 ITYPE_t

# Type for the indices and counters
ctypedef np.npy_intp SIZE_t

# Type for boolean
ctypedef bint BOOL_t

# Add pointer types here
ctypedef fused PTR_t:
    (DTYPE_t*)
    (StackRecord*)


# =============================================================================
# Stack data structure
# =============================================================================

# A record on the stack
cdef struct StackRecord:
    SIZE_t start
    SIZE_t end
    DTYPE_t* sum_total


cdef class Stack:
    """A stack data structure."""

    cdef SIZE_t capacity  # Capacity of the stack
    cdef SIZE_t top  # Top record on the stack
    cdef StackRecord* stack  # The stack

    cdef int push(self,
                  SIZE_t start,
                  SIZE_t end,
                  DTYPE_t* sum_total,
                  SIZE_t sum_stride,
                  SIZE_t n_outputs) nogil except -1

    cdef int pop(self,
                 StackRecord* stack_record,
                 SIZE_t sum_stride,
                 SIZE_t n_outputs) nogil except -1
