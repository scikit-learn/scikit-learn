from cython.operator cimport dereference as deref
from cpython.ref cimport Py_INCREF
cimport numpy as np

from ._typedefs cimport DTYPECODE, ITYPECODE

np.import_array()


cdef StdVectorSentinel _create_sentinel(vector_typed * vect_ptr):
    if vector_typed is vector[DTYPE_t]:
        return StdVectorSentinelFloat64.create_for(vect_ptr)
    else:
        return StdVectorSentinelIntP.create_for(vect_ptr)


cdef class StdVectorSentinel:
    """Wraps a reference to a vector which will be deallocated with this object.

    When created, the StdVectorSentinel swaps the reference of its internal
    vectors with the provided one (vect_ptr), thus making the StdVectorSentinel
    manage the provided one's lifetime.
    """
    cdef void* get_data(self):
        """Return pointer to data."""

    cdef int get_typenum(self):
        """Get typenum for PyArray_SimpleNewFromData."""


cdef class StdVectorSentinelFloat64(StdVectorSentinel):
    cdef vector[DTYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[DTYPE_t] * vect_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelFloat64 sentinel = StdVectorSentinelFloat64.__new__(StdVectorSentinelFloat64)
        sentinel.vec.swap(deref(vect_ptr))
        return sentinel

    cdef void* get_data(self):
        return self.vec.data()

    cdef int get_typenum(self):
        return DTYPECODE


cdef class StdVectorSentinelIntP(StdVectorSentinel):
    cdef vector[ITYPE_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[ITYPE_t] * vect_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelIntP sentinel = StdVectorSentinelIntP.__new__(StdVectorSentinelIntP)
        sentinel.vec.swap(deref(vect_ptr))
        return sentinel

    cdef void* get_data(self):
        return self.vec.data()

    cdef int get_typenum(self):
        return ITYPECODE


cdef np.ndarray vector_to_nd_array(vector_typed * vect_ptr):
    cdef:
        np.npy_intp size = deref(vect_ptr).size()
        StdVectorSentinel sentinel =  _create_sentinel(vect_ptr)
        np.ndarray arr = np.PyArray_SimpleNewFromData(
            1, &size, sentinel.get_typenum(), sentinel.get_data())

    # Makes the numpy array responsible of the life-cycle of its buffer.
    # A reference to the StdVectorSentinel will be stolen by the call to
    # `PyArray_SetBaseObject` below, so we increase its reference counter.
    # See: https://docs.python.org/3/c-api/intro.html#reference-count-details
    Py_INCREF(sentinel)
    np.PyArray_SetBaseObject(arr, sentinel)
    return arr
