from cython.operator cimport dereference as deref
from cpython.ref cimport Py_INCREF
cimport numpy as cnp

cnp.import_array()


cdef StdVectorSentinel _create_sentinel(vector_typed * vect_ptr):
    if vector_typed is vector[float64_t]:
        return StdVectorSentinelFloat64.create_for(vect_ptr)
    elif vector_typed is vector[int32_t]:
        return StdVectorSentinelInt32.create_for(vect_ptr)
    elif vector_typed is vector[int64_t]:
        return StdVectorSentinelInt64.create_for(vect_ptr)
    else:  # intp_t
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
    cdef vector[float64_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[float64_t] * vect_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelFloat64 sentinel = StdVectorSentinelFloat64.__new__(StdVectorSentinelFloat64)
        sentinel.vec.swap(deref(vect_ptr))
        return sentinel

    cdef void* get_data(self):
        return self.vec.data()

    cdef int get_typenum(self):
        return cnp.NPY_FLOAT64


cdef class StdVectorSentinelIntP(StdVectorSentinel):
    cdef vector[intp_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[intp_t] * vect_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelIntP sentinel = StdVectorSentinelIntP.__new__(StdVectorSentinelIntP)
        sentinel.vec.swap(deref(vect_ptr))
        return sentinel

    cdef void* get_data(self):
        return self.vec.data()

    cdef int get_typenum(self):
        return cnp.NPY_INTP


cdef class StdVectorSentinelInt32(StdVectorSentinel):
    cdef vector[int32_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[int32_t] * vect_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelInt32 sentinel = StdVectorSentinelInt32.__new__(StdVectorSentinelInt32)
        sentinel.vec.swap(deref(vect_ptr))
        return sentinel

    cdef void* get_data(self):
        return self.vec.data()

    cdef int get_typenum(self):
        return cnp.NPY_INT32


cdef class StdVectorSentinelInt64(StdVectorSentinel):
    cdef vector[int64_t] vec

    @staticmethod
    cdef StdVectorSentinel create_for(vector[int64_t] * vect_ptr):
        # This initializes the object directly without calling __init__
        # See: https://cython.readthedocs.io/en/latest/src/userguide/extension_types.html#instantiation-from-existing-c-c-pointers # noqa
        cdef StdVectorSentinelInt64 sentinel = StdVectorSentinelInt64.__new__(StdVectorSentinelInt64)
        sentinel.vec.swap(deref(vect_ptr))
        return sentinel

    cdef void* get_data(self):
        return self.vec.data()

    cdef int get_typenum(self):
        return cnp.NPY_INT64


cdef cnp.ndarray vector_to_nd_array(vector_typed * vect_ptr):
    cdef:
        cnp.npy_intp size = deref(vect_ptr).size()
        StdVectorSentinel sentinel = _create_sentinel(vect_ptr)
        cnp.ndarray arr = cnp.PyArray_SimpleNewFromData(
            1, &size, sentinel.get_typenum(), sentinel.get_data())

    # Makes the numpy array responsible of the life-cycle of its buffer.
    # A reference to the StdVectorSentinel will be stolen by the call to
    # `PyArray_SetBaseObject` below, so we increase its reference counter.
    # See: https://docs.python.org/3/c-api/intro.html#reference-count-details
    Py_INCREF(sentinel)
    cnp.PyArray_SetBaseObject(arr, sentinel)
    return arr
