from libc.string cimport memcpy

# Commonly used types
# These are redefinitions of the ones defined by numpy in
# https://github.com/numpy/numpy/blob/main/numpy/__init__.pxd.
# It will eventually avoid having to always include the numpy headers even when we
# would only use it for the types.
#
# When used to declare variables that will receive values from numpy arrays, it
# should match the dtype of the array. For example, to declare a variable that will
# receive values from a numpy array of dtype np.float64, the type float64_t must be
# used.
#
# TODO: Stop defining custom types locally or globally like DTYPE_t and friends and
# use these consistently throughout the codebase.
# NOTE: Extend this list as needed when converting more cython extensions.
ctypedef unsigned char uint8_t
ctypedef unsigned int uint32_t
ctypedef unsigned long long uint64_t
# Note: In NumPy 2, indexing always happens with npy_intp which is an alias for
# the Py_ssize_t type, see PEP 353.
#
# Note that on most platforms Py_ssize_t is equivalent to C99's intptr_t,
# but they can differ on architecture with segmented memory (none
# supported by scikit-learn at the time of writing).
#
# intp_t/np.intp should be used to index arrays in a platform dependent way.
# Storing arrays with platform dependent dtypes as attribute on picklable
# objects is not recommended as it requires special care when loading and
# using such datastructures on a host with different bitness. Instead one
# should rather use fixed width integer types such as int32 or uint32 when we know
# that the number of elements to index is not larger to 2 or 4 billions.
ctypedef Py_ssize_t intp_t
ctypedef float float32_t
ctypedef double float64_t
# Sparse matrices indices and indices' pointers arrays must use int32_t over
# intp_t because intp_t is platform dependent.
# When large sparse matrices are supported, indexing must use int64_t.
# See https://github.com/scikit-learn/scikit-learn/issues/23653 which tracks the
# ongoing work to support large sparse matrices.
ctypedef signed char int8_t
ctypedef signed int int32_t
ctypedef signed long long int64_t

ctypedef fused float32_or_float64_t:
    float32_t
    float64_t


cdef inline bint isnan(float32_or_float64_t x) noexcept nogil:
    """Check whether x is NaN.

    Prefer this over libc.math.isnan in hot loops: unlike that libm call,
    which some compilers/libc fail to inline, this is guaranteed to be
    inlined. See https://github.com/scikit-learn/scikit-learn/issues/34869.
    """
    cdef uint32_t bits32
    cdef uint64_t bits64
    if float32_or_float64_t is float32_t:
        memcpy(&bits32, &x, sizeof(bits32))
        return (bits32 & 0x7fffffff) > 0x7f800000
    else:
        memcpy(&bits64, &x, sizeof(bits64))
        return (bits64 & <uint64_t>0x7fffffffffffffff) > <uint64_t>0x7ff0000000000000
