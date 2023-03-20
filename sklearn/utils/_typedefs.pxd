#!python
cimport numpy as cnp

# Commonly used types
# These are redefinitions of the ones defined by numpy in
# https://github.com/numpy/numpy/blob/main/numpy/__init__.pxd
# and exposed by cython in
# https://github.com/cython/cython/blob/master/Cython/Includes/numpy/__init__.pxd.
# It will eventually avoid having to always include the numpy headers even when we
# would only use it for the types.
# TODO: don't cimport numpy in this extension.
#
# When used to declare variables that will receive values from numpy arrays, it
# should match the dtype of the array. For example, to declare a variable that will
# receive values from a numpy array of dtype np.float64, the type float64_t must be
# used.
#
# TODO: Stop defining custom types locally or globally like DTYPE_t and friends and
# use these consistently throughout the codebase.
# NOTE: Extend this list as needed when converting more cython extensions.
ctypedef unsigned char bool_t
ctypedef Py_ssize_t intp_t
ctypedef float float32_t
ctypedef double float64_t
ctypedef signed int int32_t
ctypedef signed long long int64_t


# Floating point/data type
ctypedef cnp.float64_t DTYPE_t  # WARNING: should match DTYPE in typedefs.pyx

cdef enum:
    DTYPECODE = cnp.NPY_FLOAT64

# Index/integer type.
#  WARNING: ITYPE_t must be a signed integer type or you will have a bad time!
ctypedef cnp.intp_t ITYPE_t  # WARNING: should match ITYPE in typedefs.pyx

# scipy matrices indices dtype (namely for indptr and indices arrays)
#
#   Note that indices might need to be represented as cnp.int64_t.
#   Currently, we use Cython classes which do not handle fused types
#   so we hardcode this type to cnp.int32_t, supporting all but edge
#   cases.
#
# TODO: support cnp.int64_t for this case
# See: https://github.com/scikit-learn/scikit-learn/issues/23653
ctypedef cnp.int32_t SPARSE_INDEX_TYPE_t
