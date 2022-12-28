#!python
cimport numpy as cnp

# Floating point/data type
ctypedef cnp.float64_t DTYPE_t  # WARNING: should match DTYPE in typedefs.pyx

cdef enum:
    DTYPECODE = cnp.NPY_FLOAT64
    ITYPECODE = cnp.NPY_INTP
    INT32TYPECODE = cnp.NPY_INT32
    INT64TYPECODE = cnp.NPY_INT64

# Index/integer type.
#  WARNING: ITYPE_t must be a signed integer type or you will have a bad time!
ctypedef cnp.intp_t ITYPE_t  # WARNING: should match ITYPE in typedefs.pyx
ctypedef cnp.int32_t INT32TYPE_t  # WARNING: should match INT32TYPE in typedefs.pyx
ctypedef cnp.int64_t INT64TYPE_t  # WARNING: should match INT32TYPE in typedefs.pyx

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
