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
