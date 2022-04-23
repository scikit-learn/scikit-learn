#!python
cimport numpy as np

# Floating point/data type
ctypedef np.float64_t DTYPE_t  # WARNING: should match DTYPE in typedefs.pyx

cdef enum:
    DTYPECODE = np.NPY_FLOAT64
    ITYPECODE = np.NPY_INTP
    INT32TYPECODE = np.NPY_INT32
    INT64TYPECODE = np.NPY_INT64

# Index/integer type.
#  WARNING: ITYPE_t must be a signed integer type or you will have a bad time!
ctypedef np.intp_t ITYPE_t  # WARNING: should match ITYPE in typedefs.pyx
ctypedef np.int32_t INT32TYPE_t  # WARNING: should match INT32TYPE in typedefs.pyx
ctypedef np.int64_t INT64TYPE_t  # WARNING: should match INT32TYPE in typedefs.pyx
