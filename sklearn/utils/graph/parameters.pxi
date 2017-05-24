
DTYPE = np.float64
ctypedef np.float64_t DTYPE_t

ITYPE = np.int32
ctypedef np.int32_t ITYPE_t

# Fused type for int32 and int64
ctypedef fused int32_or_int64:
    np.int32_t
    np.int64_t

# EPS is the precision of DTYPE
DEF DTYPE_EPS = 1E-15

# NULL_IDX is the index used in predecessor matrices to store a non-path
DEF NULL_IDX = -9999
