import numpy as np
cimport numpy as np

ctypedef np.npy_intp SIZE_t              # Type for indices and counters
ctypedef SIZE_t* LEVEL
ctypedef np.float64_t FLOAT

cdef inline void set_level(LEVEL level, SIZE_t start,
                           SIZE_t end, SIZE_t depth):
    level[0] = start
    level[1] = end
    level[2] = depth

#cdef class SearchLevel:
#    cdef int start, end, depth
#
#    def inline __cinit__(self, start, end, depth):
#        self.start = start
#        self.end = end
        #self.depth = depth

