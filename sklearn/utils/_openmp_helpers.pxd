# Helpers to access OpenMP threads information
#
# Those interfaces act as indirections which allows the non-support of OpenMP
# for implementations which have been written for it.

cdef extern from *:
    """
    #ifdef _OPENMP
        #define USE_OPENMP 1
        #include <omp.h>
    #else
        #define USE_OPENMP 0
        struct omp_lock_t {
            int dummy;
        };
    #endif
    """
    bint USE_OPENMP
    ctypedef struct omp_lock_t:
        pass

cdef int _openmp_thread_num() noexcept nogil
