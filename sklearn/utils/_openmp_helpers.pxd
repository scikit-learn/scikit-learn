# Helpers to access OpenMP threads information
#
# Those interfaces act as indirections which allows the non-support of OpenMP
# for implementations which have been written for it.

cdef extern from *:
    """
    #ifdef _OPENMP
        #include <omp.h>
        #define SKLEARN_OPENMP_PARALLELISM_ENABLED 1
    #else
        #include <_fake_openmp.h>
        #define SKLEARN_OPENMP_PARALLELISM_ENABLED 0
    #endif
    """
    bint SKLEARN_OPENMP_PARALLELISM_ENABLED

    ctypedef struct omp_lock_t:
        pass

    void omp_init_lock(omp_lock_t*) noexcept nogil
    void omp_destroy_lock(omp_lock_t*) noexcept nogil
    void omp_set_lock(omp_lock_t*) noexcept nogil
    void omp_unset_lock(omp_lock_t*) noexcept nogil

    int omp_get_thread_num() noexcept nogil
    int omp_get_max_threads() noexcept nogil
