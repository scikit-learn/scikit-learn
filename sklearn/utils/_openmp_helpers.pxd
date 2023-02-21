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
        typedef void omp_lock_t;

        void omp_init_lock(omp_lock_t *lock) {}
        void omp_destroy_lock(omp_lock_t *lock) {}
        void omp_set_lock(omp_lock_t *lock) {}
        void omp_unset_lock(omp_lock_t *lock) {}
    #endif
    """
    bint USE_OPENMP
    ctypedef struct omp_lock_t:
        pass

    void omp_init_lock(omp_lock_t*) nogil
    void omp_destroy_lock(omp_lock_t*) nogil
    void omp_set_lock(omp_lock_t*) nogil
    void omp_unset_lock(omp_lock_t*) nogil

cdef int _openmp_thread_num() noexcept nogil
