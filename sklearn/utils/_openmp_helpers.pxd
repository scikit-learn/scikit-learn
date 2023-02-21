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
        void omp_init_lock(struct omp_lock_t *lock) {}
        void omp_destroy_lock(struct omp_lock_t *lock) {}
        void omp_set_lock(struct omp_lock_t *lock) {}
        void omp_unset_lock(struct omp_lock_t *lock) {}
    #endif
    """
    bint USE_OPENMP
    ctypedef struct omp_lock_t:
        pass
    
    void omp_init_lock(struct omp_lock_t *lock) nogil
    void omp_destroy_lock(struct omp_lock_t *lock) nogil
    void omp_set_lock(struct omp_lock_t *lock) nogil
    void omp_unset_lock(struct omp_lock_t *lock) nogil

cdef int _openmp_thread_num() noexcept nogil
