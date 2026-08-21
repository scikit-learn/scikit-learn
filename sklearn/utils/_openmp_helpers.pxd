# Helpers to safely access OpenMP routines
#
# no-op implementations are provided for the case where OpenMP is not available.
#
# All calls to OpenMP routines should be cimported from this module.

cdef extern from *:
    """
    #ifdef _OPENMP
        #include <omp.h>
        #define SKLEARN_OPENMP_PARALLELISM_ENABLED 1
    #else
        #define SKLEARN_OPENMP_PARALLELISM_ENABLED 0
        #define omp_lock_t int
        #define omp_init_lock(l) (void)0
        #define omp_destroy_lock(l) (void)0
        #define omp_set_lock(l) (void)0
        #define omp_unset_lock(l) (void)0
        #define omp_get_thread_num() 0
        #define omp_get_max_threads() 1
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


cdef inline bint _use_threads_for_workload(
        int n_threads, long long n_work_items, int ops_per_item,
        long long min_instructions_per_thread) noexcept nogil:
    # Only parallelize a loop if each of the n_threads threads would get
    # more than min_instructions_per_thread simple instructions worth of
    # work (by default ~2000, i.e. ~1us, e.g. 400 samples at ~5 ops/sample):
    # below that, the cost of dispatching/synchronizing an OpenMP team
    # dominates the actual work, and running single-threaded is faster.
    # This is meant to feed `use_threads_if`, which compiles down to a plain
    # `if(...)` clause on the OpenMP construct, so no thread team is created
    # at all when it evaluates to false.
    #
    # min_instructions_per_thread should be obtained by the caller via
    # _min_instructions_per_thread(n_threads) (in _openmp_helpers.pyx), which
    # caches it, per power-of-2 n_threads bucket, after calibrating it once
    # on first use.
    return (n_threads > 1) and (
        n_work_items * <long long> ops_per_item
        >= min_instructions_per_thread * <long long> n_threads
    )
