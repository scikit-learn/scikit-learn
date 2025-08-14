from libc.time cimport timespec

cdef extern from *:
    """
    #include <string.h>
    #include <threads.h>
    static void __Pyx_once_flag_init(once_flag* flag) {
        once_flag other = ONCE_FLAG_INIT;
        memcpy(flag, &other, sizeof(once_flag));
    }
    """
    void once_flag_init "__Pyx_once_flag_init"(once_flag*)

cdef extern from "<threads.h>" nogil:
    ctypedef void (*_once_func_type)() noexcept nogil

    # Threads
    ctypedef struct thrd_t:
        pass
    # Declare as noexcept because you'll regret trying to throw a Python exception from it
    # and nogil because it's a new thread, so you definitely don't have a Python threadstate yet.
    ctypedef int (*thrd_start_t)(void*) noexcept nogil

    int thrd_create(thrd_t*, thrd_start_t, void*)
    int thrd_equal(thrd_t lhs, thrd_t rhs)
    thrd_t thrd_current()
    int thrd_sleep(const timespec* duration, timespec* remaining)
    void thrd_yield()
    void thrd_exit(int res)
    void thrd_detach(thrd_t thr)
    # Be very wary of deadlocks if calling thrd_join with the GIL.
    int thrd_join(thrd_t, int *res)
    enum:
        thrd_success
        thrd_nomem
        thrd_timedout
        thrd_busy
        thrd_error

    # Mutexes
    ctypedef struct mtx_t:
        pass
    int mtx_init(mtx_t* mtx, int type)
    int mtx_lock(mtx_t* mtx)
    int mtx_timed_lock(mtx_t* mtx, const timespec* timepoint)
    int mtx_trylock(mtx_t* mtx)
    int mtx_unlock(mtx_t* mtx)
    void mtx_destroy(mtx_t* mtx)
    enum:
        mtx_plain
        mtx_recursive
        mtx_timed

    # call once. Potentially this may be hard to use with Cython's code generation
    # (depending on the exact definition of "ONCE_FLAG_INIT").
    # Use the helper function "once_flag_init" to do it non-statically.
    ctypedef struct once_flag:
        pass
    once_flag ONCE_FLAG_INIT

    # We strongly recommend not calling this with the GIL held.
    void call_once(once_flag* flag, _once_func_type func) noexcept

    # condition variables
    ctypedef struct cnd_t:
        pass
    int cnd_init(cnd_t*)
    int cnd_signal(cnd_t*)
    int cnd_broadcast(cnd_t*)
    # Be very wary of calling wait/timedwait with the GIL held.
    int cnd_wait(cnd_t*, mtx_t*)
    int cnd_timedwait(cnd_t*, mtx_t*, const timespec*)
    void cnd_destroy(cnd_t*)

    # Thread-local storage
    ctypedef struct tss_t:
        pass
    enum:
        TSS_DTOR_ITERATIONS
    ctypedef void (*tss_dtor_t)(void*) noexcept nogil
    int tss_create(tss_t* key, tss_dtor_t destructor)
    void *tss_get(tss_t key)
    int tss_set(tss_t key, void* val)
    void tss_delete(tss_t key)
