# cython: preliminary_late_includes_cy28=True

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

    # Mutexes.
    # We strongly recommend not calling mtx_lock with the GIL held and being
    # very wary of calling mtx_timed_lock with the GIL (because it will block the GIL until the timeout).
    ctypedef struct mtx_t:
        pass
    int mtx_init(mtx_t* mtx, int type)
    int mtx_lock(mtx_t* mtx)  # see py_safe_mutex_lock for GIL-safe version
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
    # We define some GIL-friendly wrappers for C++ call-once which you can use
    # in these circumstances.
    # py_safe_call_once in this file avoids blocking, but the function
    # must still be nogil.
    void call_once(once_flag* flag, _once_func_type func) noexcept

    # condition variables
    ctypedef struct cnd_t:
        pass
    int cnd_init(cnd_t*)
    int cnd_signal(cnd_t*)
    int cnd_broadcast(cnd_t*)
    # Be very wary of calling wait/timedwait with the GIL held. Two py_safe_ versions exist.
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


cdef inline void _dummy_force_utility_code_inclusion() nogil:
    with nogil:
        pass

cdef extern from *:
    ctypedef struct py_safe_once_flag "once_flag":
        pass
    py_safe_once_flag PY_SAFE_ONCE_FLAG_INIT "ONCE_FLAG_INIT"

cdef extern from *:
    """
    static PyGILState_STATE __pyx_libc_threads_limited_api_ensure_gil() {
        #if CYTHON_COMPILING_IN_LIMITED_API
        if ((__PYX_LIMITED_VERSION_HEX < 0x030d0000) && __Pyx_get_runtime_version() < 0x030d0000) {
            return PyGILState_Ensure();
        }
        #endif
        return PyGILState_LOCKED;  // Unused
    }

    #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
    static void __pyx_libc_threads_limited_api_release_gil(PyGILState_STATE gil_state) {
        if (__Pyx_get_runtime_version() < 0x030d0000)
            PyGILState_Release(gil_state);
    }
    #else
    #define __pyx_libc_threads_limited_api_release_gil(ignore) (void)ignore
    #endif

    static int __pyx_libc_threads_has_gil() {
        #if CYTHON_COMPILING_IN_LIMITED_API
            if ((__PYX_LIMITED_VERSION_HEX >= 0x030d0000) || __Pyx_get_runtime_version() >= 0x030d0000) {
                // In 3.13+ we can temporarily give up the GIL to find out what the thread state was.
                PyThreadState *ts = PyThreadState_Swap(NULL);
                if (ts) {
                    PyThreadState_Swap(ts);
                    return 1;
                }
                return 0;
            }
            /* There is no way to know if we have the GIL. Therefore the only
             * thing we can safely do is make absolutely sure that we have it
             * (in __pyx_libc_threads_limited_api_ensure_gil).
             */
            return 1;
        #elif PY_VERSION_HEX >= 0x030d0000
            return PyThreadState_GetUnchecked() != NULL;
        #elif PY_VERSION_HEX >= 0x030C0000
            return _PyThreadState_UncheckedGet() != NULL;
        #else
            return PyGILState_Check();
        #endif
    }

    static int __pyx_py_safe_mtx_lock_release_lock_reacquire(mtx_t* mutex) {
        // Release the GIL, acquire the lock, then reacquire the GIL.
        // This is safe provided the user never holds the GIL while trying
        // to reacquire the lock (i.e. it's safe provided they always use
        // the py-safe wrappers).
        int result;
        Py_BEGIN_ALLOW_THREADS
        result = mtx_lock(mutex);
        Py_END_ALLOW_THREADS
        return result;
    }

    static int __pyx_py_safe_mtx_lock_impl(mtx_t* mutex) {
        int lock_result = mtx_trylock(mutex);
        if (lock_result == thrd_busy) {
            return __pyx_py_safe_mtx_lock_release_lock_reacquire(mutex);
        }
        return lock_result;
    }

    CYTHON_UNUSED static int __pyx_py_safe_mtx_lock(mtx_t* mutex) {
        PyGILState_STATE gil_state = __pyx_libc_threads_limited_api_ensure_gil();
        if (!__pyx_libc_threads_has_gil())
            return mtx_lock(mutex); /* No GIL, no problem */
        int result = __pyx_py_safe_mtx_lock_impl(mutex);
        __pyx_libc_threads_limited_api_release_gil(gil_state);
        return result;
    }

    CYTHON_UNUSED static int __pyx_py_safe_cnd_wait( cnd_t* cond, mtx_t* mutex) {
        __Pyx_UnknownThreadState thread_state = __Pyx_SaveUnknownThread();
        int result = cnd_wait(cond, mutex);
        if (__Pyx_UnknownThreadStateMayHaveHadGil(thread_state)) {
            int unlock_result, relock_result;
            unlock_result = mtx_unlock(mutex);
            __Pyx_RestoreUnknownThread(thread_state);
            PyGILState_STATE gil_state = __pyx_libc_threads_limited_api_ensure_gil();
            relock_result = __pyx_py_safe_mtx_lock_impl(mutex);
            __pyx_libc_threads_limited_api_release_gil(gil_state);
            return relock_result != thrd_success ? relock_result :
                unlock_result != thrd_success ? unlock_result :
                result;
        } else {
            __Pyx_RestoreUnknownThread(thread_state);
            return result;
        }
    }

    CYTHON_UNUSED static int __pyx_py_safe_cnd_timedwait(cnd_t* cond, mtx_t* mutex, const struct timespec* time_point) {
        __Pyx_UnknownThreadState thread_state = __Pyx_SaveUnknownThread();
        int result = cnd_timedwait(cond, mutex, time_point);
        if (__Pyx_UnknownThreadStateMayHaveHadGil(thread_state)) {
            int unlock_result, relock_result;
            unlock_result = mtx_unlock(mutex);
            __Pyx_RestoreUnknownThread(thread_state);
            PyGILState_STATE gil_state = __pyx_libc_threads_limited_api_ensure_gil();
            relock_result = __pyx_py_safe_mtx_lock_impl(mutex);
            __pyx_libc_threads_limited_api_release_gil(gil_state);
            return relock_result != thrd_success ? relock_result :
                unlock_result != thrd_success ? unlock_result :
                result;
        } else {
            __Pyx_RestoreUnknownThread(thread_state);
            return result;
        }
    }

    CYTHON_UNUSED static void __pyx_libc_threads_py_safe_call_once(once_flag* flag, void (*func)(void)) {
        __Pyx_UnknownThreadState thread_state = __Pyx_SaveUnknownThread();
        call_once(flag, func);
        __Pyx_RestoreUnknownThread(thread_state);
    }
    """
    # Acquire the lock without deadlocking on the GIL.
    # At the end, the GIL will be restored to its original state.
    int py_safe_mtx_lock "__pyx_py_safe_mtx_lock"(mtx_t* mutex) nogil

    # Wait on the condition variable without deadlocking on the GIL.
    # At the end, the GIL will be restored to its original state.
    int py_safe_cnd_wait "__pyx_py_safe_cnd_wait"(cnd_t*, mtx_t*) nogil
    int py_safe_cnd_timedwait "__pyx_py_safe_cnd_timedwait"(cnd_t*, mtx_t*, const timespec*) nogil

    # call_once making sure not to deadlock on the GIL.
    # The callable must still be nogil though.
    void py_safe_once_flag_init "__Pyx_once_flag_init"(py_safe_once_flag*)
    int py_safe_call_once "__pyx_libc_threads_py_safe_call_once"(py_safe_once_flag* flag, _once_func_type func) nogil
