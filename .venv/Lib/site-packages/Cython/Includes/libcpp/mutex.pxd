# cython: preliminary_late_includes_cy28=True

from libcpp cimport bool

cdef extern from "<mutex>" namespace "std" nogil:
    # For all these mutex classes, we strongly recommend you do not use any
    # blocking lock function while holding the GIL (try_lock should be fine though).
    cppclass mutex:
        # may not be present, and we know nothing about it
        cppclass native_handle_type:
            pass

        void lock() except+
        bool try_lock()
        void unlock()

        native_handle_type native_handle() except+

    cppclass timed_mutex:
        # may not be present, and we know nothing about it
        cppclass native_handle_type:
            pass

        void lock() except+
        bool try_lock()
        bool try_lock_for[T](const T& duration) except+
        bool try_lock_until[T](const T& time_point) except+
        void unlock()

        native_handle_type native_handle() except+

    # We strongly recommend not mixing recursive_mutex and the GIL at all.
    # Because "unlock" may not actually unlock it, it's pretty hard to reason about
    # avoiding deadlocks.
    cppclass recursive_mutex:
        # may not be present, and we know nothing about it
        cppclass native_handle_type:
            pass

        void lock() except+
        bool try_lock()
        void unlock()

        native_handle_type native_handle() except+

    # We strongly recommend not mixing timed_recursive_mutex and the GIL at all.
    # Because "unlock" may not actually unlock it, it's pretty hard to reason about
    # avoiding deadlocks.
    cppclass timed_recursive_mutex:
        # may not be present, and we know nothing about it
        cppclass native_handle_type:
            pass

        void lock() except+
        bool try_lock()
        bool try_lock_for[T](const T& duration) except+
        bool try_lock_until[T](const T& time_point) except+
        void unlock()

        native_handle_type native_handle() except+


    # tags
    cppclass defer_lock_t:
        pass
    defer_lock_t defer_lock
    cppclass try_to_lock_t:
        pass
    try_to_lock_t try_to_lock
    cppclass adopt_lock_t:
        pass
    adopt_lock_t adopt_lock

    # lock_guard is probably unusable without cpp_locals because it
    # can't be default-constructed or moved.
    # We strongly recommend you do not use this while holding the GIL.
    # A safe way to construct with the GIL is to lock the mutex with py_safe_lock
    # and then construct the lock guard with adopt_lock_t.
    cppclass lock_guard[T]:
        ctypedef T mutex_type
        lock_guard(mutex_type&) except+
        lock_guard(mutex_type&, adopt_lock_t)

    # We strongly recommend that you do not use any blocking lock function with the GIL.
    # (try_lock is fine though).
    cppclass unique_lock[T]:
        ctypedef T mutex_type
        unique_lock()
        # This covers both the plain regular constructor, the 3 versions with tags
        # and two std::chrono constructors.  The two templated chrono versions stop
        # us from declaring the overloads explicitly.
        unique_lock(mutex_type&, ...) except+
        #unique_lock(mutex_type&, defer_lock_t)
        #unique_lock(mutex_type&, try_to_lock_t) except+
        ## this feels like it should be noexcet, but cppreference implies it isn't
        #unique_lock(mutex_type&, adopt_lock_t) except+

        void lock() except+
        bool try_lock() except+
        bool try_lock_for[T](const T& duration) except+
        bool try_lock_until[T](const T& time_point) except+
        void unlock() except+

        void swap(unique_lock& other)
        # "release" is definitely not the same as unlock. Noted here because
        # DW always makes this mistake and regrets it and wants to save you
        # from the same fate.
        mutex_type* release()

        mutex_type* mutex()
        bool owns_lock()
        bool operator bool()

    # scoped lock is probably unusable without cpp_locals.
    # It's also a variadic template type so can take potentially unlimited number of
    # arguments. Cython doesn't support this, so if you want more than 26 arguments,
    # you're on your own.
    # We strongly recommend that you do not use this while holding the GIL.
    # The safe way to use it when holding the GIL is to:
    # 1. first lock the mutexes with py_safe_lock
    # 2. construct scoped_lock with adopt_lock as the first argument
    cppclass scoped_lock[A=*, B=*, C=*, D=*, E=*, F=*, G=*, H=*, I=*, J=*, K=*,
                         L=*, M=*, N=*, O=*, P=*, Q=*, R=*, S=*, T=*, U=*, V=*,
                         W=*, X=*, Y=*, Z=*]:
        scoped_lock(...) except+

    cppclass once_flag:
        pass

    bool try_lock(...) except+
    # We strongly recommend that you do not call "lock" while holding the GIL.
    # See py_safe_lock.
    void lock(...) except+

    # We can't enforce this in the interface, but you need to make sure that Callable
    # doesn't require the GIL and doesn't throw Python exceptions.
    # You should also not call this with the GIL held.
    # We strong recommend using the py_safe_call_once wrappers below if you require the GIL.
    void call_once[Callable](once_flag&, Callable& callable,  ...) except +


cdef inline void _dummy_force_utility_code_inclusion() nogil:
    with nogil:
        pass

cdef extern from *:
    # Treat this as a different type from Cython's point of view so that if users use our safe functions
    # then they can't mix them with the unsafe functions using the same flag.
    cdef cppclass py_safe_once_flag "std::once_flag":
        pass

# Cython-specific wrappers to avoid deadlock.
# If you're holding the GIL then we strongly recommend you use these.
cdef extern from *:
    """
    #include <utility>
    #include <exception>

    namespace {

    static PyGILState_STATE __pyx_libcpp_mutex_limited_api_ensure_gil() {
        #if CYTHON_COMPILING_IN_LIMITED_API
        if ((__PYX_LIMITED_VERSION_HEX < 0x030d0000) && __Pyx_get_runtime_version() < 0x030d0000) {
            return PyGILState_Ensure();
        }
        #endif
        return PyGILState_LOCKED;  // Unused
    }

    #if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030d0000
    static void __pyx_libcpp_mutex_limited_api_release_gil(PyGILState_STATE gil_state) {
        if (__Pyx_get_runtime_version() < 0x030d0000)
            PyGILState_Release(gil_state);
    }
    #else
    #define __pyx_libcpp_mutex_limited_api_release_gil(ignore) (void)ignore
    #endif

    static int __pyx_libcpp_mutex_has_gil() {
        #if CYTHON_COMPILING_IN_LIMITED_API
            if ((__PYX_LIMITED_VERSION_HEX >= 0x030d0000) || __Pyx_get_runtime_version() >= 0x030d0000) {
                // In 3.13+ we can temporarily give up the GIL to find out what the thread state was
                PyThreadState *ts = PyThreadState_Swap(NULL);
                if (ts) {
                    PyThreadState_Swap(ts);
                    return 1;
                }
                return 0;
            }
            /* There is no way to know if we have the GIL. Therefore the only
            * thing we can safely do is make absolutely sure that we have it
            * in (__pyx_libcpp_mutex_limited_api_ensure_gil).
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

    template <typename F>
    class __pyx_libcpp_mutex_cleanup_on_exit {
        F on_exit;
        bool invoke = true;

    public:
        explicit __pyx_libcpp_mutex_cleanup_on_exit(F f)
            : on_exit(f)
        {}
        __pyx_libcpp_mutex_cleanup_on_exit(__pyx_libcpp_mutex_cleanup_on_exit &&rhs)
            : on_exit(std::move(rhs.on_exit))
            , invoke(rhs.invoke)
        {
            rhs.invoke = false;
        }

        __pyx_libcpp_mutex_cleanup_on_exit(const __pyx_libcpp_mutex_cleanup_on_exit&) = delete;
        __pyx_libcpp_mutex_cleanup_on_exit& operator=(const __pyx_libcpp_mutex_cleanup_on_exit&) = delete;

        ~__pyx_libcpp_mutex_cleanup_on_exit() {
            if (invoke)
                on_exit();
        }
    };

    template<typename F>
    __pyx_libcpp_mutex_cleanup_on_exit<F> __pyx_make_libcpp_mutex_cleanup_on_exit(F f) {
        return __pyx_libcpp_mutex_cleanup_on_exit<F>(std::move(f));
    }

    template <typename Callable, typename ... Args>
    void __pyx_cpp_py_safe_call_once(std::once_flag& flag, Callable& callable, Args&&... args) {
        class PyException : public std::exception {
        public:
            using std::exception::exception;
        };

        __Pyx_UnknownThreadState thread_state = __Pyx_SaveUnknownThread();
        auto on_exit = __pyx_make_libcpp_mutex_cleanup_on_exit(
            [&]() {
                __Pyx_RestoreUnknownThread(thread_state);
            });

        try {
            std::call_once(flag,
                [&](Args& ...args) {
                    // Make sure we have the GIL
                    PyGILState_STATE gil_state;
                    int had_gil_on_call = __Pyx_UnknownThreadStateDefinitelyHadGil(thread_state);
                    if (had_gil_on_call) {
                        __Pyx_RestoreUnknownThread(thread_state);
                    } else {
                        gil_state = PyGILState_Ensure();
                    }
                    auto on_callable_exit = __pyx_make_libcpp_mutex_cleanup_on_exit(
                        [&]() {
                            if (had_gil_on_call) {
                                thread_state = __Pyx_SaveUnknownThread();
                            } else {
                                PyGILState_Release(gil_state);
                            }
                        });

                    std::forward<Callable>(callable)(std::forward<Args>(args)...);
                    if (PyErr_Occurred()) {
                        throw PyException();
                    }
                },
                std::forward<Args>(args)...
            );
        } catch (const PyException&) {
            // Do nothing
        }
    }

    CYTHON_UNUSED void __pyx_cpp_py_safe_call_object_once(std::once_flag& flag, PyObject *py_callable) {
        auto callable = [py_callable]() {
            auto result = PyObject_CallObject(py_callable, nullptr);
            Py_XDECREF(result);
        };
        __pyx_cpp_py_safe_call_once(flag, callable);
    }

    template <typename LockableT>
    void __pyx_std_lock_wrapper(LockableT& arg0) {
        // std::lock only handles 2 or more arguments.
        // So create a 1 argument version.
        arg0.lock();
    }

    template <typename Lockable0T, typename Lockable1T, typename ... Lockables>
    void __pyx_std_lock_wrapper(Lockable0T& arg0, Lockable1T& arg1, Lockables&... args) {
        std::lock(arg0, arg1, args...);
    }

    inline void __pyx_libcpp_mutex_unlock() {} // no-op

    template <typename Lockable0T, typename ... Lockables>
    void __pyx_libcpp_mutex_unlock(Lockable0T& arg0, Lockables&... locks) {
        arg0.unlock();
        __pyx_libcpp_mutex_unlock(locks...);
    }

    template <typename LockableT>
    int __pyx_std_try_lock_wrapper(LockableT& arg0) {
        // std::try_lock only handles 2 or more arguments.
        // So create a 1 argument version.
        if (arg0.try_lock()) {
            return -1;
        }
        return 0;
    }

    template <typename Lockable0T, typename Lockable1T, typename ...Lockables>
    int __pyx_std_try_lock_wrapper(Lockable0T& arg0, Lockable1T& arg1, Lockables&... args) {
        return std::try_lock(arg0, arg1, args...);
    }

    template <typename... LockTs>
    void __pyx_py_safe_std_lock_release_lock_reacquire(LockTs& ...locks) {
        // Release the GIL, acquire the lock, then reacquire the GIL.
        // This is safe provided the user never holds the GIL while trying
        // to reacquire the lock (i.e. it's safe provided they always use
        // the py-safe wrappers).
        PyThreadState *_save;
        Py_UNBLOCK_THREADS
        try {
            __pyx_std_lock_wrapper(locks...);
        } catch (...) {
            // In this case, we probably can't reason about the state of the locks but we can at least
            // make sure the GIL is consistent.
            Py_BLOCK_THREADS
            throw;
        }
        Py_BLOCK_THREADS
        return;
    }

    template <typename... LockTs>
    void __pyx_py_safe_std_lock(LockTs& ...locks) {
        PyGILState_STATE gil_state = __pyx_libcpp_mutex_limited_api_ensure_gil();
        int had_gil_on_call = __pyx_libcpp_mutex_has_gil();
        auto on_exit = __pyx_make_libcpp_mutex_cleanup_on_exit(
            [&]() {
                __pyx_libcpp_mutex_limited_api_release_gil(gil_state);
            });
        if (!had_gil_on_call) {
            // Nothing special to do - just lock and quit
            __pyx_std_lock_wrapper(locks...);
            return;
        }

        // It's a real shame there's no try_lock for the GIL, otherwise
        // we could just defer this whole thing to c++ std::lock.
        if (__pyx_std_try_lock_wrapper(locks...) == -1) {
            // success!
            return;
        }
        __pyx_py_safe_std_lock_release_lock_reacquire(locks...);
    }

    template <typename MutexT>
    std::unique_lock<MutexT> __pyx_py_safe_construct_unique_lock(MutexT& mutex) {
        std::unique_lock<MutexT> l{mutex, std::defer_lock};
        __pyx_py_safe_std_lock(l);
        return l;
    }
    } // namespace
    """
    # Call a Python callable once, ensuring that the GIL is released before locking, and the callable
    # is called with the GIL held. The callable may be a Python object or C function.
    # If using a generic callable, the callable can throw either Python
    # or C++ exceptions. Both are treated like a C++ exception.
    # The GIL state on exit is the same as on entry.
    void py_safe_call_object_once "__pyx_cpp_py_safe_call_object_once" (py_safe_once_flag&, object callable) except +* nogil
    void py_safe_call_once "__pyx_cpp_py_safe_call_once" [Callable](py_safe_once_flag&, Callable callable, ...) except +* nogil

    # Call std::lock on the lockable objects with the GIL held, avoiding deadlock with the GIL.
    # (Unlike the standard library version, it works with a single argument).
    # The GIL state on exit is the same as on entry.
    void py_safe_lock "__pyx_py_safe_std_lock" (...) except+ nogil

    # construct a unique lock with the GIL held, avoiding deadlocks with the GIL
    unique_lock[MutexT] py_safe_construct_unique_lock "__pyx_py_safe_construct_unique_lock" [MutexT](MutexT& m) except+ nogil
