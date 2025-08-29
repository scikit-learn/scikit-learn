from libcpp cimport bool
from libcpp.exception cimport exception_ptr

cdef extern from "<future>":
    pass  # include first

cdef extern from *:
    """
    static void __Pyx_CppExn2PyErr();
    CYTHON_UNUSED static void __pyx_future_error_handler() {
        try {
            throw;
        } catch (const std::future_error &e) {
            PyObject *args = Py_BuildValue(
                "sis",
                "std::future_error",
                e.code().value(),
                e.what()
            );
            if (!args) return;
            Py_INCREF(PyExc_RuntimeError);
            // It'd be nice to raise something more specific, but I can't easily do that
            // so I've tried to at least make them identifiable by the arguments.
            PyErr_Restore(PyExc_RuntimeError, args, NULL);
        } catch (...) {
            // Forward to the normal Cython exception conversion;
            __Pyx_CppExn2PyErr();
        }
    }
    """
    void __pyx_future_error_handler() except+

cdef extern from "<future>" namespace "std" nogil:
    cdef enum class future_status:
        ready
        timeout
        deferred

    cdef cppclass future[T]:
        future() noexcept
        future(future& other) noexcept  # slightly mistyped move constructor

        future& operator=(future& other) noexcept  # slightly mistyped move assignment

        shared_future[T] share() noexcept

        # Not strictly right when T is a reference or void
        T get() except +__pyx_future_error_handler

        bool valid() noexcept

        void wait() except+
        # For wait_for and wait_until we don't have chrono yet
        future_status wait_for(...) except+
        future_status wait_until(...) except+

    cdef cppclass shared_future[T]:
        shared_future() noexcept
        shared_future(const shared_future& other) noexcept  # noexcept only from c++17 strictly
        shared_future(future[T]& other) noexcept  # slightly mistyped move

        shared_future& operator=(const shared_future& other) noexcept  # noexcept only from c++17 strictly
        shared_future& operator=(future& other) noexcept  # slightly mistyped move

        # Not strictly right when T is a reference or void
        const T& get() except +__pyx_future_error_handler

        bool valid() noexcept

        void wait() except+
        # For wait_for and wait_until we don't have chrono yet
        future_status wait_for(...) except+
        future_status wait_until(...) except+

    cdef cppclass promise[T]:
        promise() noexcept
        # don't expose the allocator version
        promise(promise& other) noexcept  # slightly mistyped move

        promise& operator=(promise& other) noexcept  # slightly mistyped move

        void swap(promise& other) noexcept

        future[T] get_future() except+

        # Note that the set_* functions all have specialized error handling.
        # A std::future_error will be converted to a RuntimeError with 3 arguments.
        # The first argument is the string "std::future_error",
        # the second is the error code as an int for comparison with future_errc,
        # and the third is the message.
        void set_value(T value) except +__pyx_future_error_handler
        void set_value() except +__pyx_future_error_handler  # only for void specialization
        void set_value_at_thread_exit(T value) except +__pyx_future_error_handler
        void set_value_at_thread_exit() except +__pyx_future_error_handler  # only for void specialization

        void set_exception(exception_ptr p) except +__pyx_future_error_handler
        void set_exception_at_thread_exit(exception_ptr p) except +__pyx_future_error_handler

    cdef enum class future_errc:
        broken_promise
        future_already_retrieved
        promise_already_satisfied
        no_state
