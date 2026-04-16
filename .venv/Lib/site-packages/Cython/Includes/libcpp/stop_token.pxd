from libcpp cimport bool

cdef extern from * nogil:
    """
    // <thread> is defined in C++ and may define __cpp_lib_jthread
    // Including it in this way lets us cimport it in the C++11 condition_variable wrappings
    #include <thread>
    #if defined(__cpp_lib_jthread)
    #include <stop_token>
    #endif
    """

cdef extern from * namespace "std" nogil:
    cdef cppclass stop_token:
        bool stop_requested() noexcept
        bool stop_possible() noexcept

    cdef cppclass nostopstate_t:
        pass

    nostopstate_t nostopstate

    cdef cppclass stop_source:
        stop_source() except+
        stop_source(nostopstate_t) noexcept

        bool request_stop() noexcept
        void swap(stop_source& other) noexcept

        stop_token get_token() noexcept
        bool stop_requested() noexcept
        bool stop_possible() noexcept

    # stop_callback is not copyable or moveable which currently means it must
    # be heap-allocated in Cython (although use with cpp_locals should eventually be supported too)
    cdef cppclass stop_callback[Callback]:
        # in principle the second argument is a template argument with a "std::constructable_from" constraint,
        # but for Cython's purposes it probably makes sense to assume no conversion
        stop_callback(stop_token st, Callback cb) noexcept

cdef extern from *:
    """
    #if defined(__cpp_lib_jthread)
    #include <optional>
    #include <utility>

    namespace {
        using __pyx_func_ptr_stop_callback = std::stop_callback<void (*)()>;

        class __pyx_python_stop_callback_holder {
            class callable_py_object_holder {
                PyObject *o;

              public:
                explicit callable_py_object_holder(PyObject *o)
                    : o(o)
                {
                    Py_INCREF(o);
                }

                callable_py_object_holder(callable_py_object_holder&& rhs)
                    : o(std::exchange(rhs.o, nullptr))
                {
                }

                ~callable_py_object_holder() {
                    if (o) {
                        PyGILState_STATE state = PyGILState_Ensure();
                        Py_DECREF(o);
                        PyGILState_Release(state);
                    }
                }

                callable_py_object_holder(const callable_py_object_holder&) = delete;
                callable_py_object_holder& operator=(const callable_py_object_holder&) = delete;

                void operator()() const {
                    PyGILState_STATE state = PyGILState_Ensure();
                    PyObject *result = PyObject_CallObject(o, NULL);
                    if (!result) {
                        PyObject *s = PyUnicode_FromString("python_stop_callback_holder callback");
                        PyErr_WriteUnraisable(s);
                        Py_XDECREF(s);
                    } else {
                        Py_DECREF(result);
                    }
                    PyGILState_Release(state);
                }
            };

            std::optional<std::stop_callback<callable_py_object_holder>> callback;

          public:
            __pyx_python_stop_callback_holder() = default;
            __pyx_python_stop_callback_holder(std::stop_token token, PyObject *callable) {
                initialize(std::move(token), callable);
            }
            __pyx_python_stop_callback_holder(const __pyx_python_stop_callback_holder&) = delete;

            __pyx_python_stop_callback_holder& operator=(const __pyx_python_stop_callback_holder&) = delete;

            void initialize(std::stop_token token, PyObject *callable) {
                callback.emplace(std::move(token), callable_py_object_holder(callable));
            }
        };
    }
    #endif // defined(__cpp_lib_jthread)
    """
    # This is provided as a convenience mainly as a reminder to use nogil functions!
    ctypedef stop_callback[void (*)() noexcept nogil] func_ptr_stop_callback "__pyx_func_ptr_stop_callback"

    # A fairly thin wrapper to let you create a stop callback with a Python object.
    # For most uses, it should be created empty and then filled with "initialize"
    cdef cppclass python_stop_callback_holder "__pyx_python_stop_callback_holder":
        python_stop_callback_holder()
        python_stop_callback_holder(stop_token token, object callable)
        void initialize(stop_token token, object callable)
