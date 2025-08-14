from libcpp cimport bool

cdef extern from "<functional>" namespace "std" nogil:
    cdef cppclass function[T]:
        function() except +
        function(T*) except +
        function(function&) except +
        function(void*) except +

        function operator=(T*)
        function operator=(function&)
        function operator=(void*)
        function operator=[U](U)

        bool operator bool()

    # Comparisons
    cdef cppclass greater[T=*]:
        # https://github.com/cython/cython/issues/3193
        greater() except +
        bool operator()(const T& lhs, const T& rhs) except +

    cdef cppclass reference_wrapper[T]:
        reference_wrapper()
        reference_wrapper(T)
        T& get() const
