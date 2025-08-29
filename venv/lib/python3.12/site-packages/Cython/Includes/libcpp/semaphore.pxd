from libcpp cimport bool

cdef extern from "<semaphore>" namespace "std" nogil:
    # Notes on templating:
    # Cython doesn't currently support non-class template types so it's
    # declared here as a class template.
    # counting_semaphore has a default max value anyway, so you typically
    # don't need to use the templates. If you do though, then there's a
    # few tricks - e.g. you can define a class with a cname of "3" for
    # example.

    # Note on thread-safety:
    # You should use these classes without the GIL. It's very easy
    # to deadlock with the GIL, and we're assuming that anyone going
    # low enough level to be using semaphores doesn't want the overhead
    # of Cython handling GIL safety.

    cdef cppclass counting_semaphore[LeastMaxValue=*]:
        counting_semaphore(ptrdiff_t desired)

        void release() except+
        void release(ptrdiff_t update) except+
        void acquire() except+
        bool try_acquire()

        # unfortunately, we don't currently define chrono types
        bool try_acquire_for[T](T duration) except+
        bool try_acquire_until[T](T timepoint) except+

        ptrdiff_t max()

    cdef cppclass binary_semaphore:
        binary_semaphore(ptrdiff_t desired)

        void release() except+
        void release(ptrdiff_t update) except+
        void acquire() except+
        bool try_acquire()

        # unfortunately, we don't currently define chrono types
        bool try_acquire_for[T](T timepoint) except+
        bool try_acquire_until[T](T duration) except+

        ptrdiff_t max()
