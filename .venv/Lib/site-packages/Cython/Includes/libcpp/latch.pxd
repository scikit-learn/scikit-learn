from libcpp cimport bool

cdef extern from "<latch>" namespace "std" nogil:
    # Note on thread safety:
    # For any of the blocking functions, you should be very careful
    # to ensure that you are not deadlocked on the GIL.
    cdef cppclass latch:
        latch(ptrdiff_t expected)

        void count_down() except+
        void count_down(ptrdiff_t n) except+
        bool try_wait()
        void wait() except+
        void arrive_and_wait() except+
        void arrive_and_wait(ptrdiff_t n) except+

        ptrdiff_t max()
