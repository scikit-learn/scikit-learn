cdef extern from "<barrier>" namespace "std" nogil:
    # Note on thread safety:
    # For any of the blocking functions here you should be very
    # careful to ensure that you are not deadlocked on the GIL.
    cdef cppclass barrier[CompletionFunction = *]:
        cppclass arrival_token:
            pass

        barrier(ptrdiff_t expected) except+
        # We *really* advise that CompletionFunction is nogil (but because it's a template, can't enforce it)
        barrier(ptrdiff_t expected, CompletionFunction f) except+

        arrival_token arrive() except+
        arrival_token arrive(ptrdiff_t n) except+

        void wait(arrival_token t) except+

        void arrive_and_wait() except+

        void arrive_and_drop() except+

        ptrdiff_t max()
