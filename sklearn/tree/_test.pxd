cdef class BaseCriterion:
    """Abstract interface for criterion."""
    cdef double A 

cdef class Criterion(BaseCriterion):
    cdef int init(self) except -1

cdef class BaseSplitter:
    """Abstract interface for splitter."""

    # Internal structures
    cdef public BaseCriterion criterion  # Impurity criterion
    cdef public double max_features      # Number of features to test

cdef class Splitter(BaseSplitter):
    cdef const double[:, ::1] X

    cdef int init(
        self,
        object X
    ) except -1

    # used to test the program in Python
    cpdef int test(self, object X)