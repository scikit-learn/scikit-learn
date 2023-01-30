cdef class BaseCriterion:
    def __cinit__(self):
        self.A = 0.

cdef class Criterion(BaseCriterion):
    cdef int init(self) except -1:
        self.A = 1.


cdef class BaseSplitter:
    def __cinit__(self, BaseCriterion criterion, double max_features):
        self.criterion = criterion
        self.max_features = max_features

cdef class Splitter(BaseSplitter):
    def __cinit__(self, BaseCriterion criterion, double max_features):
        self.criterion = criterion
        self.max_features = max_features

    cdef int init(
        self,
        object X
    ) except -1:
        self.X = X
        <Criterion>(self.criterion).init()
        return 0

    cpdef int test(self, object X):
        return self.init(X)