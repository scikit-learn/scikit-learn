cdef LearningRate get_learning_rate(int learning_rate)

cdef class LearningRate(object):
    cdef double step(self, double eta0, double alpha, int t, int power_t)

cdef class Constant(LearningRate):
    cdef double step(self, double eta0, double alpha, int t, int power_t)

cdef class Optimal(LearningRate):
    cdef double step(self, double eta0, double alpha, int t, int power_t)

cdef class InvScaling(LearningRate):
    cdef double step(self, double eta0, double alpha, int t, int power_t)
