cdef LearningRate get_learning_rate(str learning_rate)

cdef class LearningRate(object):
    cdef double step(self, double eta0)

cdef class Constant(LearningRate):
    cdef double step(self, double eta0)
