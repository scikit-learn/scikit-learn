cdef LearningRate get_learning_rate(str learning_rate):
    if learning_rate == "constant":
        return Constant()

cdef class LearningRate(object):
    cdef double step(self, double eta0):
        pass

cdef class Constant(LearningRate):
    cdef double step(self, double eta0):
        return eta0
