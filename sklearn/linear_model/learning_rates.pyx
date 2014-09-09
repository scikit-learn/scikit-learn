
cdef extern from "stdlib.h":
    int strcmp(const char*, const char*)

cdef LearningRate get_learning_rate(char* learning_rate):
    if(strcmp(learning_rate, "constant")):
        return Constant()

cdef class LearningRate:
    cdef double step(self, double eta0):
        pass

cdef class Constant(LearningRate):
    cdef double step(self, double eta0):
        return eta0
