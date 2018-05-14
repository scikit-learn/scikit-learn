# distutils: language = c++

# Cython interface file for wrapping the object
#
#

from libcpp.vector cimport vector

# c++ interface to cython
cdef extern from "src/libaode/aode_helper.h":
  cdef cppclass AODEHelper:
        AODEHelper(vector[vector[vector[int]]], vector[vector[vector[vector[vector[int]]]]], vector[int], int, int, int, int, int) except +
        vector[int] calculate_predictions(vector[vector[int]], int)
        vector[vector[double]] calculate_probabilities(vector[vector[int]], int)
        vector[vector[double]] calculate_naive_bayes(vector[vector[int]], int)

# cython wrapper class
cdef class AODE_helper:
    cdef AODEHelper *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self, vector[vector[vector[int]]] ft, vector[vector[vector[vector[vector[int]]]]] jft, vector[int] uniq, int na, int nc, int ns, int m, int m_lmt):
        self.thisptr = new AODEHelper(ft, jft, uniq, na, nc, ns, m, m_lmt)
    def __dealloc__(self):
        del self.thisptr
    def calculate_predictions(self, preds, num):
        return self.thisptr.calculate_predictions(preds, num)
    def calculate_probabilities(self, preds, num):
        return self.thisptr.calculate_probabilities(preds, num)


def get_tables(X, int num_samples, vector[int] y_indexes, int num_attr, x_indexes, vector[int] num_uniq_x_vals, int num_classes, vector[vector[vector[int]]] freq_table, vector[vector[vector[vector[vector[int]]]]] joint_freq_table):
    
    cdef int y, x, ii, xi, xj
    cdef int this_y_idx, this_class_idx

    for y in range(num_samples):
        this_y_idx = y_indexes[y]
        for x in range(num_attr):
            this_x_idx = x_indexes[x][X[y,x]]
            freq_table[this_y_idx][x][this_x_idx] += 1

    for ii in range(num_samples):
        this_class_idx = y_indexes[ii]
        for xj in range(num_attr):
            this_xj_idx = x_indexes[xj][X[ii,xj]]
            for xi in range(0,xj):
                this_xi_idx = x_indexes[xi][X[ii,xi]]
                joint_freq_table[this_class_idx][xj][xi][this_xj_idx][this_xi_idx] += 1

    return (freq_table, joint_freq_table)
