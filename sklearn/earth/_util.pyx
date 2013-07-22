# distutils: language = c
# cython: cdivision = True
# cython: boundscheck = False
# cython: wraparound = False
# cython: profile = False

import numpy as np
from libc.math cimport sqrt, log

cdef FLOAT_t log2(FLOAT_t x):
    return log(x) / log(2.0)

cpdef apply_weights_2d(cnp.ndarray[FLOAT_t, ndim=2] B, cnp.ndarray[FLOAT_t, ndim=1] weights):
    cdef INDEX_t i
    cdef INDEX_t j
    cdef INDEX_t m = B.shape[0]
    cdef INDEX_t n = B.shape[1]
    for i in range(m):
        for j in range(n):
            B[i,j] *= sqrt(weights[i])
        
cpdef apply_weights_1d(cnp.ndarray[FLOAT_t, ndim=1] y, cnp.ndarray[FLOAT_t, ndim=1] weights):
    cdef INDEX_t i
    cdef INDEX_t m = y.shape[0]
    for i in range(m):
        y[i] *= sqrt(weights[i])

cpdef FLOAT_t gcv(FLOAT_t mse, INDEX_t basis_size, INDEX_t data_size, FLOAT_t penalty):
    return mse * gcv_adjust(basis_size, data_size, penalty)

cpdef FLOAT_t gcv_adjust(INDEX_t basis_size, INDEX_t data_size, FLOAT_t penalty):
    return 1.0 / ((1 - ((basis_size + penalty*(basis_size - 1))/data_size)) ** 2)

cpdef str_pad(string, length):
    if len(string) >= length:
        return string[0:length]
    pad = length - len(string)
    return string + ' '*pad

cpdef ascii_table(header, data):
    '''
    header - list of strings representing the header row
    data - list of lists of strings representing data rows
    '''
    m = len(data)
    n = len(header)
    column_widths = [len(head) for head in header]
    for i, row in enumerate(data):
        for j, col in enumerate(row):
            if len(col) > column_widths[j]:
                column_widths[j] = len(col)
    
    for j in range(n):
        column_widths[j] += 1
    
    result = ''
    for j, col_width in enumerate(column_widths):
        result += '-'*col_width + '-'
    result += '\n'
    for j, head in enumerate(header):
        result += str_pad(head,column_widths[j]) + ' '
    result += '\n'
    for j, col_width in enumerate(column_widths):
        result += '-'*col_width + '-'
    for i, row in enumerate(data):
        result += '\n'
        for j, item in enumerate(row):
            result += str_pad(item,column_widths[j]) + ' '
    result += '\n'
    for j, col_width in enumerate(column_widths):
        result += '-'*col_width + '-'
    return result
    
    
    
    
    
    
