# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False
#
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
#
# Modified from Milk
# Copyright (C) 2010, Luis Pedro Coelho <lpc@cmu.edu>
# License: MIT

import numpy as np
import sys
from time import time

cimport numpy as np
cimport cython

#def myfunc(np.ndarray[np.float64_t, ndim=2] A):
def set_entropy(np.ndarray[np.int_t, ndim=1] data, N, clen):
    '''Set entropy function
    '''
    cdef Py_ssize_t j
    cdef int value
    cdef np.ndarray[np.int_t, ndim=1] counts = np.zeros((clen), dtype=np.int)
    for j from 0 <= j < N:
        value = data[j]
        if (value >= clen):
            raise RuntimeError("tree_model._tree.set_entropy: label value too large. aborting");
        counts[value] += 1.;

    '''
    Here is the formula:
    
    H = - \sum px \log(px)
      = - \sum (cx/N) \log( cx / N)
      = - 1/N \sum { cx [ \log cx - \log N ] }
      = - 1/N { (\sum cx \log cx ) - ( \sum cx \log N ) }
      = - 1/N { (\sum cx \log cx ) - N \log N }
      = ( - 1/N \sum cx \log cx ) + \log N
    '''
    
    cdef double entropy = 0.
    cdef double cx
    for j from 0 <= j < clen:
        cx = counts[j]
        if cx:
            entropy += cx * np.log(cx);
    entropy /= -N;
    entropy += np.log(N);
    return entropy
