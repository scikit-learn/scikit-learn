# encoding: utf-8
# cython: profile=True
from collections import defaultdict

import numpy as np

cimport numpy as np
cimport cython

np.import_array()

from libc.stdlib cimport malloc, free


def find_pivot(set nodes, set new_candidates, dict neighbors, long lim):
    cdef long pivot = -1
    cdef long pivot_degree = -1
    cdef long n
    cdef long degree
    for n in nodes:
        degree = len(new_candidates & neighbors[n])
        if degree > pivot_degree:
            pivot = n
            pivot_degree = degree
        if degree == lim:
            break
    return pivot, pivot_degree
