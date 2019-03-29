# cython: language_level=3

def check_build():
    return

import os 

if os.getenv('SKLEARN_OPENMP'):
    cimport openmp

def tst():
    if os.getenv('SKLEARN_OPENMP'):
        openmp.omp_get_max_threads()
