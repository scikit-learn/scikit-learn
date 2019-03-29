# cython: language_level=3

def check_build():
    return

cimport openmp

def tst():
    openmp.omp_get_max_threads()
