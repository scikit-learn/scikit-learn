# cython: language_level=3

# the build will fail if OpenMP support is not handled correctly.

IF SKLEARN_OPENMP:
    cimport openmp


def check_openmp():
    IF SKLEARN_OPENMP:
        openmp.omp_get_max_threads()
