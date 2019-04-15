# the build will fail if OpenMP support is not handled correctly.

IF SKLEARN_OPENMP_SUPPORTED:
    cimport openmp


def check_openmp():
    IF SKLEARN_OPENMP_SUPPORTED:
        openmp.omp_get_max_threads()
