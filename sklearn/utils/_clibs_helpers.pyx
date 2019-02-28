cimport openmp
from cython.parallel import prange


cpdef int openmp_num_threads(n_threads=None):
    """Determine the effective number of threads used for parallel OpenMP calls

    For n_threads=None, returns openmp.omp_get_max_threads().
    For n_threads > 0, use this as the maximal number of threads for parallel
    OpenMP calls and for n_threads < 0, use the maximal number of threads minus
    - n_threads + 1.
    Raise a ValueError for n_threads=0.
    """
    if n_threads == 0:
        raise ValueError("n_threads=0 is invalid")
    elif n_threads < 0:
        return max(1, openmp.omp_get_max_threads() + n_threads + 1)
    else:
        return n_threads


def check_num_threads(int n=100):
    """Return the default number of threads used for parallel calls to openMP
    """
    cdef double n_sum = 0
    cdef int i, num_threads

    for i in prange(n, nogil=True):
        num_threads = openmp.omp_get_num_threads()
        n_sum += i

    # Assert that the parallel sum result is correct
    assert n_sum == n * (n - 1) / 2

    return num_threads
