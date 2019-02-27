cimport openmp
from cython.parallel import prange


cpdef int openmp_n_threads(n_jobs=None):
    """Determine the number of threads for OpenMP

    The behavior differs from joblib's effective_n_jobs for n_jobs=None in
    which case openmp_n_threads returns openmp.omp_get_max_threads().
    """
    if n_jobs is None:
        return openmp.omp_get_max_threads()
    elif n_jobs == 0:
        raise ValueError("n_jobs=0 is invalid")
    elif n_jobs < 0:
        return openmp.omp_get_max_threads() + n_jobs + 1
    else:
        return n_jobs


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
