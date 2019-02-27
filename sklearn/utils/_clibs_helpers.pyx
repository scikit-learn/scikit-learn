cimport openmp
from cython.parallel import prange

from ._joblib import effective_n_jobs


cpdef int effective_n_jobs_openmp(n_jobs=None):
    """Determine the number of threads for OpenMP

    The behavior differs from joblib's effective_n_jobs for n_jobs=None in
    which case effective_n_jobs_openmp returns openmp.omp_get_max_threads().
    """
    if n_jobs is None:
        return openmp.omp_get_max_threads()

    return effective_n_jobs(n_jobs)


def check_num_threads(int N=100):
    """Return the default number of threads used for parallel calls to openMP
    """
    cdef double Ysum = 0
    cdef int i, num_threads

    for i in prange(N, nogil=True):
        num_threads = openmp.omp_get_num_threads()
        Ysum += i

    # Assert that the parallel sum result is correct
    assert Ysum == N * (N - 1) / 2

    return num_threads
