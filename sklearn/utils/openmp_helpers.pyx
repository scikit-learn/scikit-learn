cimport openmp
from joblib import effective_n_jobs


cpdef _openmp_effective_n_threads(n_threads=None):
    """Determine the effective number of threads used for parallel OpenMP calls

    - For ``n_threads = None``, returns the minimum between
      openmp.omp_get_max_threads() and joblib.effective_n_jobs(-1).
    - For ``n_threads > 0``, use this as the maximal number of threads for
      parallel OpenMP calls.
    - For ``n_threads < 0``, use the maximal number of threads minus
      ``|n_threads + 1|``.
    - Raise a ValueError for ``n_threads = 0``.
    """
    if n_threads == 0:
        raise ValueError("n_threads = 0 is invalid")

    max_threads = min(openmp.omp_get_max_threads(), effective_n_jobs(-1))

    if n_threads is None:
        return max_threads
    elif n_threads < 0:
        return max(1, max_threads + n_threads + 1)

    return n_threads
