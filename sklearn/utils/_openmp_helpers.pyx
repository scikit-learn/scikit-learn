IF SKLEARN_OPENMP_SUPPORTED:
    cimport openmp
    from joblib import effective_n_jobs


cpdef _openmp_effective_n_threads(n_threads=None):
    """Determine the effective number of threads used for parallel OpenMP calls

    - For ``n_threads = None``, returns the minimum between
      openmp.omp_get_max_threads() and joblib.effective_n_jobs(-1).
      The result of ``omp_get_max_threads`` can be influenced by environement
      variable ``OMP_NUM_THREADS`` or at runtime by ``omp_set_num_threads``.
    - For ``n_threads > 0``, use this as the maximal number of threads for
      parallel OpenMP calls.
    - For ``n_threads < 0``, use the maximal number of threads minus
      ``|n_threads + 1|``.
    - Raise a ValueError for ``n_threads = 0``.

    If scikit-learn is built without OpenMP support, always return 1.
    """
    if n_threads == 0:
        raise ValueError("n_threads = 0 is invalid")

    IF SKLEARN_OPENMP_SUPPORTED:
        max_n_threads = min(openmp.omp_get_max_threads(), effective_n_jobs(-1))

        if n_threads is None:
            return max_n_threads
        elif n_threads < 0:
            return max(1, max_n_threads + n_threads + 1)

        return n_threads
    ELSE:
        # OpenMP not supported => sequential mode
        return 1

    
