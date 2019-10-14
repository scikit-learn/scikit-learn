IF SKLEARN_OPENMP_SUPPORTED:
    cimport openmp
    from . import cpu_count


cpdef _openmp_effective_n_threads(n_threads=None):
    """Determine the effective number of threads used for parallel OpenMP calls

    - For ``n_threads = None``, returns the minimum between the number of cpus
      and ``openmp.omp_get_max_threads()``.
      The result of ``omp_get_max_threads`` can be influenced by environment
      variable ``OMP_NUM_THREADS`` or at runtime by ``omp_set_num_threads``.
    - For ``n_threads > 0``, use this as the maximal number of threads for
      parallel OpenMP calls.
    - For ``n_threads < 0``, use the maximal number of threads minus
      ``|n_threads + 1|``. In particular ``n_threads = -1`` will use as many
      threads as there are available cores on the machine.
    - Raise a ValueError for ``n_threads = 0``.

    If scikit-learn is built without OpenMP support, always return 1.
    """
    if n_threads == 0:
        raise ValueError("n_threads = 0 is invalid")

    IF SKLEARN_OPENMP_SUPPORTED:
        max_n_threads = min(openmp.omp_get_max_threads(), cpu_count())

        if n_threads is None:
            return max_n_threads
        elif n_threads < 0:
            return max(1, max_n_threads + n_threads + 1)

        return n_threads
    ELSE:
        # OpenMP not supported => sequential mode
        return 1

    
