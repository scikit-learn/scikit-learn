IF SKLEARN_OPENMP_PARALLELISM_ENABLED:
    import os
    cimport openmp
    from joblib import cpu_count


def _openmp_parallelism_enabled():
    """Determines whether scikit-learn has been built with OpenMP
    
    It allows to retrieve at runtime the information gathered at compile time.
    """
    # SKLEARN_OPENMP_PARALLELISM_ENABLED is resolved at compile time during
    # cythonization. It is defined via the `compile_time_env` kwarg of the
    # `cythonize` call and behaves like the `-D` option of the C preprocessor.
    return SKLEARN_OPENMP_PARALLELISM_ENABLED


cpdef _openmp_effective_n_threads(n_threads=None):
    """Determine the effective number of threads to be used for OpenMP calls

    - For ``n_threads = None``,
      - if the ``OMP_NUM_THREADS`` environment variable is set, return
        ``openmp.omp_get_max_threads()``
      - otherwise, return the minimum between ``openmp.omp_get_max_threads()``
        and the number of cpus, taking cgroups quotas into account. Cgroups 
        quotas can typically be set by tools such as Docker.
      The result of ``omp_get_max_threads`` can be influenced by environment
      variable ``OMP_NUM_THREADS`` or at runtime by ``omp_set_num_threads``.

    - For ``n_threads > 0``, return this as the maximal number of threads for
      parallel OpenMP calls.

    - For ``n_threads < 0``, return the maximal number of threads minus
      ``|n_threads + 1|``. In particular ``n_threads = -1`` will use as many
      threads as there are available cores on the machine.

    - Raise a ValueError for ``n_threads = 0``.

    If scikit-learn is built without OpenMP support, always return 1.
    """
    if n_threads == 0:
        raise ValueError("n_threads = 0 is invalid")

    IF SKLEARN_OPENMP_PARALLELISM_ENABLED:
        if os.getenv("OMP_NUM_THREADS"):
            # Fall back to user provided number of threads making it possible
            # to exceed the number of cpus.
            max_n_threads = openmp.omp_get_max_threads()
        else:
            max_n_threads = min(openmp.omp_get_max_threads(), cpu_count())

        if n_threads is None:
            return max_n_threads
        elif n_threads < 0:
            return max(1, max_n_threads + n_threads + 1)

        return n_threads
    ELSE:
        # OpenMP disabled at build-time => sequential mode
        return 1

    
