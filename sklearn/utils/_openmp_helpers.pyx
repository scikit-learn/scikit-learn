import os
from joblib import cpu_count
from cython.parallel import prange
from ctypes import addressof

from.parallel import _get_threadpool_controller


# Module level cache for cpu_count as we do not expect this to change during
# the lifecycle of a Python program. This dictionary is keyed by
# only_physical_cores.
_CPU_COUNTS = {}


def _openmp_parallelism_enabled():
    """Determines whether scikit-learn has been built with OpenMP

    It allows to retrieve at runtime the information gathered at compile time.
    """
    # SKLEARN_OPENMP_PARALLELISM_ENABLED is resolved at compile time and defined
    # in _openmp_helpers.pxd as a boolean. This function exposes it to Python.
    return SKLEARN_OPENMP_PARALLELISM_ENABLED


cpdef _openmp_effective_n_threads(n_threads=None, only_physical_cores=True):
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

    Passing the `only_physical_cores=False` flag makes it possible to use extra
    threads for SMT/HyperThreading logical cores. It has been empirically
    observed that using as many threads as available SMT cores can slightly
    improve the performance in some cases, but can severely degrade
    performance other times. Therefore it is recommended to use
    `only_physical_cores=True` unless an empirical study has been conducted to
    assess the impact of SMT on a case-by-case basis (using various input data
    shapes, in particular small data shapes).

    If scikit-learn is built without OpenMP support, always return 1.
    """
    if n_threads == 0:
        raise ValueError("n_threads = 0 is invalid")

    if not SKLEARN_OPENMP_PARALLELISM_ENABLED:
        # OpenMP disabled at build-time => sequential mode
        return 1

    if os.getenv("OMP_NUM_THREADS"):
        # Fall back to user provided number of threads making it possible
        # to exceed the number of cpus.
        max_n_threads = omp_get_max_threads()
    else:
        try:
            n_cpus = _CPU_COUNTS[only_physical_cores]
        except KeyError:
            n_cpus = cpu_count(only_physical_cores=only_physical_cores)
            _CPU_COUNTS[only_physical_cores] = n_cpus
        max_n_threads = min(omp_get_max_threads(), n_cpus)

    if n_threads is None:
        return max_n_threads
    elif n_threads < 0:
        return max(1, max_n_threads + n_threads + 1)

    return n_threads


ctypedef void (*openblas_dojob_callback)(int, void*, int) noexcept nogil
ctypedef void (*openblas_threads_callback)(int, openblas_dojob_callback, int, size_t, void*, int)
ctypedef void (*openblas_set_threads_callback_function_type)(openblas_threads_callback)


# Callback for OpenBLAS to make it use an OpenMP backend, took from
# https://github.com/OpenMathLib/OpenBLAS/pull/4577#issue-2204960832
cdef void openblas_openmp_callback(
    int sync,
    openblas_dojob_callback dojob,
    int numjobs,
    size_t jobdata_elsize,
    void *jobdata,
    int dojob_data
) noexcept nogil:
    cdef int i
    cdef void *element_adrr

    for i in prange(numjobs, nogil=True):
        element_adrr = <void *>(((<char *>jobdata) + (<unsigned>i) * jobdata_elsize))
        dojob(i, element_adrr, dojob_data)


def set_openblas_openmp_callback():
    controller = _get_threadpool_controller()
    openblas_controllers = controller.select(internal_api="openblas").lib_controllers

    cdef openblas_set_threads_callback_function_type f_ptr 

    for ct in openblas_controllers:
        lib = ct.dynlib

        # openblas_set_threads_callback_function is available since v0.3.28
        if hasattr(lib, "openblas_set_threads_callback_function"):
            func = lib.openblas_set_threads_callback_function

            # cast to the correct function pointer type
            f_ptr = (<openblas_set_threads_callback_function_type*><size_t>addressof(func))[0]
            f_ptr(openblas_openmp_callback)

            print("OpenBLAS OpenMP backend callback set")
