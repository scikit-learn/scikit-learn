import os
from time import perf_counter

from cython.parallel import prange
from joblib import cpu_count


# Module level cache for cpu_count as we do not expect this to change during
# the lifecycle of a Python program. This dictionary is keyed by
# only_physical_cores.
_CPU_COUNTS = {}

# Fallback amount of work (in ~simple-instruction units) that must be
# available per thread before bothering to parallelize a loop with OpenMP;
# see `_use_threads_for_workload` in `_openmp_helpers.pxd`. Used until
# `_calibrate_min_instructions_per_thread` has run once, and whenever it
# fails to produce a sensible measurement.
_DEFAULT_MIN_INSTRUCTIONS_PER_THREAD = 2000

# None means "not calibrated yet". Calibration always ends up setting this
# to an int (falling back to the default above on any issue), so it only
# ever runs once per process.
_MIN_INSTRUCTIONS_PER_THREAD = None


def _min_instructions_per_thread():
    """Minimum amount of work (in ~simple-instruction units) that must be
    available per thread before bothering to parallelize a loop with
    OpenMP; see ``_use_threads_for_workload`` in ``_openmp_helpers.pxd``.

    The actual cost of dispatching an OpenMP parallel region varies a lot
    across machines and OpenMP runtimes/configurations (e.g. whether idle
    threads spin or sleep between parallel regions). Rather than guess at
    it, this is measured empirically, once per process, by
    ``_calibrate_min_instructions_per_thread`` the first time this function
    is called; every call after that just returns the cached result.
    """
    global _MIN_INSTRUCTIONS_PER_THREAD
    if _MIN_INSTRUCTIONS_PER_THREAD is None:
        _MIN_INSTRUCTIONS_PER_THREAD = _calibrate_min_instructions_per_thread()
    return _MIN_INSTRUCTIONS_PER_THREAD


def _bench_prange_dispatch(int n_threads, Py_ssize_t repeats, Py_ssize_t idle_gap):
    """Time ``repeats`` dispatches of a near-empty OpenMP parallel region.

    Each dispatch parallelizes a loop with exactly one iteration per thread
    (as close to "no actual work" as an OpenMP parallel region gets), so the
    measured time is essentially the fixed cost of spinning up/joining the
    thread team. ``idle_gap`` trivial sequential iterations are inserted
    between dispatches, giving worker threads a short idle period before the
    next one, similar to the gaps between the parallel regions this is
    calibrating for.
    """
    cdef:
        Py_ssize_t _r, _i, j
        Py_ssize_t _sink = 0
    t0 = perf_counter()
    with nogil:
        for _r in range(repeats):
            for _i in prange(n_threads, schedule='static', num_threads=n_threads):
                pass
            for j in range(idle_gap):
                _sink = j
    return perf_counter() - t0


def _calibrate_min_instructions_per_thread():
    """Measure the actual cost of dispatching an OpenMP parallel region on
    this machine, and convert it to the "simple instructions" unit that
    ``_use_threads_for_workload`` uses.

    ``n_threads``, ``repeats`` and ``idle_gap`` below were empirically found
    to give stable, representative measurements: half the physical cores is
    enough to exercise real thread dispatch/synchronization overhead without
    over-subscribing on small machines; more repeats or a larger idle gap
    did not noticeably change the result.

    Falls back to ``_DEFAULT_MIN_INSTRUCTIONS_PER_THREAD`` if OpenMP is
    disabled or the measurement fails or looks unreasonable for any reason
    (this must never raise, since it runs implicitly on first use).
    """
    if not SKLEARN_OPENMP_PARALLELISM_ENABLED:
        return _DEFAULT_MIN_INSTRUCTIONS_PER_THREAD

    try:
        try:
            n_physical_cores = _CPU_COUNTS[True]
        except KeyError:
            n_physical_cores = cpu_count(only_physical_cores=True)
            _CPU_COUNTS[True] = n_physical_cores
        n_threads = max(2, n_physical_cores // 2)
        repeats = 100
        idle_gap = 10

        # The very first OpenMP parallel region dispatched in a process pays
        # for one-time thread pool creation (tens of us), on top of the
        # steady-state dispatch cost this is trying to measure. Warm it up,
        # untimed, before measuring.
        _bench_prange_dispatch(n_threads, repeats, idle_gap)

        # Take several independent measurements and keep the minimum:
        # external interference (CPU frequency scaling, other processes
        # scheduler preemption, ...) can only ever slow a measurement down,
        # never speed it up below the true steady-state dispatch cost,
        # so the minimum across samples is the one least
        # contaminated by such noise (the same reasoning `timeit` uses).
        best_dispatch_seconds = None
        for _ in range(7):
            elapsed = _bench_prange_dispatch(n_threads, repeats, idle_gap)
            dispatch_seconds = elapsed / repeats
            if best_dispatch_seconds is None or dispatch_seconds < best_dispatch_seconds:
                best_dispatch_seconds = dispatch_seconds

        # The default (2000) was picked as the rough equivalent
        # of ~1us of simple work, i.e. ~2e9 simple-instructions-per-second;
        # convert the measured dispatch overhead to that unit:
        measured = round(best_dispatch_seconds * 2e9)
        return max(1, measured)
    except Exception:
        return _DEFAULT_MIN_INSTRUCTIONS_PER_THREAD


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
