import os
from time import perf_counter

from cython.parallel import prange
from joblib import cpu_count


# Module level cache for cpu_count as we do not expect this to change during
# the lifecycle of a Python program. This dictionary is keyed by
# only_physical_cores.
_CPU_COUNTS = {}


# Cache of calibrated per-thread dispatch-overhead values (in seconds),
# keyed by thread-count "bucket" (see `_omp_prange_dispatch_overhead`).
# Populated lazily: calibration always ends up setting a bucket's entry to
# a float (falling back to the default below on any issue), so a given
# bucket only ever gets calibrated once per process.
_OVERHEAD_PER_THREAD_CACHE = {}

# Fallback per-thread OpenMP dispatch overhead (in seconds), used until
# `_calibrate_overhead_per_thread` has run once for a given bucket, and
# whenever it fails to produce a sensible measurement.
_DEFAULT_OVERHEAD_PER_THREAD = 1e-6


def _omp_prange_dispatch_overhead(n_threads):
    """Estimated total dispatch overhead (in seconds) of an OpenMP parallel
    region run with a team of ``n_threads`` threads on this machine; used by
    ``_use_threads_for_workload`` below to decide whether parallelizing is
    worth it.

    The actual cost of dispatching an OpenMP parallel region varies a lot
    across machines and OpenMP runtimes/configurations (e.g. whether idle
    threads spin or sleep between parallel regions). Rather than guess at
    it, this is measured empirically by ``_calibrate_overhead_per_thread``,
    the first time a given team size is requested; every call after that
    just scales the cached per-bucket value by ``n_threads``.

    Dispatch overhead does not scale linearly with the size of the thread
    team (e.g. synchronizing a team spread across several NUMA nodes costs
    more per thread than a small, single-socket one), so a single
    process-wide calibration would be systematically wrong for some team
    sizes. Calibrating separately for every distinct ``n_threads`` value
    would fix that, but would also mean most calls pay for a fresh, uncached
    calibration. As a middle ground, ``n_threads`` is rounded down to the
    closest power of 2: this keeps the number of distinct calibrations (and
    thus their one-off cost) small.
    """
    cdef int bucket = 1 << ((max(1, n_threads)).bit_length() - 1)

    global _OVERHEAD_PER_THREAD_CACHE
    if bucket not in _OVERHEAD_PER_THREAD_CACHE:
        _OVERHEAD_PER_THREAD_CACHE[bucket] = (
            _calibrate_overhead_per_thread(bucket)
        )
    return _OVERHEAD_PER_THREAD_CACHE[bucket] * n_threads


def _use_threads_for_workload(workload_size, n_threads):
    """Decide whether a workload is worth parallelizing over ``n_threads``.

    ``workload_size`` is the estimated size of the workload, in ~simple-instructions
    (callers derive it as roughly ``n_work_items * ops_per_item``;
     we assume ~2e9 such simple instructions per second).
    This compares the estimated serial runtime against the estimated parallel runtime
    and only recommends threading if that is actually faster.
    """
    if n_threads == 1:
        return False
    t_serial = workload_size / 2e9
    t_parallelized = (
        t_serial / n_threads
        + _omp_prange_dispatch_overhead(n_threads)
    )

    return t_parallelized < t_serial


def _optimal_n_threads_for_workload(workload_size, max_n_threads):
    t_serial = workload_size / 2e9
    best_t = t_serial
    best_n_threads = 1
    for n_threads in range(2, max_n_threads + 1):
        t_parallelized = (
            t_serial / n_threads
            + _omp_prange_dispatch_overhead(n_threads)
        )
        if t_parallelized > t_serial:
            break
        elif t_parallelized < best_t:
            best_t = best_t
            best_n_threads = n_threads

    return best_n_threads


def _bench_prange_dispatch(int n_threads, Py_ssize_t n_repeats):
    """Time ``n_repeats`` dispatches of a near-empty OpenMP parallel region."""
    cdef:
        Py_ssize_t _r, _i, j
        Py_ssize_t _sink = 0
    t0 = perf_counter()
    with nogil:
        for _r in range(n_repeats):
            for _i in prange(n_threads, schedule='static', num_threads=n_threads):
                pass
            # Trivial sequential iterations giving worker threads a short idle period:
            for j in range(5):
                _sink = j
            # This is to make sure we measure the active wait behavior, see gh-34764
    return perf_counter() - t0


def _calibrate_overhead_per_thread(n_threads):
    """Measure the actual cost, in seconds, of dispatching an OpenMP
    parallel region with a team of ``n_threads`` threads on this machine.

    ``n_repeats`` and ``n_measures`` below were empirically found to give
    stable, representative measurements while keeping this cheap: this runs
    synchronously on first use (e.g. the first ``fit``/``predict`` call in a
    process), and larger ``n_threads`` teams are themselves more expensive to
    dispatch, so both ``repeats`` and the number of draws are kept modest
    rather than maximizing robustness on noisy machines. On a busy/many-core
    machine the estimate could in principle keep trending downward well past
    this many draws, but measurements across a wide range of team sizes
    showed the minimum stabilizing well before that, so this trades a bit of
    that robustness for bounding the worst-case (largest-team) calibration
    cost to roughly a hundred milliseconds rather than close to a second.

    Falls back to ``_DEFAULT_OVERHEAD_PER_THREAD`` if OpenMP is
    disabled or the measurement fails.
    """
    if not SKLEARN_OPENMP_PARALLELISM_ENABLED:
        return _DEFAULT_OVERHEAD_PER_THREAD

    try:
        n_repeats = 30
        n_measures = 15

        # The very first OpenMP parallel region dispatched in a process pays
        # for one-time thread pool creation (tens of us), on top of the
        # steady-state dispatch cost this is trying to measure. Warm it up,
        # untimed, before measuring.
        _bench_prange_dispatch(n_threads, n_repeats)

        # Take several independent measurements and keep the minimum:
        # external interference (CPU frequency scaling, other processes
        # scheduler preemption, ...) can only ever slow a measurement down,
        # never speed it up below the true steady-state dispatch cost,
        # so the minimum across samples is the one least
        # contaminated by such noise (the same reasoning `timeit` uses).
        return min(
            _bench_prange_dispatch(n_threads, n_repeats) / n_repeats
            for _ in range(n_measures)
        ) / n_threads
    except Exception:
        return _DEFAULT_OVERHEAD_PER_THREAD


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
