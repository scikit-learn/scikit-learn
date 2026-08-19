import os
import re
import subprocess
import sys
from joblib import cpu_count


# Module level cache for cpu_count as we do not expect this to change during
# the lifecycle of a Python program. This dictionary is keyed by
# only_physical_cores.
_CPU_COUNTS = {}

# Module level cache for _openmp_uses_active_wait, since computing it spawns
# a subprocess. None means "not computed yet".
_ACTIVE_WAIT_CACHE = None


def _min_instructions_per_thread():
    """Minimum amount of work (in ~simple-instruction units) that must be
    available per thread before bothering to parallelize a loop with
    OpenMP; see ``_use_threads_for_workload`` in ``_openmp_helpers.pxd``.

    Meant to be called once, at import time, by modules that use
    ``_use_threads_for_workload``, and cached in a module-level constant
    there. Can be overridden with the
    ``SKLEARN_MIN_INSTRUCTIONS_PER_THREAD`` environment variable, for
    experimentation.

    Defaults to 2000 if OpenMP threads actively wait (see
    ``_openmp_uses_active_wait``), since they can then rejoin a parallel
    region cheaply. Without active waiting, spinning up threads for small
    workloads is comparatively more costly, so the default is raised to
    20000.
    """
    default = 2000 if _openmp_uses_active_wait() else 20_000
    return int(os.environ.get("SKLEARN_MIN_INSTRUCTIONS_PER_THREAD", default))


# Patterns matched against the block an OpenMP runtime prints on stderr when
# OMP_DISPLAY_ENV is set: they identify a *resolved* spin-count setting (as
# opposed to what was requested through environment variables, which the
# runtime may ignore, normalize, or default differently) that means idle
# threads do not spin. Note that ``OMP_WAIT_POLICY`` is not a reliable
# signal by itself: e.g. libgomp can report it as ``PASSIVE`` while
# ``GOMP_SPINCOUNT`` was explicitly set to a large value, which still makes
# threads spin in practice. Values may be prefixed by a device scope (e.g.
# "[host] ") and wrapped in quotes, hence the loose whitespace/quoting.
_PASSIVE_WAIT_PATTERNS = tuple(
    re.compile(pattern)
    for pattern in (
        # GNU libgomp: spinning at most once is effectively passive.
        r"GOMP_SPINCOUNT\s*=\s*'?[01]'?\D",
        # Intel's and LLVM's OpenMP runtimes.
        r"KMP_BLOCKTIME\s*=\s*'?0'?\D",
    )
)


def _openmp_uses_active_wait():
    """Whether idle OpenMP threads actually spin (busy-wait) rather than
    sleep while waiting for work.

    A thread that actively waits can rejoin a parallel region with much
    lower latency than one that has gone to sleep, which makes it worth
    spinning up more threads for smaller workloads. This is controlled by
    ``GOMP_SPINCOUNT`` (GNU's ``libgomp``) and ``KMP_BLOCKTIME`` (Intel's
    and LLVM's OpenMP runtimes).

    Reading these directly from ``os.environ`` is unreliable: which
    variable applies depends on which OpenMP runtime is actually loaded,
    and each has its own default (that can change across versions/vendors)
    when unset. Instead, this asks the runtime for its own resolved
    settings: it spawns a short-lived subprocess with
    ``OMP_DISPLAY_ENV=VERBOSE``, which makes the OpenMP runtime print the
    settings it actually resolved to on stderr, and scans that output for
    ``_PASSIVE_WAIT_PATTERNS`` above.

    If the subprocess fails, times out, or none of these patterns are
    found, this conservatively assumes active waiting is enabled.

    Spawning a subprocess is too costly to redo on every call, so the
    result is cached at the module level after the first call.
    """
    global _ACTIVE_WAIT_CACHE

    if _ACTIVE_WAIT_CACHE is not None:
        return _ACTIVE_WAIT_CACHE

    if not SKLEARN_OPENMP_PARALLELISM_ENABLED:
        _ACTIVE_WAIT_CACHE = False
        return _ACTIVE_WAIT_CACHE

    env = os.environ.copy()
    env["OMP_DISPLAY_ENV"] = "VERBOSE"
    # Forces the OpenMP runtime to actually initialize (and thus dump its
    # resolved environment) in this fresh subprocess: GNU libgomp
    # initializes eagerly on load, but the LLVM/Intel libomp runtime only
    # initializes lazily, on the first real OpenMP call.
    import_script = (
        "from sklearn.utils._openmp_helpers import _openmp_effective_n_threads; "
        "_openmp_effective_n_threads()"
    )

    active_wait = True  # conservative fallback
    try:
        completed = subprocess.run(
            [sys.executable, "-c", import_script],
            check=False,
            capture_output=True,
            env=env,
            text=True,
            timeout=30,
        )
    except Exception:
        completed = None

    if completed is not None:
        output = completed.stdout + "\n" + completed.stderr
        if any(pattern.search(output) for pattern in _PASSIVE_WAIT_PATTERNS):
            active_wait = False

    _ACTIVE_WAIT_CACHE = active_wait
    return _ACTIVE_WAIT_CACHE


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
