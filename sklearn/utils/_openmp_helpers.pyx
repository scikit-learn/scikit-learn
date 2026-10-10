import multiprocessing.spawn
import os
import re
import subprocess
import sys


# Module level cache for cpu_count as we do not expect this to change during
# the lifecycle of a Python program. This dictionary is keyed by
# only_physical_cores.
_CPU_COUNTS = {}

# Module level cache for _openmp_uses_active_wait, since computing it spawns
# a subprocess. None means "not computed yet".
_ACTIVE_WAIT_CACHE = None

# Set in the environment of the subprocess spawned by _openmp_uses_active_wait.
_ACTIVE_WAIT_PROBE_ENV_VAR = "_SKLEARN_OPENMP_ACTIVE_WAIT_PROBE"


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
    """Whether idle OpenMP threads actually spin (busy-wait) rather than sleep while
    waiting for work.

    A thread that actively waits can rejoin a parallel region with much lower latency
    than one that has gone to sleep, which makes it worth spinning up more threads for
    smaller workloads. This is controlled by ``GOMP_SPINCOUNT`` (GNU's ``libgomp``) and
    ``KMP_BLOCKTIME`` (Intel's and LLVM's OpenMP runtimes).

    This asks the runtime for its own resolved settings: it spawns a short-lived
    subprocess with ``OMP_DISPLAY_ENV=VERBOSE``, which makes the OpenMP runtime print
    the settings scans that output for ``_PASSIVE_WAIT_PATTERNS`` above.

    If no Python interpreter is available to run the subprocess (e.g. in a frozen
    application), or if it fails, times out, or none of these patterns are found, this
    conservatively assumes active waiting is enabled.

    Spawning a subprocess is too costly to redo on every call, so the result is cached
    at the module level after the first call. The subprocess itself is kept as cheap as
    possible: rather than ``import sklearn`` (which pulls in the whole package, ~0.3s on
    its own), it loads only this compiled extension module directly from its file,
    stubbing out its parent packages in ``sys.modules`` (with the ``__path__`` they
    have in this process, so genuine subpackage imports such as ``sklearn._cyutility``
    resolve to the same files) to skip running ``sklearn/__init__.py``.
    """
    global _ACTIVE_WAIT_CACHE

    if _ACTIVE_WAIT_CACHE is not None:
        return _ACTIVE_WAIT_CACHE

    if not SKLEARN_OPENMP_PARALLELISM_ENABLED:
        _ACTIVE_WAIT_CACHE = False
        return _ACTIVE_WAIT_CACHE

    # Only a Python interpreter can run the script below. A frozen application
    # (PyInstaller, cx_Freeze, ...) would relaunch itself instead, and so would
    # an application embedding Python (uWSGI, Blender, ...) unless it was pointed
    # to an interpreter with `multiprocessing.set_executable`, which this honors
    # like multiprocessing does.
    if getattr(sys, "frozen", False):
        executable = None
    else:
        executable = multiprocessing.spawn.get_executable()
    # If the executable is not an interpreter after all, the application it
    # launches must not probe again in turn, or an application fitting a HGB
    # model on startup would keep relaunching itself.
    if os.environ.get(_ACTIVE_WAIT_PROBE_ENV_VAR):
        executable = None

    env = os.environ.copy()
    env["OMP_DISPLAY_ENV"] = "VERBOSE"
    env[_ACTIVE_WAIT_PROBE_ENV_VAR] = "1"
    # Forces the OpenMP runtime to actually initialize (and thus dump its
    # resolved environment) in this fresh subprocess: GNU libgomp
    # initializes eagerly on load, but the LLVM/Intel libomp runtime only
    # initializes lazily, on the first real OpenMP call (triggered here by
    # _force_openmp_runtime_init).
    #
    # The parent packages' __path__ is taken from this process rather than
    # looked up in the subprocess, whose sys.path may not find the same
    # scikit-learn (e.g. when it was made importable by editing sys.path).
    parts = __name__.split(".")
    parent_paths = [
        (name, list(sys.modules[name].__path__))
        for name in (".".join(parts[:i]) for i in range(1, len(parts)))
    ]
    import_script = (
        "import sys, types, importlib.util\n"
        f"for name, path in {parent_paths!r}:\n"
        "    sys.modules[name] = types.ModuleType(name)\n"
        "    sys.modules[name].__path__ = path\n"
        f"spec = importlib.util.spec_from_file_location({__name__!r}, {__file__!r})\n"
        "mod = importlib.util.module_from_spec(spec)\n"
        f"sys.modules[{__name__!r}] = mod\n"
        "spec.loader.exec_module(mod)\n"
        "mod._force_openmp_runtime_init()\n"
    )

    active_wait = True  # conservative fallback
    completed = None
    if executable:
        try:
            completed = subprocess.run(
                [os.fsdecode(executable), "-c", import_script],
                check=False,
                capture_output=True,
                env=env,
                text=True,
                timeout=5,
            )
        except Exception:
            pass

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


cpdef _force_openmp_runtime_init():
    """Make a real OpenMP call, to force the runtime to initialize.

    Used by ``_openmp_uses_active_wait``, in a fresh subprocess, to make
    lazily-initializing runtimes (LLVM/Intel libomp) resolve and display
    their settings. Kept free of any dependency beyond this compiled
    extension module, so that subprocess can be spawned without paying for
    a full ``import sklearn``.
    """
    omp_get_max_threads()


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
            from joblib import cpu_count

            n_cpus = cpu_count(only_physical_cores=only_physical_cores)
            _CPU_COUNTS[only_physical_cores] = n_cpus
        max_n_threads = min(omp_get_max_threads(), n_cpus)

    if n_threads is None:
        return max_n_threads
    elif n_threads < 0:
        return max(1, max_n_threads + n_threads + 1)

    return n_threads
