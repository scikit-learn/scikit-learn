import os
import re
import subprocess
import sys
import warnings


# Module level cache for cpu_count as we do not expect this to change during
# the lifecycle of a Python program. This dictionary is keyed by
# only_physical_cores.
_CPU_COUNTS = {}

# Module level cache for _openmp_uses_active_wait. Its computation is costly
# because it spawns a subprocess. None means "not computed yet".
_ACTIVE_WAIT_CACHE = None

# Lets users tell the answer to _openmp_uses_active_wait instead of detecting it.
_WAIT_POLICY_ENV_VAR = "SKLEARN_OMP_WAIT_POLICY"


# Used by _openmp_uses_active_wait to parse the settings that the OpenMP runtime
# prints on stderr when OMP_DISPLAY_ENV=VERBOSE. Each pattern matches a value
# meaning that idle threads sleep right away instead of spinning, e.g.
# "[host] KMP_BLOCKTIME='0ms'" (libomp) or "GOMP_SPINCOUNT = '1'" (libgomp).
# OMP_WAIT_POLICY is not checked: libgomp can report it as PASSIVE while
# a large GOMP_SPINCOUNT still makes threads spin.
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
    """Whether idle OpenMP threads spin (active wait) rather than sleep.

    OpenMP keeps its threads alive for the whole process. Between two parallel
    regions (the body of a ``prange`` loop or of a ``with parallel()`` block), only
    the main thread works: the other threads either spin on their CPU core for a
    while ("active wait") or go to sleep ("passive wait"). Handing work to spinning
    threads at the start of the next parallel region is fast, while waking up
    sleeping threads is slow. So with passive wait, multi-threading only pays off on
    larger workloads.

    How long idle threads spin before sleeping is set by ``GOMP_SPINCOUNT`` for GNU's
    ``libgomp`` (a number of iterations) and by ``KMP_BLOCKTIME`` for Intel's and
    LLVM's ``libomp`` (a duration). Their defaults depend on the runtime, its version,
    its build options and even the CPU, and ``OMP_WAIT_POLICY`` changes them too. So
    rather than guessing from environment variables, this starts a short-lived
    subprocess with ``OMP_DISPLAY_ENV=VERBOSE``, in which the runtime prints the
    values it uses, and looks for one meaning passive wait (see
    ``_PASSIVE_WAIT_PATTERNS``).

    The result is cached for the process, and written to ``SKLEARN_OMP_WAIT_POLICY``
    in ``os.environ`` so that processes started later (e.g. joblib workers) skip the
    detection.

    If ``SKLEARN_OMP_WAIT_POLICY`` is set to ``"active"`` or ``"passive"``
    (case-insensitive), its value is used without any detection. It does not change
    the wait policy of the OpenMP runtime itself.

    Active wait is assumed when the detection is not possible, fails or finds nothing:
    in a frozen application (no Python interpreter to run the subprocess), when the
    subprocess fails or times out (with a warning pointing to
    ``SKLEARN_OMP_WAIT_POLICY``), or when its output has no value meaning passive
    wait. The latter is not a failure: some runtimes, such as MSVC's ``vcomp`` (which
    waits actively by default), do not support ``OMP_DISPLAY_ENV``.
    """
    global _ACTIVE_WAIT_CACHE

    if not SKLEARN_OPENMP_PARALLELISM_ENABLED:
        return False

    wait_policy = os.environ.get(_WAIT_POLICY_ENV_VAR)
    if wait_policy is not None:
        if wait_policy.lower() not in ("active", "passive"):
            raise ValueError(
                f"{_WAIT_POLICY_ENV_VAR} should be 'active' or 'passive', got "
                f"{wait_policy!r}."
            )
        return wait_policy.lower() == "active"

    if _ACTIVE_WAIT_CACHE is not None:
        return _ACTIVE_WAIT_CACHE

    # Only a Python interpreter can run the script below. A frozen application
    # (PyInstaller, cx_Freeze, ...) would relaunch itself instead.
    if getattr(sys, "frozen", False):
        _ACTIVE_WAIT_CACHE = True  # conservative fallback
        return _ACTIVE_WAIT_CACHE

    env = os.environ.copy()
    env["OMP_DISPLAY_ENV"] = "VERBOSE"
    # Protect against accidental recursive inspections.
    env[_WAIT_POLICY_ENV_VAR] = "active"
    # The subprocess loads only this compiled module, directly from its file:
    # `import sklearn` would pull in the whole package (~0.3 s on its own). Its
    # parent packages are stubbed in sys.modules to skip sklearn/__init__.py,
    # with the __path__ they have in this process, so that genuine subpackage
    # imports (e.g. sklearn._cyutility) resolve to the same files even when the
    # subprocess's sys.path would not find this scikit-learn (e.g. when it was
    # made importable by editing sys.path).
    #
    # It then calls _force_openmp_runtime_init: libgomp initializes (and prints
    # its settings) as soon as it is loaded, but libomp only on the first OpenMP
    # call.
    #
    # TODO: move this detection to threadpoolctl. For libomp, it could read the
    # value directly with kmp_get_blocktime, without any subprocess:
    # https://github.com/joblib/threadpoolctl/issues/270
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

    # Assume active wait by default and use that as a conservative fallback
    # when runtime detection fails. Note that on Windows, the `vcomp`
    # implementation does not offer a verbose inspection option but its
    # default behavior is aligned with that assumption.
    active_wait = True
    failure = None
    timeout = 5
    try:
        completed = subprocess.run(
            [sys.executable, "-c", import_script],
            check=False,
            capture_output=True,
            env=env,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        failure = f"timed out after {timeout} s"
    except Exception as exc:
        failure = f"failed with {type(exc).__name__}: {exc}"
    else:
        if completed.returncode != 0:
            # The last line of a traceback names the exception. Not indexed with
            # [-1], which is invalid under the wraparound=False directive.
            last_line = completed.stderr.strip().rpartition("\n")[2]
            failure = f"exited with code {completed.returncode}: {last_line}"
        else:
            output = completed.stdout + "\n" + completed.stderr
            if any(pattern.search(output) for pattern in _PASSIVE_WAIT_PATTERNS):
                active_wait = False

    if failure is not None:
        warnings.warn(
            "Could not detect whether idle OpenMP threads spin or sleep, assuming "
            f"they spin (the detection subprocess {failure}). Set the "
            f"{_WAIT_POLICY_ENV_VAR} environment variable to 'active' or 'passive' "
            "to skip this detection.",
            UserWarning,
        )

    _ACTIVE_WAIT_CACHE = active_wait
    # Avoid redundant runtime inspection in joblib workers spawned in the
    # future, assuming this policy state will not change during the life of
    # this Python process.
    os.environ[_WAIT_POLICY_ENV_VAR] = "active" if active_wait else "passive"
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
            # lazy important so that just importing `_force_openmp_runtime_init` above
            # doesn't require to import joblib
            from joblib import cpu_count

            n_cpus = cpu_count(only_physical_cores=only_physical_cores)
            _CPU_COUNTS[only_physical_cores] = n_cpus
        max_n_threads = min(omp_get_max_threads(), n_cpus)

    if n_threads is None:
        return max_n_threads
    elif n_threads < 0:
        return max(1, max_n_threads + n_threads + 1)

    return n_threads
