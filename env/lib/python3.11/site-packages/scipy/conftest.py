# Pytest customization
import json
import os
import warnings
import tempfile
from contextlib import contextmanager

import numpy as np
import numpy.testing as npt
import pytest
import hypothesis

from scipy._lib._fpumode import get_fpu_mode
from scipy._lib._testutils import FPUModeChangeWarning
from scipy._lib._array_api import SCIPY_ARRAY_API, SCIPY_DEVICE
from scipy._lib import _pep440

try:
    from scipy_doctest.conftest import dt_config
    HAVE_SCPDT = True
except ModuleNotFoundError:
    HAVE_SCPDT = False


def pytest_configure(config):
    config.addinivalue_line("markers",
        "slow: Tests that are very slow.")
    config.addinivalue_line("markers",
        "xslow: mark test as extremely slow (not run unless explicitly requested)")
    config.addinivalue_line("markers",
        "xfail_on_32bit: mark test as failing on 32-bit platforms")
    try:
        import pytest_timeout  # noqa:F401
    except Exception:
        config.addinivalue_line(
            "markers", 'timeout: mark a test for a non-default timeout')
    try:
        # This is a more reliable test of whether pytest_fail_slow is installed
        # When I uninstalled it, `import pytest_fail_slow` didn't fail!
        from pytest_fail_slow import parse_duration  # type: ignore[import-not-found] # noqa:F401,E501
    except Exception:
        config.addinivalue_line(
            "markers", 'fail_slow: mark a test for a non-default timeout failure')
    config.addinivalue_line("markers",
        "skip_xp_backends(*backends, reasons=None, np_only=False, cpu_only=False): "
        "mark the desired skip configuration for the `skip_xp_backends` fixture.")


def pytest_runtest_setup(item):
    mark = item.get_closest_marker("xslow")
    if mark is not None:
        try:
            v = int(os.environ.get('SCIPY_XSLOW', '0'))
        except ValueError:
            v = False
        if not v:
            pytest.skip("very slow test; "
                        "set environment variable SCIPY_XSLOW=1 to run it")
    mark = item.get_closest_marker("xfail_on_32bit")
    if mark is not None and np.intp(0).itemsize < 8:
        pytest.xfail(f'Fails on our 32-bit test platform(s): {mark.args[0]}')

    # Older versions of threadpoolctl have an issue that may lead to this
    # warning being emitted, see gh-14441
    with npt.suppress_warnings() as sup:
        sup.filter(pytest.PytestUnraisableExceptionWarning)

        try:
            from threadpoolctl import threadpool_limits

            HAS_THREADPOOLCTL = True
        except Exception:  # observed in gh-14441: (ImportError, AttributeError)
            # Optional dependency only. All exceptions are caught, for robustness
            HAS_THREADPOOLCTL = False

        if HAS_THREADPOOLCTL:
            # Set the number of openmp threads based on the number of workers
            # xdist is using to prevent oversubscription. Simplified version of what
            # sklearn does (it can rely on threadpoolctl and its builtin OpenMP helper
            # functions)
            try:
                xdist_worker_count = int(os.environ['PYTEST_XDIST_WORKER_COUNT'])
            except KeyError:
                # raises when pytest-xdist is not installed
                return

            if not os.getenv('OMP_NUM_THREADS'):
                max_openmp_threads = os.cpu_count() // 2  # use nr of physical cores
                threads_per_worker = max(max_openmp_threads // xdist_worker_count, 1)
                try:
                    threadpool_limits(threads_per_worker, user_api='blas')
                except Exception:
                    # May raise AttributeError for older versions of OpenBLAS.
                    # Catch any error for robustness.
                    return


@pytest.fixture(scope="function", autouse=True)
def check_fpu_mode(request):
    """
    Check FPU mode was not changed during the test.
    """
    old_mode = get_fpu_mode()
    yield
    new_mode = get_fpu_mode()

    if old_mode != new_mode:
        warnings.warn(f"FPU mode changed from {old_mode:#x} to {new_mode:#x} during "
                      "the test",
                      category=FPUModeChangeWarning, stacklevel=0)


# Array API backend handling
xp_available_backends = {'numpy': np}

if SCIPY_ARRAY_API and isinstance(SCIPY_ARRAY_API, str):
    # fill the dict of backends with available libraries
    try:
        import array_api_strict
        xp_available_backends.update({'array_api_strict': array_api_strict})
        if _pep440.parse(array_api_strict.__version__) < _pep440.Version('2.0'):
            raise ImportError("array-api-strict must be >= version 2.0")
        array_api_strict.set_array_api_strict_flags(
            api_version='2023.12'
        )
    except ImportError:
        pass

    try:
        import torch  # type: ignore[import-not-found]
        xp_available_backends.update({'pytorch': torch})
        # can use `mps` or `cpu`
        torch.set_default_device(SCIPY_DEVICE)
    except ImportError:
        pass

    try:
        import cupy  # type: ignore[import-not-found]
        xp_available_backends.update({'cupy': cupy})
    except ImportError:
        pass

    try:
        import jax.numpy  # type: ignore[import-not-found]
        xp_available_backends.update({'jax.numpy': jax.numpy})
        jax.config.update("jax_enable_x64", True)
        jax.config.update("jax_default_device", jax.devices(SCIPY_DEVICE)[0])
    except ImportError:
        pass

    # by default, use all available backends
    if SCIPY_ARRAY_API.lower() not in ("1", "true"):
        SCIPY_ARRAY_API_ = json.loads(SCIPY_ARRAY_API)

        if 'all' in SCIPY_ARRAY_API_:
            pass  # same as True
        else:
            # only select a subset of backend by filtering out the dict
            try:
                xp_available_backends = {
                    backend: xp_available_backends[backend]
                    for backend in SCIPY_ARRAY_API_
                }
            except KeyError:
                msg = f"'--array-api-backend' must be in {xp_available_backends.keys()}"
                raise ValueError(msg)

if 'cupy' in xp_available_backends:
    SCIPY_DEVICE = 'cuda'

array_api_compatible = pytest.mark.parametrize("xp", xp_available_backends.values())

skip_xp_invalid_arg = pytest.mark.skipif(SCIPY_ARRAY_API,
    reason = ('Test involves masked arrays, object arrays, or other types '
              'that are not valid input when `SCIPY_ARRAY_API` is used.'))


@pytest.fixture
def skip_xp_backends(xp, request):
    """
    Skip based on the ``skip_xp_backends`` marker.

    Parameters
    ----------
    *backends : tuple
        Backends to skip, e.g. ``("array_api_strict", "torch")``.
        These are overriden when ``np_only`` is ``True``, and are not
        necessary to provide for non-CPU backends when ``cpu_only`` is ``True``.
    reasons : list, optional
        A list of reasons for each skip. When ``np_only`` is ``True``,
        this should be a singleton list. Otherwise, this should be a list
        of reasons, one for each corresponding backend in ``backends``.
        If unprovided, default reasons are used. Note that it is not possible
        to specify a custom reason with ``cpu_only``. Default: ``None``.
    np_only : bool, optional
        When ``True``, the test is skipped for all backends other
        than the default NumPy backend. There is no need to provide
        any ``backends`` in this case. To specify a reason, pass a
        singleton list to ``reasons``. Default: ``False``.
    cpu_only : bool, optional
        When ``True``, the test is skipped on non-CPU devices.
        There is no need to provide any ``backends`` in this case,
        but any ``backends`` will also be skipped on the CPU.
        Default: ``False``.
    """
    if "skip_xp_backends" not in request.keywords:
        return
    backends = request.keywords["skip_xp_backends"].args
    kwargs = request.keywords["skip_xp_backends"].kwargs
    np_only = kwargs.get("np_only", False)
    cpu_only = kwargs.get("cpu_only", False)
    if np_only:
        reasons = kwargs.get("reasons", ["do not run with non-NumPy backends."])
        reason = reasons[0]
        if xp.__name__ != 'numpy':
            pytest.skip(reason=reason)
        return
    if cpu_only:
        reason = "do not run with `SCIPY_ARRAY_API` set and not on CPU"
        if SCIPY_ARRAY_API and SCIPY_DEVICE != 'cpu':
            if xp.__name__ == 'cupy':
                pytest.skip(reason=reason)
            elif xp.__name__ == 'torch':
                if 'cpu' not in xp.empty(0).device.type:
                    pytest.skip(reason=reason)
            elif xp.__name__ == 'jax.numpy':
                for d in xp.empty(0).devices():
                    if 'cpu' not in d.device_kind:
                        pytest.skip(reason=reason)

    if backends is not None:
        reasons = kwargs.get("reasons", False)
        for i, backend in enumerate(backends):
            if xp.__name__ == backend:
                if not reasons:
                    reason = f"do not run with array API backend: {backend}"
                else:
                    reason = reasons[i]
                pytest.skip(reason=reason)


# Following the approach of NumPy's conftest.py...
# Use a known and persistent tmpdir for hypothesis' caches, which
# can be automatically cleared by the OS or user.
hypothesis.configuration.set_hypothesis_home_dir(
    os.path.join(tempfile.gettempdir(), ".hypothesis")
)

# We register two custom profiles for SciPy - for details see
# https://hypothesis.readthedocs.io/en/latest/settings.html
# The first is designed for our own CI runs; the latter also
# forces determinism and is designed for use via scipy.test()
hypothesis.settings.register_profile(
    name="nondeterministic", deadline=None, print_blob=True,
)
hypothesis.settings.register_profile(
    name="deterministic",
    deadline=None, print_blob=True, database=None, derandomize=True,
    suppress_health_check=list(hypothesis.HealthCheck),
)

# Profile is currently set by environment variable `SCIPY_HYPOTHESIS_PROFILE`
# In the future, it would be good to work the choice into dev.py.
SCIPY_HYPOTHESIS_PROFILE = os.environ.get("SCIPY_HYPOTHESIS_PROFILE",
                                          "deterministic")
hypothesis.settings.load_profile(SCIPY_HYPOTHESIS_PROFILE)


############################################################################
# doctesting stuff

if HAVE_SCPDT:

    # FIXME: populate the dict once
    @contextmanager
    def warnings_errors_and_rng(test=None):
        """Temporarily turn (almost) all warnings to errors.

        Filter out known warnings which we allow.
        """
        known_warnings = dict()

        # these functions are known to emit "divide by zero" RuntimeWarnings
        divide_by_zero = [
            'scipy.linalg.norm', 'scipy.ndimage.center_of_mass',
        ]
        for name in divide_by_zero:
            known_warnings[name] = dict(category=RuntimeWarning,
                                        message='divide by zero')

        # Deprecated stuff in scipy.signal and elsewhere
        deprecated = [
            'scipy.signal.cwt', 'scipy.signal.morlet', 'scipy.signal.morlet2',
            'scipy.signal.ricker',
            'scipy.integrate.simpson',
            'scipy.interpolate.interp2d',
        ]
        for name in deprecated:
            known_warnings[name] = dict(category=DeprecationWarning)

        from scipy import integrate
        # the funcions are known to emit IntergrationWarnings
        integration_w = ['scipy.special.ellip_normal',
                         'scipy.special.ellip_harm_2',
        ]
        for name in integration_w:
            known_warnings[name] = dict(category=integrate.IntegrationWarning,
                                        message='The occurrence of roundoff')

        # scipy.stats deliberately emits UserWarnings sometimes
        user_w = ['scipy.stats.anderson_ksamp', 'scipy.stats.kurtosistest',
                  'scipy.stats.normaltest', 'scipy.sparse.linalg.norm']
        for name in user_w:
            known_warnings[name] = dict(category=UserWarning)

        # additional one-off warnings to filter
        dct = {
            'scipy.sparse.linalg.norm':
                dict(category=UserWarning, message="Exited at iteration"),
            # tutorials
            'linalg.rst':
                dict(message='the matrix subclass is not',
                     category=PendingDeprecationWarning),
            'stats.rst':
                dict(message='The maximum number of subdivisions',
                     category=integrate.IntegrationWarning),
        }
        known_warnings.update(dct)

        # these legitimately emit warnings in examples
        legit = set('scipy.signal.normalize')

        # Now, the meat of the matter: filter warnings,
        # also control the random seed for each doctest.

        # XXX: this matches the refguide-check behavior, but is a tad strange:
        # makes sure that the seed the old-fashioned np.random* methods is
        # *NOT* reproducible but the new-style `default_rng()` *IS* repoducible.
        # Should these two be either both repro or both not repro?

        from scipy._lib._util import _fixed_default_rng
        import numpy as np
        with _fixed_default_rng():
            np.random.seed(None)
            with warnings.catch_warnings():
                if test and test.name in known_warnings:
                    warnings.filterwarnings('ignore',
                                            **known_warnings[test.name])
                    yield
                elif test and test.name in legit:
                    yield
                else:
                    warnings.simplefilter('error', Warning)
                    yield


    dt_config.user_context_mgr = warnings_errors_and_rng
    dt_config.skiplist = set([
        'scipy.linalg.LinAlgError',     # comes from numpy
        'scipy.fftpack.fftshift',       # fftpack stuff is also from numpy
        'scipy.fftpack.ifftshift',
        'scipy.fftpack.fftfreq',
        'scipy.special.sinc',           # sinc is from numpy
        'scipy.optimize.show_options',  # does not have much to doctest
        'scipy.signal.normalize',       # manipulates warnings (XXX temp skip)
        'scipy.sparse.linalg.norm',     # XXX temp skip
    ])

    # these are affected by NumPy 2.0 scalar repr: rely on string comparison
    if np.__version__ < "2":
        dt_config.skiplist.update(set([
            'scipy.io.hb_read',
            'scipy.io.hb_write',
            'scipy.sparse.csgraph.connected_components',
            'scipy.sparse.csgraph.depth_first_order',
            'scipy.sparse.csgraph.shortest_path',
            'scipy.sparse.csgraph.floyd_warshall',
            'scipy.sparse.csgraph.dijkstra',
            'scipy.sparse.csgraph.bellman_ford',
            'scipy.sparse.csgraph.johnson',
            'scipy.sparse.csgraph.yen',
            'scipy.sparse.csgraph.breadth_first_order',
            'scipy.sparse.csgraph.reverse_cuthill_mckee',
            'scipy.sparse.csgraph.structural_rank',
            'scipy.sparse.csgraph.construct_dist_matrix',
            'scipy.sparse.csgraph.reconstruct_path',
            'scipy.ndimage.value_indices',
            'scipy.stats.mstats.describe',
    ]))

    # help pytest collection a bit: these names are either private
    # (distributions), or just do not need doctesting.
    dt_config.pytest_extra_ignore = [
        "scipy.stats.distributions",
        "scipy.optimize.cython_optimize",
        "scipy.test",
        "scipy.show_config",
    ]

    dt_config.pytest_extra_xfail = {
        # name: reason
        "io.rst": "",
        "ND_regular_grid.rst": "ReST parser limitation",
        "extrapolation_examples.rst": "ReST parser limitation",
        "sampling_pinv.rst": "__cinit__ unexpected argument",
        "sampling_srou.rst": "nan in scalar_power",
        "probability_distributions.rst": "integration warning",
    }

    # tutorials
    dt_config.pseudocode = set(['integrate.nquad(func,'])
    dt_config.local_resources = {'io.rst': ["octave_a.mat"]}
############################################################################
