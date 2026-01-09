# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import builtins
import platform
import sys
from contextlib import suppress
from functools import wraps
from os import environ
from unittest import SkipTest

import joblib
import numpy as np
import pytest
from _pytest.doctest import DoctestItem
from scipy.datasets import face
from threadpoolctl import threadpool_limits

from sklearn._min_dependencies import PYTEST_MIN_VERSION
from sklearn.datasets import (
    fetch_20newsgroups,
    fetch_20newsgroups_vectorized,
    fetch_california_housing,
    fetch_covtype,
    fetch_kddcup99,
    fetch_lfw_pairs,
    fetch_lfw_people,
    fetch_olivetti_faces,
    fetch_rcv1,
    fetch_species_distributions,
)
from sklearn.utils._testing import get_pytest_filterwarning_lines
from sklearn.utils.fixes import (
    _IS_32BIT,
    np_base_version,
    parse_version,
    sp_version,
)

try:
    import pytest_run_parallel  # noqa:F401

    PARALLEL_RUN_AVAILABLE = True
except ImportError:
    PARALLEL_RUN_AVAILABLE = False


try:
    from scipy_doctest.conftest import dt_config
except ModuleNotFoundError:
    dt_config = None

if parse_version(pytest.__version__) < parse_version(PYTEST_MIN_VERSION):
    raise ImportError(
        f"Your version of pytest is too old. Got version {pytest.__version__}, you"
        f" should have pytest >= {PYTEST_MIN_VERSION} installed."
    )


def raccoon_face_or_skip():
    # SciPy requires network access to get data
    run_network_tests = environ.get("SKLEARN_SKIP_NETWORK_TESTS", "1") == "0"
    if not run_network_tests:
        raise SkipTest("test is enabled when SKLEARN_SKIP_NETWORK_TESTS=0")
    try:
        import pooch  # noqa: F401
    except ImportError:
        raise SkipTest("test requires pooch to be installed")

    return face(gray=True)


dataset_fetchers = {
    "fetch_20newsgroups_fxt": fetch_20newsgroups,
    "fetch_20newsgroups_vectorized_fxt": fetch_20newsgroups_vectorized,
    "fetch_california_housing_fxt": fetch_california_housing,
    "fetch_covtype_fxt": fetch_covtype,
    "fetch_kddcup99_fxt": fetch_kddcup99,
    "fetch_lfw_pairs_fxt": fetch_lfw_pairs,
    "fetch_lfw_people_fxt": fetch_lfw_people,
    "fetch_olivetti_faces_fxt": fetch_olivetti_faces,
    "fetch_rcv1_fxt": fetch_rcv1,
    "fetch_species_distributions_fxt": fetch_species_distributions,
}

dataset_fetchers["raccoon_face_fxt"] = raccoon_face_or_skip

_SKIP32_MARK = pytest.mark.skipif(
    environ.get("SKLEARN_RUN_FLOAT32_TESTS", "0") != "1",
    reason="Set SKLEARN_RUN_FLOAT32_TESTS=1 to run float32 dtype tests",
)


# Global fixtures
@pytest.fixture(params=[pytest.param(np.float32, marks=_SKIP32_MARK), np.float64])
def global_dtype(request):
    yield request.param


def _fetch_fixture(f):
    """Fetch dataset (download if missing and requested by environment)."""
    download_if_missing = environ.get("SKLEARN_SKIP_NETWORK_TESTS", "1") == "0"

    @wraps(f)
    def wrapped(*args, **kwargs):
        kwargs["download_if_missing"] = download_if_missing
        try:
            return f(*args, **kwargs)
        except OSError as e:
            if str(e) != "Data not found and `download_if_missing` is False":
                raise
            pytest.skip("test is enabled when SKLEARN_SKIP_NETWORK_TESTS=0")

    return pytest.fixture(lambda: wrapped)


# Adds fixtures for fetching data
fetch_20newsgroups_fxt = _fetch_fixture(fetch_20newsgroups)
fetch_20newsgroups_vectorized_fxt = _fetch_fixture(fetch_20newsgroups_vectorized)
fetch_california_housing_fxt = _fetch_fixture(fetch_california_housing)
fetch_covtype_fxt = _fetch_fixture(fetch_covtype)
fetch_kddcup99_fxt = _fetch_fixture(fetch_kddcup99)
fetch_lfw_pairs_fxt = _fetch_fixture(fetch_lfw_pairs)
fetch_lfw_people_fxt = _fetch_fixture(fetch_lfw_people)
fetch_olivetti_faces_fxt = _fetch_fixture(fetch_olivetti_faces)
fetch_rcv1_fxt = _fetch_fixture(fetch_rcv1)
fetch_species_distributions_fxt = _fetch_fixture(fetch_species_distributions)
raccoon_face_fxt = pytest.fixture(raccoon_face_or_skip)


def pytest_collection_modifyitems(config, items):
    """Called after collect is completed.

    Parameters
    ----------
    config : pytest config
    items : list of collected items
    """
    run_network_tests = environ.get("SKLEARN_SKIP_NETWORK_TESTS", "1") == "0"
    skip_network = pytest.mark.skip(
        reason="test is enabled when SKLEARN_SKIP_NETWORK_TESTS=0"
    )

    # download datasets during collection to avoid thread unsafe behavior
    # when running pytest in parallel with pytest-xdist
    dataset_features_set = set(dataset_fetchers)
    datasets_to_download = set()

    for item in items:
        if isinstance(item, DoctestItem) and "fetch_" in item.name:
            fetcher_function_name = item.name.split(".")[-1]
            dataset_fetchers_key = f"{fetcher_function_name}_fxt"
            dataset_to_fetch = set([dataset_fetchers_key]) & dataset_features_set
        elif not hasattr(item, "fixturenames"):
            continue
        else:
            item_fixtures = set(item.fixturenames)
            dataset_to_fetch = item_fixtures & dataset_features_set

        if not dataset_to_fetch:
            continue

        if run_network_tests:
            datasets_to_download |= dataset_to_fetch
        else:
            # network tests are skipped
            item.add_marker(skip_network)

    # Only download datasets on the first worker spawned by pytest-xdist
    # to avoid thread unsafe behavior. If pytest-xdist is not used, we still
    # download before tests run.
    worker_id = environ.get("PYTEST_XDIST_WORKER", "gw0")
    if worker_id == "gw0" and run_network_tests:
        for name in datasets_to_download:
            with suppress(SkipTest):
                dataset_fetchers[name]()

    for item in items:
        # Known failure on with GradientBoostingClassifier on ARM64
        if (
            item.name.endswith("GradientBoostingClassifier")
            and platform.machine() == "aarch64"
        ):
            marker = pytest.mark.xfail(
                reason=(
                    "know failure. See "
                    "https://github.com/scikit-learn/scikit-learn/issues/17797"
                )
            )
            item.add_marker(marker)

    skip_doctests = False
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        skip_doctests = True
        reason = "matplotlib is required to run the doctests"

    if _IS_32BIT:
        reason = "doctest are only run when the default numpy int is 64 bits."
        skip_doctests = True
    elif sys.platform.startswith("win32"):
        reason = (
            "doctests are not run for Windows because numpy arrays "
            "repr is inconsistent across platforms."
        )
        skip_doctests = True

    if np_base_version < parse_version("2"):
        # TODO: configure numpy to output scalar arrays as regular Python scalars
        # once possible to improve readability of the tests docstrings.
        # https://numpy.org/neps/nep-0051-scalar-representation.html#implementation
        reason = "Due to NEP 51 numpy scalar repr has changed in numpy 2"
        skip_doctests = True

    if sp_version < parse_version("1.14"):
        reason = "Scipy sparse matrix repr has changed in scipy 1.14"
        skip_doctests = True

    # Normally doctest has the entire module's scope. Here we set globs to an empty dict
    # to remove the module's scope:
    # https://docs.python.org/3/library/doctest.html#what-s-the-execution-context
    for item in items:
        if isinstance(item, DoctestItem):
            item.dtest.globs = {}

    if skip_doctests:
        skip_marker = pytest.mark.skip(reason=reason)

        for item in items:
            if isinstance(item, DoctestItem):
                # work-around an internal error with pytest if adding a skip
                # mark to a doctest in a contextmanager, see
                # https://github.com/pytest-dev/pytest/issues/8796 for more
                # details.
                if item.name != "sklearn._config.config_context":
                    item.add_marker(skip_marker)
    try:
        import PIL  # noqa: F401

        pillow_installed = True
    except ImportError:
        pillow_installed = False

    if not pillow_installed:
        skip_marker = pytest.mark.skip(reason="pillow (or PIL) not installed!")
        for item in items:
            if item.name in [
                "sklearn.feature_extraction.image.PatchExtractor",
                "sklearn.feature_extraction.image.extract_patches_2d",
            ]:
                item.add_marker(skip_marker)


@pytest.fixture(scope="function")
def pyplot():
    """Setup and teardown fixture for matplotlib.

    This fixture checks if we can import matplotlib. If not, the tests will be
    skipped. Otherwise, we close the figures before and after running the
    functions.

    Returns
    -------
    pyplot : module
        The ``matplotlib.pyplot`` module.
    """
    pyplot = pytest.importorskip("matplotlib.pyplot")
    pyplot.close("all")
    yield pyplot
    pyplot.close("all")


def pytest_generate_tests(metafunc):
    """Parametrization of global_random_seed fixture

    based on the SKLEARN_TESTS_GLOBAL_RANDOM_SEED environment variable.

    The goal of this fixture is to prevent tests that use it to be sensitive
    to a specific seed value while still being deterministic by default.

    See the documentation for the SKLEARN_TESTS_GLOBAL_RANDOM_SEED
    variable for instructions on how to use this fixture.

    https://scikit-learn.org/dev/computing/parallelism.html#sklearn-tests-global-random-seed

    """
    # When using pytest-xdist this function is called in the xdist workers.
    # We rely on SKLEARN_TESTS_GLOBAL_RANDOM_SEED environment variable which is
    # set in before running pytest and is available in xdist workers since they
    # are subprocesses.
    RANDOM_SEED_RANGE = list(range(100))  # All seeds in [0, 99] should be valid.
    random_seed_var = environ.get("SKLEARN_TESTS_GLOBAL_RANDOM_SEED")

    default_random_seeds = [42]

    if random_seed_var is None:
        random_seeds = default_random_seeds
    elif random_seed_var == "all":
        random_seeds = RANDOM_SEED_RANGE
    else:
        if "-" in random_seed_var:
            start, stop = random_seed_var.split("-")
            random_seeds = list(range(int(start), int(stop) + 1))
        else:
            random_seeds = [int(random_seed_var)]

        if min(random_seeds) < 0 or max(random_seeds) > 99:
            raise ValueError(
                "The value(s) of the environment variable "
                "SKLEARN_TESTS_GLOBAL_RANDOM_SEED must be in the range [0, 99] "
                f"(or 'all'), got: {random_seed_var}"
            )

    if "global_random_seed" in metafunc.fixturenames:
        metafunc.parametrize("global_random_seed", random_seeds)


def pytest_addoption(parser, pluginmanager):
    if not PARALLEL_RUN_AVAILABLE:
        parser.addini("thread_unsafe_fixtures", "list of stuff")


def pytest_configure(config):
    # Use matplotlib agg backend during the tests including doctests
    try:
        import matplotlib

        matplotlib.use("agg")
    except ImportError:
        pass

    allowed_parallelism = joblib.cpu_count(only_physical_cores=True)
    xdist_worker_count = environ.get("PYTEST_XDIST_WORKER_COUNT")
    if xdist_worker_count is not None:
        # Set the number of OpenMP and BLAS threads based on the number of workers
        # xdist is using to prevent oversubscription.
        allowed_parallelism = max(allowed_parallelism // int(xdist_worker_count), 1)
    threadpool_limits(allowed_parallelism)

    if environ.get("SKLEARN_WARNINGS_AS_ERRORS", "0") != "0":
        # This seems like the only way to programmatically change the config
        # filterwarnings. This was suggested in
        # https://github.com/pytest-dev/pytest/issues/3311#issuecomment-373177592
        for line in get_pytest_filterwarning_lines():
            config.addinivalue_line("filterwarnings", line)

    if not PARALLEL_RUN_AVAILABLE:
        config.addinivalue_line(
            "markers",
            "parallel_threads(n): run the given test function in parallel "
            "using `n` threads.",
        )
        config.addinivalue_line(
            "markers",
            "thread_unsafe: mark the test function as single-threaded",
        )
        config.addinivalue_line(
            "markers",
            "iterations(n): run the given test function `n` times in each thread",
        )
        config.addinivalue_line(
            "markers",
            "iterations(n): run the given test function `n` times in each thread",
        )


@pytest.fixture
def hide_available_pandas(monkeypatch):
    """Pretend pandas was not installed."""
    import_orig = builtins.__import__

    def mocked_import(name, *args, **kwargs):
        if name == "pandas":
            raise ImportError()
        return import_orig(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", mocked_import)


if dt_config is not None:
    # Strict mode to differentiate between 3.14 and np.float64(3.14)
    dt_config.strict_check = True
    # dt_config.rtol = 0.01
