import os
from os import environ
from functools import wraps

import pytest
from threadpoolctl import threadpool_limits

from sklearn.utils._openmp_helpers import _openmp_effective_n_threads
from sklearn.datasets import fetch_20newsgroups
from sklearn.datasets import fetch_20newsgroups_vectorized
from sklearn.datasets import fetch_california_housing
from sklearn.datasets import fetch_covtype
from sklearn.datasets import fetch_kddcup99
from sklearn.datasets import fetch_olivetti_faces
from sklearn.datasets import fetch_rcv1


dataset_fetchers = {
    'fetch_20newsgroups_fxt': fetch_20newsgroups,
    'fetch_20newsgroups_vectorized_fxt': fetch_20newsgroups_vectorized,
    'fetch_california_housing_fxt': fetch_california_housing,
    'fetch_covtype_fxt': fetch_covtype,
    'fetch_kddcup99_fxt': fetch_kddcup99,
    'fetch_olivetti_faces_fxt': fetch_olivetti_faces,
    'fetch_rcv1_fxt': fetch_rcv1,
}


def _fetch_fixture(f):
    """Fetch dataset (download if missing and requested by environment)."""
    download_if_missing = environ.get('SKLEARN_SKIP_NETWORK_TESTS', '1') == '0'

    @wraps(f)
    def wrapped(*args, **kwargs):
        kwargs['download_if_missing'] = download_if_missing
        try:
            return f(*args, **kwargs)
        except IOError:
            pytest.skip("test is enabled when SKLEARN_SKIP_NETWORK_TESTS=0")
    return pytest.fixture(lambda: wrapped)


# Adds fixtures for fetching data
fetch_20newsgroups_fxt = _fetch_fixture(fetch_20newsgroups)
fetch_20newsgroups_vectorized_fxt = \
    _fetch_fixture(fetch_20newsgroups_vectorized)
fetch_california_housing_fxt = _fetch_fixture(fetch_california_housing)
fetch_covtype_fxt = _fetch_fixture(fetch_covtype)
fetch_kddcup99_fxt = _fetch_fixture(fetch_kddcup99)
fetch_olivetti_faces_fxt = _fetch_fixture(fetch_olivetti_faces)
fetch_rcv1_fxt = _fetch_fixture(fetch_rcv1)


def pytest_collection_modifyitems(config, items):
    """Called after collect is completed.

    Parameters
    ----------
    config : pytest config
    items : list of collected items
    """
    run_network_tests = environ.get('SKLEARN_SKIP_NETWORK_TESTS', '1') == '0'
    skip_network = pytest.mark.skip(
        reason="test is enabled when SKLEARN_SKIP_NETWORK_TESTS=0")

    # download datasets during collection to avoid thread unsafe behavior
    # when running pytest in parallel with pytest-xdist
    dataset_features_set = set(dataset_fetchers)
    datasets_to_download = set()

    for item in items:
        if not hasattr(item, "fixturenames"):
            continue
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
            dataset_fetchers[name]()


@pytest.fixture(scope='function')
def pyplot():
    """Setup and teardown fixture for matplotlib.

    This fixture checks if we can import matplotlib. If not, the tests will be
    skipped. Otherwise, we setup matplotlib backend and close the figures
    after running the functions.

    Returns
    -------
    pyplot : module
        The ``matplotlib.pyplot`` module.
    """
    matplotlib = pytest.importorskip('matplotlib')
    matplotlib.use('agg')
    pyplot = pytest.importorskip('matplotlib.pyplot')
    yield pyplot
    pyplot.close('all')


def pytest_runtest_setup(item):
    """Set the number of openmp threads based on the number of workers
    xdist is using to prevent oversubscription.

    Parameters
    ----------
    item : pytest item
        item to be processed
    """
    try:
        xdist_worker_count = int(os.environ['PYTEST_XDIST_WORKER_COUNT'])
    except KeyError:
        # raises when pytest-xdist is not installed
        return

    openmp_threads = _openmp_effective_n_threads()
    threads_per_worker = max(openmp_threads // xdist_worker_count, 1)
    threadpool_limits(threads_per_worker, user_api='openmp')
