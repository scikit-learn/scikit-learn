import os

import pytest
from threadpoolctl import threadpool_limits

from sklearn.utils._openmp_helpers import _openmp_effective_n_threads


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
    matplotlib.use('agg', warn=False, force=True)
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
