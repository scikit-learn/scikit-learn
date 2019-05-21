import gc
import sys

import pytest


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


is_pypy = '__pypy__' in sys.builtin_module_names


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_call(item):
    yield
    if is_pypy:
        gc.collect()
