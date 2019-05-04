import pytest


@pytest.fixture(scope='function')
def pyplot():
    """Setup and teardown fixture for matplotlib.

    This fixture checks if we can import matplotlib. If not, the tests will be
    skipped. Otherwise, we setup matplotlib backend and close the figures
    after running the functions.

    Returns
    -------
    plt : module
        The ``matplotlib.pyplot`` module.
    """
    pyplot = pytest.importorskip('matplotlib.pyplot')
    import matplotlib
    matplotlib.use('agg', warn=False, force=True)
    yield pyplot
    pyplot.close('all')
