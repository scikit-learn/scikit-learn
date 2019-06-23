import pytest
import numpy as np
import warnings


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
if np.iinfo(np.int).dtype != np.int64:
    warnings.warn(UserWarning("Your default numpy int is a 32 bit one, it may be because you are "
                      "running a 32 bit Python (or more complicated if you're on Windows "
                      "or Mac. This causes some docstring tests to fail on your machine...."))


