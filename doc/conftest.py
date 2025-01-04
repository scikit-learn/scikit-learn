import os
import warnings
from os import environ
from os.path import exists, join

import pytest
from _pytest.doctest import DoctestItem

from sklearn.datasets import get_data_home
from sklearn.datasets._base import _pkl_filepath
from sklearn.datasets._twenty_newsgroups import CACHE_NAME
from sklearn.utils._testing import SkipTest, check_skip_network
from sklearn.utils.fixes import np_base_version, parse_version, sp_version


def setup_labeled_faces():
    data_home = get_data_home()
    if not exists(join(data_home, "lfw_home")):
        raise SkipTest("Skipping dataset loading doctests")


def setup_rcv1():
    check_skip_network()
    # skip the test in rcv1.rst if the dataset is not already loaded
    rcv1_dir = join(get_data_home(), "RCV1")
    if not exists(rcv1_dir):
        raise SkipTest("Download RCV1 dataset to run this test.")


def setup_twenty_newsgroups():
    cache_path = _pkl_filepath(get_data_home(), CACHE_NAME)
    if not exists(cache_path):
        raise SkipTest("Skipping dataset loading doctests")


def setup_working_with_text_data():
    check_skip_network()
    cache_path = _pkl_filepath(get_data_home(), CACHE_NAME)
    if not exists(cache_path):
        raise SkipTest("Skipping dataset loading doctests")


def setup_loading_other_datasets():
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("Skipping loading_other_datasets.rst, pandas not installed")

    # checks SKLEARN_SKIP_NETWORK_TESTS to see if test should run
    run_network_tests = environ.get("SKLEARN_SKIP_NETWORK_TESTS", "1") == "0"
    if not run_network_tests:
        raise SkipTest(
            "Skipping loading_other_datasets.rst, tests can be "
            "enabled by setting SKLEARN_SKIP_NETWORK_TESTS=0"
        )


def setup_compose():
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("Skipping compose.rst, pandas not installed")


def setup_impute():
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("Skipping impute.rst, pandas not installed")


def setup_grid_search():
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("Skipping grid_search.rst, pandas not installed")


def setup_preprocessing():
    try:
        import pandas  # noqa

        if parse_version(pandas.__version__) < parse_version("1.1.0"):
            raise SkipTest("Skipping preprocessing.rst, pandas version < 1.1.0")
    except ImportError:
        raise SkipTest("Skipping preprocessing.rst, pandas not installed")


def setup_unsupervised_learning():
    try:
        import skimage  # noqa
    except ImportError:
        raise SkipTest("Skipping unsupervised_learning.rst, scikit-image not installed")
    # ignore deprecation warnings from scipy.misc.face
    warnings.filterwarnings(
        "ignore", "The binary mode of fromstring", DeprecationWarning
    )


def skip_if_matplotlib_not_installed(fname):
    try:
        import matplotlib  # noqa
    except ImportError:
        basename = os.path.basename(fname)
        raise SkipTest(f"Skipping doctests for {basename}, matplotlib not installed")


def skip_if_cupy_not_installed(fname):
    try:
        import cupy  # noqa
    except ImportError:
        basename = os.path.basename(fname)
        raise SkipTest(f"Skipping doctests for {basename}, cupy not installed")


def pytest_runtest_setup(item):
    fname = item.fspath.strpath
    # normalize filename to use forward slashes on Windows for easier handling
    # later
    fname = fname.replace(os.sep, "/")

    is_index = fname.endswith("datasets/index.rst")
    if fname.endswith("datasets/labeled_faces.rst") or is_index:
        setup_labeled_faces()
    elif fname.endswith("datasets/rcv1.rst") or is_index:
        setup_rcv1()
    elif fname.endswith("datasets/twenty_newsgroups.rst") or is_index:
        setup_twenty_newsgroups()
    elif fname.endswith("modules/compose.rst") or is_index:
        setup_compose()
    elif fname.endswith("datasets/loading_other_datasets.rst"):
        setup_loading_other_datasets()
    elif fname.endswith("modules/impute.rst"):
        setup_impute()
    elif fname.endswith("modules/grid_search.rst"):
        setup_grid_search()
    elif fname.endswith("modules/preprocessing.rst"):
        setup_preprocessing()
    elif fname.endswith("statistical_inference/unsupervised_learning.rst"):
        setup_unsupervised_learning()

    rst_files_requiring_matplotlib = [
        "modules/partial_dependence.rst",
        "modules/tree.rst",
    ]
    for each in rst_files_requiring_matplotlib:
        if fname.endswith(each):
            skip_if_matplotlib_not_installed(fname)

    if fname.endswith("array_api.rst"):
        skip_if_cupy_not_installed(fname)


def pytest_configure(config):
    # Use matplotlib agg backend during the tests including doctests
    try:
        import matplotlib

        matplotlib.use("agg")
    except ImportError:
        pass


def pytest_collection_modifyitems(config, items):
    """Called after collect is completed.

    Parameters
    ----------
    config : pytest config
    items : list of collected items
    """
    skip_doctests = False
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
                item.add_marker(skip_marker)
