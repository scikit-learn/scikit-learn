import os
from os.path import exists
from os.path import join
from os import environ
import warnings

from sklearn.utils import IS_PYPY
from sklearn.utils._testing import SkipTest
from sklearn.utils._testing import check_skip_network
from sklearn.datasets import get_data_home
from sklearn.datasets._base import _pkl_filepath
from sklearn.datasets._twenty_newsgroups import CACHE_NAME


def setup_labeled_faces():
    data_home = get_data_home()
    if not exists(join(data_home, 'lfw_home')):
        raise SkipTest("Skipping dataset loading doctests")


def setup_rcv1():
    check_skip_network()
    # skip the test in rcv1.rst if the dataset is not already loaded
    rcv1_dir = join(get_data_home(), "RCV1")
    if not exists(rcv1_dir):
        raise SkipTest("Download RCV1 dataset to run this test.")


def setup_twenty_newsgroups():
    data_home = get_data_home()
    cache_path = _pkl_filepath(get_data_home(), CACHE_NAME)
    if not exists(cache_path):
        raise SkipTest("Skipping dataset loading doctests")


def setup_working_with_text_data():
    if IS_PYPY and os.environ.get('CI', None):
        raise SkipTest('Skipping too slow test with PyPy on CI')
    check_skip_network()
    cache_path = _pkl_filepath(get_data_home(), CACHE_NAME)
    if not exists(cache_path):
        raise SkipTest("Skipping dataset loading doctests")


def setup_loading_other_datasets():
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("Skipping loading_other_datasets.rst, "
                       "pandas not installed")

    # checks SKLEARN_SKIP_NETWORK_TESTS to see if test should run
    run_network_tests = environ.get("SKLEARN_SKIP_NETWORK_TESTS", '1') == "0"
    if not run_network_tests:
        raise SkipTest("Skipping loading_other_datasets.rst, tests can be "
                       "enabled by settting SKLEARN_SKIP_NETWORK_TESTS=0")


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
    except ImportError:
        raise SkipTest("Skipping preprocessing.rst, pandas not installed")


def setup_unsupervised_learning():
    try:
        import skimage  # noqa
    except ImportError:
        raise SkipTest("Skipping unsupervised_learning.rst, scikit-image "
                       "not installed")
    # ignore deprecation warnings from scipy.misc.face
    warnings.filterwarnings('ignore', 'The binary mode of fromstring',
                            DeprecationWarning)


def pytest_runtest_setup(item):
    fname = item.fspath.strpath
    is_index = fname.endswith('datasets/index.rst')
    if fname.endswith('datasets/labeled_faces.rst') or is_index:
        setup_labeled_faces()
    elif fname.endswith('datasets/rcv1.rst') or is_index:
        setup_rcv1()
    elif fname.endswith('datasets/twenty_newsgroups.rst') or is_index:
        setup_twenty_newsgroups()
    elif fname.endswith('tutorial/text_analytics/working_with_text_data.rst')\
            or is_index:
        setup_working_with_text_data()
    elif fname.endswith('modules/compose.rst') or is_index:
        setup_compose()
    elif IS_PYPY and fname.endswith('modules/feature_extraction.rst'):
        raise SkipTest('FeatureHasher is not compatible with PyPy')
    elif fname.endswith('datasets/loading_other_datasets.rst'):
        setup_loading_other_datasets()
    elif fname.endswith('modules/impute.rst'):
        setup_impute()
    elif fname.endswith('modules/grid_search.rst'):
        setup_grid_search()
    elif fname.endswith('modules/preprocessing.rst'):
        setup_preprocessing()
    elif fname.endswith('statistical_inference/unsupervised_learning.rst'):
        setup_unsupervised_learning()
