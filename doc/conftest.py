import os
from os.path import exists
from os.path import join
import warnings

import numpy as np

from sklearn.utils import IS_PYPY
from sklearn.utils.testing import SkipTest
from sklearn.utils.testing import check_skip_network
from sklearn.datasets import get_data_home
from sklearn.datasets.base import _pkl_filepath
from sklearn.datasets.twenty_newsgroups import CACHE_NAME
from sklearn.utils.testing import install_mldata_mock
from sklearn.utils.testing import uninstall_mldata_mock


def setup_labeled_faces():
    data_home = get_data_home()
    if not exists(join(data_home, 'lfw_home')):
        raise SkipTest("Skipping dataset loading doctests")


def setup_mldata():
    # setup mock urllib2 module to avoid downloading from mldata.org
    install_mldata_mock({
        'mnist-original': {
            'data': np.empty((70000, 784)),
            'label': np.repeat(np.arange(10, dtype='d'), 7000),
        },
        'iris': {
            'data': np.empty((150, 4)),
        },
        'datasets-uci-iris': {
            'double0': np.empty((150, 4)),
            'class': np.empty((150,)),
        },
    })


def teardown_mldata():
    uninstall_mldata_mock()


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


def setup_unsupervised_learning():
    # ignore deprecation warnings from scipy.misc.face
    warnings.filterwarnings('ignore', 'The binary mode of fromstring',
                            DeprecationWarning)


def pytest_runtest_setup(item):
    fname = item.fspath.strpath
    is_index = fname.endswith('datasets/index.rst')
    if fname.endswith('datasets/labeled_faces.rst') or is_index:
        setup_labeled_faces()
    elif fname.endswith('datasets/mldata.rst') or is_index:
        setup_mldata()
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
    elif fname.endswith('modules/impute.rst'):
        setup_impute()
    elif fname.endswith('statistical_inference/unsupervised_learning.rst'):
        setup_unsupervised_learning()


def pytest_runtest_teardown(item):
    fname = item.fspath.strpath
    if fname.endswith('datasets/mldata.rst'):
        teardown_mldata()
