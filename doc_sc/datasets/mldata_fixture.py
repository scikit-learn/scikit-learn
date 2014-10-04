"""Fixture module to skip the datasets loading when offline

Mock urllib2 access to mldata.org and create a temporary data folder.
"""

from os import makedirs
from os.path import join
import numpy as np
import tempfile
import shutil

from sklearn import datasets
from sklearn.utils.testing import install_mldata_mock
from sklearn.utils.testing import uninstall_mldata_mock


def globs(globs):
    # Create a temporary folder for the data fetcher
    global custom_data_home
    custom_data_home = tempfile.mkdtemp()
    makedirs(join(custom_data_home, 'mldata'))
    globs['custom_data_home'] = custom_data_home
    return globs


def setup_module():
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


def teardown_module():
    uninstall_mldata_mock()
    shutil.rmtree(custom_data_home)
