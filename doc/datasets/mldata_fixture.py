"""Fixture module to skip the datasets loading when offline

Mock urllib2 access to mldata.org and create a temporary data folder.
"""

import numpy as np

from sklearn.utils.testing import install_mldata_mock
from sklearn.utils.testing import uninstall_mldata_mock


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
