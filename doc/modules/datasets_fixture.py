"""Fixture module to skip the datasets loading when offline

Doctests are skipped if the datasets have not already been dowloaded
and cached in the past.
"""
from os.path import exists
from os.path import join
from os import makedirs
from nose import SkipTest
from scikits.learn import datasets
from scikits.learn.datasets import get_data_home
from scikits.learn.utils.testing import mock_urllib2
import tempfile
import scipy as sp
import shutil


def globs(globs):
    # setup mock urllib2 module to avoid downloading from mldata.org
    mock_datasets = {
        'mnist-original': {
            'data': sp.empty((70000, 784)),
            'label': sp.repeat(sp.arange(10, dtype='d'), 7000),
        },
        'iris': {
            'data': sp.empty((150, 4)),
        },
        'datasets-uci-iris': {
            'double0': sp.empty((150, 4)),
            'class': sp.empty((150,)),
        },
    }

    global custom_data_home
    custom_data_home = tempfile.mkdtemp()
    makedirs(join(custom_data_home, 'mldata'))
    globs['custom_data_home'] = custom_data_home

    global _urllib2_ref
    _urllib2_ref = datasets.mldata.urllib2
    globs['_urllib2_ref'] = _urllib2_ref
    globs['mock_urllib2'] = mock_urllib2(mock_datasets)
    return globs


def setup_module(module):
    data_home = get_data_home()
    if (not exists(join(data_home, 'lfw_home'))
        or not exists(join(data_home, '20news_home'))):
        raise SkipTest("Skipping dataset loading doctests")


def teardown_module(module):
    datasets.mldata.urllib2 = _urllib2_ref
    shutil.rmtree(custom_data_home)
