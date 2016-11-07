"""Fixture module to skip the datasets loading when offline

Doctests are skipped if the datasets have not already been dowloaded
and cached in the past.
"""
from os.path import exists
from os.path import join
from sklearn.datasets import get_data_home
from sklearn.utils.testing import SkipTest


def setup_module(module):
    data_home = get_data_home()
    if not exists(join(data_home, 'lfw_home')):
        raise SkipTest("Skipping dataset loading doctests")
