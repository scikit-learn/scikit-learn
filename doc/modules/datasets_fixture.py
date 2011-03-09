"""Fixture module to skip the datasets loading when offline

Doctests are skipped if the datasets have not already been dowloaded
and cached in the past.
"""
from os.path import exists
from nose import SkipTest
from scikits.learn.datasets import get_data_home

def setup_module(module):
    if not exists(get_data_home()):
        raise SkipTest("Skipping dataset loading doctests")
