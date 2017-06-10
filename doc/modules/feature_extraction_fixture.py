"""Fixture module to skip the feature_extraction docs when pandas is not
installed

"""
from sklearn.utils.testing import SkipTest


def setup(module):
    try:
        import pandas  # noqa
    except ImportError:
        raise SkipTest("pandas not installed")
