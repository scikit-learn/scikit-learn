# Even if empty this file is useful so that when running from the root folder
# ./sklearn is added to sys.path by pytest. See
# https://docs.pytest.org/en/latest/pythonpath.html for more details.  For
# example, this allows to build extensions in place and run pytest
# doc/modules/clustering.rst and use sklearn from the local folder rather than
# the one from site-packages.

import platform
import sys
from distutils.version import LooseVersion

import pytest
from _pytest.doctest import DoctestItem

PYTEST_MIN_VERSION = '3.3.0'

if LooseVersion(pytest.__version__) < PYTEST_MIN_VERSION:
    raise('Your version of pytest is too old, you should have at least '
          'pytest >= {} installed.'.format(PYTEST_MIN_VERSION))

def pytest_collection_modifyitems(config, items):

    # FeatureHasher is not compatible with PyPy
    if platform.python_implementation() == 'PyPy':
        skip_marker = pytest.mark.skip(
            reason='FeatureHasher is not compatible with PyPy')
        for item in items:
            if item.name == 'sklearn.feature_extraction.hashing.FeatureHasher':
                item.add_marker(skip_marker)

    # numpy changed the str/repr formatting of numpy arrays in 1.14. We want to
    # run doctests only for numpy >= 1.14. We want to skip the doctest for
    # python 2 due to unicode.
    skip_doctests = True
    if sys.version_info[:1] >= (3,):
        skip_doctests = False
    try:
        import numpy as np
        if LooseVersion(np.__version__) >= LooseVersion('1.14'):
            skip_doctests = False
    except ImportError:
        pass

    if skip_doctests:
        skip_marker = pytest.mark.skip(
            reason='doctests are only run for numpy >= 1.14 and python >= 3')

        for item in items:
            if isinstance(item, DoctestItem):
                item.add_marker(skip_marker)
