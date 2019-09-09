# Even if empty this file is useful so that when running from the root folder
# ./sklearn is added to sys.path by pytest. See
# https://docs.pytest.org/en/latest/pythonpath.html for more details.  For
# example, this allows to build extensions in place and run pytest
# doc/modules/clustering.rst and use sklearn from the local folder rather than
# the one from site-packages.

import platform
from distutils.version import LooseVersion

import pytest
from _pytest.doctest import DoctestItem

from sklearn import set_config
from sklearn.utils import _IS_32BIT

PYTEST_MIN_VERSION = '3.3.0'

if LooseVersion(pytest.__version__) < PYTEST_MIN_VERSION:
    raise ImportError('Your version of pytest is too old, you should have '
                      'at least pytest >= {} installed.'
                      .format(PYTEST_MIN_VERSION))


def pytest_addoption(parser):
    parser.addoption("--skip-network", action="store_true", default=False,
                     help="skip network tests")


def pytest_collection_modifyitems(config, items):

    # FeatureHasher is not compatible with PyPy
    if platform.python_implementation() == 'PyPy':
        skip_marker = pytest.mark.skip(
            reason='FeatureHasher is not compatible with PyPy')
        for item in items:
            if item.name in (
                    'sklearn.feature_extraction.hashing.FeatureHasher',
                    'sklearn.feature_extraction.text.HashingVectorizer'):
                item.add_marker(skip_marker)

    # Skip tests which require internet if the flag is provided
    if config.getoption("--skip-network"):
        skip_network = pytest.mark.skip(
            reason="test requires internet connectivity")
        for item in items:
            if "network" in item.keywords:
                item.add_marker(skip_network)

    # numpy changed the str/repr formatting of numpy arrays in 1.14. We want to
    # run doctests only for numpy >= 1.14.
    skip_doctests = False
    try:
        import numpy as np
        if LooseVersion(np.__version__) < LooseVersion('1.14'):
            reason = 'doctests are only run for numpy >= 1.14'
            skip_doctests = True
        elif _IS_32BIT:
            reason = ('doctest are only run when the default numpy int is '
                      '64 bits.')
            skip_doctests = True
    except ImportError:
        pass

    if skip_doctests:
        skip_marker = pytest.mark.skip(reason=reason)

        for item in items:
            if isinstance(item, DoctestItem):
                item.add_marker(skip_marker)


def pytest_configure(config):
    import sys
    sys._is_pytest_session = True


def pytest_unconfigure(config):
    import sys
    del sys._is_pytest_session


def pytest_runtest_setup(item):
    if isinstance(item, DoctestItem):
        set_config(print_changed_only=True)


def pytest_runtest_teardown(item, nextitem):
    if isinstance(item, DoctestItem):
        set_config(print_changed_only=False)
