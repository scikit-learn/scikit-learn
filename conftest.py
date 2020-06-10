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

from sklearn.utils import _IS_32BIT
from sklearn.externals import _pilutil
from sklearn.datasets import fetch_20newsgroups
from sklearn.datasets import fetch_20newsgroups_vectorized
from sklearn.datasets import fetch_california_housing
from sklearn.datasets import fetch_covtype
from sklearn.datasets import fetch_kddcup99
from sklearn.datasets import fetch_olivetti_faces
from sklearn.datasets import fetch_rcv1


PYTEST_MIN_VERSION = '3.3.0'

if LooseVersion(pytest.__version__) < PYTEST_MIN_VERSION:
    raise ImportError('Your version of pytest is too old, you should have '
                      'at least pytest >= {} installed.'
                      .format(PYTEST_MIN_VERSION))


dataset_fetchers = {
    'fetch_20newsgroups_fxt': fetch_20newsgroups,
    'fetch_20newsgroups_vectorized_fxt': fetch_20newsgroups_vectorized,
    'fetch_california_housing_fxt': fetch_california_housing,
    'fetch_covtype_fxt': fetch_covtype,
    'fetch_kddcup99_fxt': fetch_kddcup99,
    'fetch_olivetti_faces_fxt': fetch_olivetti_faces,
    'fetch_rcv1_fxt': fetch_rcv1,
}


# fetching a dataset with this fixture will never download if missing
def _fetch_fixture(f):
    def wrapped(*args, **kwargs):
        kwargs['download_if_missing'] = False
        try:
            return f(*args, **kwargs)
        except IOError:
            pytest.skip("test requires -m 'not skipnetwork' to run")
    return pytest.fixture(lambda: wrapped)


fetch_20newsgroups_fxt = _fetch_fixture(fetch_20newsgroups)
fetch_20newsgroups_vectorized_fxt = \
    _fetch_fixture(fetch_20newsgroups_vectorized)
fetch_california_housing_fxt = _fetch_fixture(fetch_california_housing)
fetch_covtype_fxt = _fetch_fixture(fetch_covtype)
fetch_kddcup99_fxt = _fetch_fixture(fetch_kddcup99)
fetch_olivetti_faces_fxt = _fetch_fixture(fetch_olivetti_faces)
fetch_rcv1_fxt = _fetch_fixture(fetch_rcv1)


def pytest_collection_modifyitems(config, items):

    # FeatureHasher is not compatible with PyPy
    if platform.python_implementation() == 'PyPy':
        skip_marker = pytest.mark.skip(
            reason='FeatureHasher is not compatible with PyPy')
        for item in items:
            if item.name.endswith(('_hash.FeatureHasher',
                                   'text.HashingVectorizer')):
                item.add_marker(skip_marker)

    run_network_tests = 'not skipnetwork' in config.getoption("markexpr")
    skip_network = pytest.mark.skip(
        reason="test requires -m 'not skipnetwork' to run")

    # download datasets during collection to avoid thread unsafe behavior
    # when running pytest in parallel with pytest-xdist
    dataset_features_set = set(dataset_fetchers)
    datasets_to_download = set()

    for item in items:
        item_keywords = set(item.keywords)
        dataset_to_fetch = item_keywords & dataset_features_set
        if not dataset_to_fetch:
            continue

        if run_network_tests:
            datasets_to_download |= dataset_to_fetch
        else:
            # network tests are skipped
            item.add_marker(skip_network)

    # download datasets that are needed to avoid thread unsafe behavior
    # by pytest-xdist
    if run_network_tests:
        for name in datasets_to_download:
            dataset_fetchers[name]()

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
        elif sys.platform.startswith("win32"):
            reason = ("doctests are not run for Windows because numpy arrays "
                      "repr is inconsistent across platforms.")
            skip_doctests = True
    except ImportError:
        pass

    if skip_doctests:
        skip_marker = pytest.mark.skip(reason=reason)

        for item in items:
            if isinstance(item, DoctestItem):
                item.add_marker(skip_marker)
    elif not _pilutil.pillow_installed:
        skip_marker = pytest.mark.skip(reason="pillow (or PIL) not installed!")
        for item in items:
            if item.name in [
                    "sklearn.feature_extraction.image.PatchExtractor",
                    "sklearn.feature_extraction.image.extract_patches_2d"]:
                item.add_marker(skip_marker)


def pytest_runtest_setup(item):
    run_network_tests = 'not skipnetwork' in item.config.getoption("markexpr")
    if "network" in item.keywords and not run_network_tests:
        pytest.skip("test requires -m 'not skipnetwork' to run")


def pytest_configure(config):
    import sys
    sys._is_pytest_session = True
    # declare our custom markers to avoid PytestUnknownMarkWarning
    config.addinivalue_line(
        "markers",
        "network: mark a test for execution if network available."
    )


def pytest_unconfigure(config):
    import sys
    del sys._is_pytest_session
