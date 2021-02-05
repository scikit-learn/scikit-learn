# Even if empty this file is useful so that when running from the root folder
# ./sklearn is added to sys.path by pytest. See
# https://docs.pytest.org/en/latest/pythonpath.html for more details.  For
# example, this allows to build extensions in place and run pytest
# doc/modules/clustering.rst and use sklearn from the local folder rather than
# the one from site-packages.

import platform
import sys

import pytest
from _pytest.doctest import DoctestItem

from sklearn.utils import _IS_32BIT
from sklearn.externals import _pilutil
from sklearn._min_dependencies import PYTEST_MIN_VERSION
from sklearn.utils.fixes import np_version, parse_version

if parse_version(pytest.__version__) < parse_version(PYTEST_MIN_VERSION):
    raise ImportError('Your version of pytest is too old, you should have '
                      'at least pytest >= {} installed.'
                      .format(PYTEST_MIN_VERSION))


def pytest_collection_modifyitems(config, items):
    for item in items:
        # FeatureHasher is not compatible with PyPy
        if (item.name.endswith(('_hash.FeatureHasher',
                                'text.HashingVectorizer'))
                and platform.python_implementation() == 'PyPy'):
            marker = pytest.mark.skip(
                reason='FeatureHasher is not compatible with PyPy')
            item.add_marker(marker)
        # Known failure on with GradientBoostingClassifier on ARM64
        elif (item.name.endswith('GradientBoostingClassifier')
                and platform.machine() == 'aarch64'):

            marker = pytest.mark.xfail(
                reason=(
                    'know failure. See '
                    'https://github.com/scikit-learn/scikit-learn/issues/17797'  # noqa
                )
            )
            item.add_marker(marker)

    # numpy changed the str/repr formatting of numpy arrays in 1.14. We want to
    # run doctests only for numpy >= 1.14.
    skip_doctests = False
    try:
        if np_version < parse_version('1.14'):
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
