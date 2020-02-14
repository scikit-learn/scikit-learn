"""Test the california_housing loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""

import builtins
from os import environ
import pytest

from sklearn.datasets import fetch_california_housing
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def fetch(*args, **kwargs):
    return fetch_california_housing(*args, **kwargs)


def _is_california_housing_dataset_not_available():
    # Do not download data, unless explicitly requested via environment var
    download_if_missing = environ.get('SKLEARN_SKIP_NETWORK_TESTS', '1') == '0'
    try:
        fetch(download_if_missing=download_if_missing)
        return False
    except IOError:
        return True


@pytest.mark.skipif(
    _is_california_housing_dataset_not_available(),
    reason='Download California Housing dataset to run this test'
)
def test_fetch():
    data = fetch()
    assert((20640, 8) == data.data.shape)
    assert((20640, ) == data.target.shape)

    # test return_X_y option
    fetch_func = partial(fetch)
    check_return_X_y(data, fetch_func)


@pytest.mark.skipif(
    _is_california_housing_dataset_not_available(),
    reason='Download California Housing dataset to run this test'
)
def test_fetch_asframe():
    pd = pytest.importorskip('pandas')
    bunch = fetch(as_frame=True)
    frame = bunch.frame
    assert hasattr(bunch, 'frame') is True
    assert frame.shape == (20640, 9)
    assert isinstance(bunch.data, pd.DataFrame)
    assert isinstance(bunch.target, pd.Series)


@pytest.fixture
def hide_available_pandas(monkeypatch):
    """ Pretend pandas was not installed. """
    import_orig = builtins.__import__

    def mocked_import(name, *args, **kwargs):
        if name == 'pandas':
            raise ImportError()
        return import_orig(name, *args, **kwargs)

    monkeypatch.setattr(builtins, '__import__', mocked_import)


@pytest.mark.skipif(
    _is_california_housing_dataset_not_available(),
    reason='Download California Housing dataset to run this test'
)
@pytest.mark.usefixtures('hide_available_pandas')
def test_pandas_dependency_message():
    # Check that pandas is imported lazily and that an informative error
    # message is raised when pandas is missing:
    expected_msg = ('fetch_california_housing with as_frame=True'
                    ' requires pandas')
    with pytest.raises(ImportError, match=expected_msg):
        fetch_california_housing(as_frame=True)
