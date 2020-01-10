"""Test the california_housing loader.

Skipped if california_housing is not already downloaded to data_home.
"""

import pytest

from sklearn.datasets import fetch_california_housing
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def fetch(*args, **kwargs):
    return fetch_california_housing(*args, download_if_missing=False, **kwargs)


def _is_california_housing_dataset_not_available():
    try:
        fetch_california_housing(download_if_missing=False)
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
    assert isinstance(bunch.target, pd.DataFrame)


@pytest.mark.skipif(
    _is_california_housing_dataset_not_available(),
    reason='Download California Housing dataset to run this test'
)
def test_pandas_dependency_message():
    try:
        import pandas  # noqa
        pytest.skip("This test requires pandas to be not installed")
    except ImportError:
        # Check that pandas is imported lazily and that an informative error
        # message is raised when pandas is missing:
        expected_msg = ('fetch_california_housing with as_frame=True'
                        ' requires pandas')
        with pytest.raises(ImportError, match=expected_msg):
            fetch_california_housing(as_frame=True)
