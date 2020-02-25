"""Test the california_housing loader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""
import pytest

from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


def test_fetch(fetch_california_housing):
    data = fetch_california_housing()
    assert((20640, 8) == data.data.shape)
    assert((20640, ) == data.target.shape)

    # test return_X_y option
    fetch_func = partial(fetch_california_housing)
    check_return_X_y(data, fetch_func)


def test_fetch_asframe(fetch_california_housing):
    pd = pytest.importorskip('pandas')
    bunch = fetch_california_housing(as_frame=True)
    frame = bunch.frame
    assert hasattr(bunch, 'frame') is True
    assert frame.shape == (20640, 9)
    assert isinstance(bunch.data, pd.DataFrame)
    assert isinstance(bunch.target, pd.Series)


def test_pandas_dependency_message(fetch_california_housing):
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
