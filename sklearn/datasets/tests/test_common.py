"""Test loaders for common functionality.
"""
import pytest
import numpy as np


def check_pandas_dependency_message(fetch_func):
    try:
        import pandas  # noqa
        pytest.skip("This test requires pandas to be not installed")
    except ImportError:
        # Check that pandas is imported lazily and that an informative error
        # message is raised when pandas is missing:
        expected_msg = ('{} with as_frame=True requires pandas'
                        .format(fetch_func.__name__))
        with pytest.raises(ImportError, match=expected_msg):
            fetch_func(as_frame=True)


def check_return_X_y(bunch, fetch_func_partial):
    X_y_tuple = fetch_func_partial(return_X_y=True)
    assert isinstance(X_y_tuple, tuple)
    assert X_y_tuple[0].shape == bunch.data.shape
    assert X_y_tuple[1].shape == bunch.target.shape


def check_as_frame(bunch, fetch_func_partial,
                   expected_data_dtype=None, expected_target_dtype=None):
    pd = pytest.importorskip('pandas')
    frame_bunch = fetch_func_partial(as_frame=True)
    assert hasattr(frame_bunch, 'frame')
    assert isinstance(frame_bunch.frame, pd.DataFrame)
    assert isinstance(frame_bunch.data, pd.DataFrame)
    assert frame_bunch.data.shape == bunch.data.shape
    if frame_bunch.target.ndim > 1:
        assert isinstance(frame_bunch.target, pd.DataFrame)
    else:
        assert isinstance(frame_bunch.target, pd.Series)
    assert frame_bunch.target.shape[0] == bunch.target.shape[0]
    if expected_data_dtype is not None:
        assert np.all(frame_bunch.data.dtypes == expected_data_dtype)
    if expected_target_dtype is not None:
        assert np.all(frame_bunch.target.dtypes == expected_target_dtype)
