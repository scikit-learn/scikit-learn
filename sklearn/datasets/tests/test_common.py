"""Test loaders for common functionality.
"""
import pytest


def check_return_X_y(bunch, fetch_func_partial):
    X_y_tuple = fetch_func_partial(return_X_y=True)
    assert isinstance(X_y_tuple, tuple)
    assert X_y_tuple[0].shape == bunch.data.shape
    assert X_y_tuple[1].shape == bunch.target.shape


def check_as_frame(bunch, fetch_func_partial):
    pd = pytest.importorskip('pandas')
    frame_bunch = fetch_func_partial(as_frame=True)
    assert hasattr(frame_bunch, 'frame') is True
    assert isinstance(frame_bunch.frame, pd.DataFrame)
    assert isinstance(frame_bunch.data, pd.DataFrame)
    assert frame_bunch.data.shape == bunch.data.shape
    assert isinstance(frame_bunch.target, pd.DataFrame)
    assert frame_bunch.target.shape == bunch.target.shape
