"""Test loaders for common functionality.
"""
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_true


def check_return_X_y(bunch, X_y_tuple):
    assert_true(isinstance(X_y_tuple, tuple))
    assert_array_equal(X_y_tuple[0].shape, bunch.data.shape)
    assert_array_equal(X_y_tuple[1].shape, bunch.target.shape)
