"""Test loaders for common functionality.
"""


def check_return_X_y(bunch, X_y_tuple):
    assert(isinstance(X_y_tuple, tuple))
    assert(X_y_tuple[0].shape == bunch.data.shape)
    assert(X_y_tuple[1].shape == bunch.target.shape)
