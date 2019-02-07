"""Test loaders for common functionality.
"""


def check_return_X_y(bunch, fetch_func_partial):
    X_y_tuple = fetch_func_partial(return_X_y=True)
    assert(isinstance(X_y_tuple, tuple))
    assert(X_y_tuple[0].shape == bunch.data.shape)
    assert(X_y_tuple[1].shape == bunch.target.shape)
