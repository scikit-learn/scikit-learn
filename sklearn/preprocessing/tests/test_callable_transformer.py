from nose.tools import assert_equal
import numpy as np

from ..callable_transformer import CallableTransformer


def _make_func(args_store, kwargs_store, func=lambda X, *a, **k: X):
    def _func(X, *args, **kwargs):
        args_store.append(X)
        args_store.extend(args)
        kwargs_store.update(kwargs)
        return func(X)

    return _func


def test_delegate_to_func():
    # (args|kwargs)_store will hold the positional and keyword arguments
    # passed to the function inside the CallableTransformer.
    args_store = []
    kwargs_store = {}
    X = np.arange(10).reshape((5, 2))
    np.testing.assert_array_equal(
        CallableTransformer(_make_func(args_store, kwargs_store)).transform(X),
        X,
        'transform should have returned X unchanged',
    )

    # The function should only have recieved X and y, where y is None.
    assert_equal(
        args_store,
        [X, None],
        'Incorrect positional arguments passed to func: {args}'.format(
            args=args_store,
        ),
    )
    assert_equal(
        kwargs_store,
        {},
        'Unexpected keyword arguments passed to func: {args}'.format(
            args=kwargs_store,
        ),
    )


def test_argument_closure():
    # (args|kwargs)_store will hold the positional and keyword arguments
    # passed to the function inside the CallableTransformer.
    args_store = []
    kwargs_store = {}
    args = (object(), object())
    kwargs = {'a': object(), 'b': object()}
    X = np.arange(10).reshape((5, 2))
    np.testing.assert_array_equal(
        CallableTransformer(
            _make_func(args_store, kwargs_store),
            args=args,
            kwargs=kwargs,
        ).transform(X),
        X,
        'transform should have returned X unchanged',
    )

    # The function should have been passed X, y (None), and the args
    # that were passed to the CallableTransformer.
    assert_equal(
        args_store,
        [X, None] + list(args),
        'Incorrect positional arguments passed to func: {args}'.format(
            args=args_store,
        ),
    )
    assert_equal(
        kwargs_store,
        kwargs,
        'Incorrect keyword arguments passed to func: {args}'.format(
            args=kwargs_store,
        ),
    )
