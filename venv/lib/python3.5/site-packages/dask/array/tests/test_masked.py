import random
from itertools import product

import numpy as np
import pytest

import dask.array as da
from dask.base import tokenize
from dask.array.utils import assert_eq

pytest.importorskip("dask.array.ma")


def test_tokenize_masked_array():
    m = np.ma.masked_array([1, 2, 3], mask=[True, True, False], fill_value=10)
    m2 = np.ma.masked_array([1, 2, 3], mask=[True, True, False], fill_value=0)
    m3 = np.ma.masked_array([1, 2, 3], mask=False, fill_value=10)
    assert tokenize(m) == tokenize(m)
    assert tokenize(m2) == tokenize(m2)
    assert tokenize(m3) == tokenize(m3)
    assert tokenize(m) != tokenize(m2)
    assert tokenize(m) != tokenize(m3)


def test_from_array_masked_array():
    m = np.ma.masked_array([1, 2, 3], mask=[True, True, False], fill_value=10)
    dm = da.from_array(m, chunks=(2,), asarray=False)
    assert_eq(dm, m)


functions = [
    lambda x: x,
    lambda x: da.expm1(x),
    lambda x: 2 * x,
    lambda x: x / 2,
    lambda x: x**2,
    lambda x: x + x,
    lambda x: x * x,
    lambda x: x[0],
    lambda x: x[:, 1],
    lambda x: x[:1, None, 1:3],
    lambda x: x.T,
    lambda x: da.transpose(x, (1, 2, 0)),
    lambda x: x.sum(),
    lambda x: x.dot(np.arange(x.shape[-1])),
    lambda x: x.dot(np.eye(x.shape[-1])),
    lambda x: da.tensordot(x, np.ones(x.shape[:2]), axes=[(0, 1), (0, 1)]),
    lambda x: x.sum(axis=0),
    lambda x: x.max(axis=0),
    lambda x: x.sum(axis=(1, 2)),
    lambda x: x.astype(np.complex128),
    lambda x: x.map_blocks(lambda x: x * 2),
    lambda x: x.round(1),
    lambda x: x.reshape((x.shape[0] * x.shape[1], x.shape[2])),
    lambda x: abs(x),
    lambda x: x > 0.5,
    lambda x: x.rechunk((4, 4, 4)),
    lambda x: x.rechunk((2, 2, 1)),
]


@pytest.mark.parametrize('func', functions)
def test_basic(func):
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.4] = 0

    y = da.ma.masked_equal(x, 0)

    xx = func(x)
    yy = func(y)

    assert_eq(xx, da.ma.filled(yy, 0))

    if yy.shape:
        zz = yy.compute()
        assert isinstance(zz, np.ma.masked_array)


def test_tensordot():
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.4] = 0
    y = da.random.random((4, 3, 2), chunks=(2, 2, 1))
    y[y < 0.4] = 0

    xx = da.ma.masked_equal(x, 0)
    yy = da.ma.masked_equal(y, 0)

    assert_eq(da.tensordot(x, y, axes=(2, 0)),
              da.ma.filled(da.tensordot(xx, yy, axes=(2, 0)), 0))
    assert_eq(da.tensordot(x, y, axes=(1, 1)),
              da.ma.filled(da.tensordot(xx, yy, axes=(1, 1)), 0))
    assert_eq(da.tensordot(x, y, axes=((1, 2), (1, 0))),
              da.ma.filled(da.tensordot(xx, yy, axes=((1, 2), (1, 0))), 0))


@pytest.mark.parametrize('func', functions)
def test_mixed_concatenate(func):
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    y = da.random.random((2, 3, 4), chunks=(1, 2, 2))

    y[y < 0.4] = 0
    yy = da.ma.masked_equal(y, 0)

    d = da.concatenate([x, y], axis=0)
    s = da.concatenate([x, yy], axis=0)

    dd = func(d)
    ss = func(s)
    assert_eq(dd, ss)


@pytest.mark.parametrize('func', functions)
def test_mixed_random(func):
    d = da.random.random((4, 3, 4), chunks=(1, 2, 2))
    d[d < 0.4] = 0

    fn = lambda x: np.ma.masked_equal(x, 0) if random.random() < 0.5 else x
    s = d.map_blocks(fn)

    dd = func(d)
    ss = func(s)

    assert_eq(dd, ss)


def test_mixed_output_type():
    y = da.random.random((10, 10), chunks=(5, 5))
    y[y < 0.4] = 0

    y = da.ma.masked_equal(y, 0)
    x = da.zeros((10, 1), chunks=(5, 1))

    z = da.concatenate([x, y], axis=1)
    assert z.shape == (10, 11)
    zz = z.compute()
    assert isinstance(zz, np.ma.masked_array)


def test_creation_functions():
    x = np.array([-2, -1, 0, 1, 2] * 20).reshape((10, 10))
    y = np.array([-2, 0, 1, 1, 0] * 2)
    dx = da.from_array(x, chunks=5)
    dy = da.from_array(y, chunks=4)

    sol = np.ma.masked_greater(x, y)
    for (a, b) in product([dx, x], [dy, y]):
        assert_eq(da.ma.masked_greater(a, b), sol)

    # These are all the same as masked_greater, just check for correct op
    assert_eq(da.ma.masked_greater(dx, 0), np.ma.masked_greater(x, 0))
    assert_eq(da.ma.masked_greater_equal(dx, 0), np.ma.masked_greater_equal(x, 0))
    assert_eq(da.ma.masked_less(dx, 0), np.ma.masked_less(x, 0))
    assert_eq(da.ma.masked_less_equal(dx, 0), np.ma.masked_less_equal(x, 0))
    assert_eq(da.ma.masked_equal(dx, 0), np.ma.masked_equal(x, 0))
    assert_eq(da.ma.masked_not_equal(dx, 0), np.ma.masked_not_equal(x, 0))

    # masked_where
    assert_eq(da.ma.masked_where(False, dx), np.ma.masked_where(False, x))
    assert_eq(da.ma.masked_where(dx > 2, dx), np.ma.masked_where(x > 2, x))

    with pytest.raises(IndexError):
        da.ma.masked_where((dx > 2)[:, 0], dx)

    assert_eq(da.ma.masked_inside(dx, -1, 1), np.ma.masked_inside(x, -1, 1))
    assert_eq(da.ma.masked_outside(dx, -1, 1), np.ma.masked_outside(x, -1, 1))
    assert_eq(da.ma.masked_values(dx, -1), np.ma.masked_values(x, -1))

    # masked_equal and masked_values in numpy sets the fill_value to `value`,
    # which can sometimes be an array. This is hard to support in dask, so we
    # forbid it. Check that this isn't supported:
    with pytest.raises(ValueError):
        da.ma.masked_equal(dx, dy)

    with pytest.raises(ValueError):
        da.ma.masked_values(dx, dy)

    y = x.astype('f8')
    y[0, 0] = y[7, 5] = np.nan
    dy = da.from_array(y, chunks=5)

    assert_eq(da.ma.masked_invalid(dy), np.ma.masked_invalid(y))

    my = np.ma.masked_greater(y, 0)
    dmy = da.ma.masked_greater(dy, 0)

    assert_eq(da.ma.fix_invalid(dmy, fill_value=0),
              np.ma.fix_invalid(my, fill_value=0))


def test_filled():
    x = np.array([-2, -1, 0, 1, 2] * 20).reshape((10, 10))
    dx = da.from_array(x, chunks=5)

    mx = np.ma.masked_equal(x, 0)
    mdx = da.ma.masked_equal(dx, 0)

    assert_eq(da.ma.filled(mdx), np.ma.filled(mx))
    assert_eq(da.ma.filled(mdx, -5), np.ma.filled(mx, -5))


def assert_eq_ma(a, b):
    res = a.compute()
    assert type(res) == type(b)
    if hasattr(res, 'mask'):
        np.testing.assert_equal(res.mask, b.mask)
        a = da.ma.filled(a)
        b = np.ma.filled(b)
    assert_eq(a, b, equal_nan=True)


@pytest.mark.parametrize('dtype', ('i8', 'f8'))
@pytest.mark.parametrize('reduction', ['sum', 'prod', 'mean', 'var', 'std',
                                       'min', 'max', 'any', 'all'])
def test_reductions(dtype, reduction):
    x = (np.random.RandomState(42).rand(11, 11) * 10).astype(dtype)
    dx = da.from_array(x, chunks=(4, 4))
    mx = np.ma.masked_greater(x, 5)
    mdx = da.ma.masked_greater(dx, 5)

    dfunc = getattr(da, reduction)
    func = getattr(np, reduction)

    assert_eq_ma(dfunc(mdx), func(mx))
    assert_eq_ma(dfunc(mdx, axis=0), func(mx, axis=0))
    assert_eq_ma(dfunc(mdx, keepdims=True, split_every=4),
                 func(mx, keepdims=True))
    assert_eq_ma(dfunc(mdx, axis=0, split_every=2), func(mx, axis=0))
    assert_eq_ma(dfunc(mdx, axis=0, keepdims=True, split_every=2),
                 func(mx, axis=0, keepdims=True))
    assert_eq_ma(dfunc(mdx, axis=1, split_every=2), func(mx, axis=1))
    assert_eq_ma(dfunc(mdx, axis=1, keepdims=True, split_every=2),
                 func(mx, axis=1, keepdims=True))


@pytest.mark.parametrize('reduction', ['argmin', 'argmax'])
def test_arg_reductions(reduction):
    x = np.random.random((10, 10, 10))
    dx = da.from_array(x, chunks=(3, 4, 5))
    mx = np.ma.masked_greater(x, 0.4)
    dmx = da.ma.masked_greater(dx, 0.4)

    dfunc = getattr(da, reduction)
    func = getattr(np, reduction)

    assert_eq_ma(dfunc(dmx), func(mx))
    assert_eq_ma(dfunc(dmx, 0), func(mx, 0))
    assert_eq_ma(dfunc(dmx, 1), func(mx, 1))
    assert_eq_ma(dfunc(dmx, 2), func(mx, 2))


def test_cumulative():
    x = np.random.RandomState(0).rand(20, 24, 13)
    dx = da.from_array(x, chunks=(6, 5, 4))
    mx = np.ma.masked_greater(x, 0.4)
    dmx = da.ma.masked_greater(dx, 0.4)

    for axis in [0, 1, 2]:
        assert_eq_ma(dmx.cumsum(axis=axis), mx.cumsum(axis=axis))
        assert_eq_ma(dmx.cumprod(axis=axis), mx.cumprod(axis=axis))


def test_accessors():
    x = np.random.random((10, 10))
    dx = da.from_array(x, chunks=(3, 4))
    mx = np.ma.masked_greater(x, 0.4)
    dmx = da.ma.masked_greater(dx, 0.4)

    assert_eq(da.ma.getmaskarray(dmx), np.ma.getmaskarray(mx))
    assert_eq(da.ma.getmaskarray(dx), np.ma.getmaskarray(x))
    assert_eq(da.ma.getdata(dmx), np.ma.getdata(mx))
    assert_eq(da.ma.getdata(dx), np.ma.getdata(x))


def test_masked_array():
    x = np.random.random((10, 10)).astype('f4')
    dx = da.from_array(x, chunks=(3, 4))
    f1 = da.from_array(np.array(1), chunks=())

    fill_values = [(None, None), (0.5, 0.5), (1, f1)]
    for data, (df, f) in product([x, dx], fill_values):
        assert_eq(da.ma.masked_array(data, fill_value=df),
                  np.ma.masked_array(x, fill_value=f))
        assert_eq(da.ma.masked_array(data, mask=data > 0.4, fill_value=df),
                  np.ma.masked_array(x, mask=x > 0.4, fill_value=f))
        assert_eq(da.ma.masked_array(data, mask=data > 0.4, fill_value=df),
                  np.ma.masked_array(x, mask=x > 0.4, fill_value=f))
        assert_eq(da.ma.masked_array(data, fill_value=df, dtype='f8'),
                  np.ma.masked_array(x, fill_value=f, dtype='f8'))

    with pytest.raises(ValueError):
        da.ma.masked_array(dx, fill_value=dx)

    with pytest.raises(np.ma.MaskError):
        da.ma.masked_array(dx, mask=dx[:3, :3])


def test_set_fill_value():
    x = np.random.randint(0, 10, (10, 10))
    dx = da.from_array(x, chunks=(3, 4))
    mx = np.ma.masked_greater(x, 3)
    dmx = da.ma.masked_greater(dx, 3)

    da.ma.set_fill_value(dmx, -10)
    np.ma.set_fill_value(mx, -10)
    assert_eq_ma(dmx, mx)

    da.ma.set_fill_value(dx, -10)
    np.ma.set_fill_value(x, -10)
    assert_eq_ma(dx, x)

    with pytest.raises(TypeError):
        da.ma.set_fill_value(dmx, 1e20)

    with pytest.raises(ValueError):
        da.ma.set_fill_value(dmx, dx)
