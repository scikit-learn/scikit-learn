import random
from distutils.version import LooseVersion

import numpy as np
import pytest

import dask.array as da
from dask.array.utils import assert_eq

sparse = pytest.importorskip('sparse')

if LooseVersion(np.__version__) < '1.11.2':
    pytestmark = pytest.mark.skip


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
    x[x < 0.8] = 0

    y = x.map_blocks(sparse.COO.from_numpy)

    xx = func(x)
    yy = func(y)

    assert_eq(xx, yy)

    if yy.shape:
        zz = yy.compute()
        if not isinstance(zz, sparse.COO):
            assert (zz != 1).sum() > np.prod(zz.shape) / 2  # mostly dense


def test_tensordot():
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    x[x < 0.8] = 0
    y = da.random.random((4, 3, 2), chunks=(2, 2, 1))
    y[y < 0.8] = 0

    xx = x.map_blocks(sparse.COO.from_numpy)
    yy = y.map_blocks(sparse.COO.from_numpy)

    assert_eq(da.tensordot(x, y, axes=(2, 0)),
              da.tensordot(xx, yy, axes=(2, 0)))
    assert_eq(da.tensordot(x, y, axes=(1, 1)),
              da.tensordot(xx, yy, axes=(1, 1)))
    assert_eq(da.tensordot(x, y, axes=((1, 2), (1, 0))),
              da.tensordot(xx, yy, axes=((1, 2), (1, 0))))


@pytest.mark.parametrize('func', functions)
def test_mixed_concatenate(func):
    x = da.random.random((2, 3, 4), chunks=(1, 2, 2))

    y = da.random.random((2, 3, 4), chunks=(1, 2, 2))
    y[y < 0.8] = 0
    yy = y.map_blocks(sparse.COO.from_numpy)

    d = da.concatenate([x, y], axis=0)
    s = da.concatenate([x, yy], axis=0)

    dd = func(d)
    ss = func(s)

    assert_eq(dd, ss)


@pytest.mark.parametrize('func', functions)
def test_mixed_random(func):
    d = da.random.random((4, 3, 4), chunks=(1, 2, 2))
    d[d < 0.7] = 0

    fn = lambda x: sparse.COO.from_numpy(x) if random.random() < 0.5 else x
    s = d.map_blocks(fn)

    dd = func(d)
    ss = func(s)

    assert_eq(dd, ss)


def test_mixed_output_type():
    y = da.random.random((10, 10), chunks=(5, 5))
    y[y < 0.8] = 0
    y = y.map_blocks(sparse.COO.from_numpy)

    x = da.zeros((10, 1), chunks=(5, 1))

    z = da.concatenate([x, y], axis=1)

    assert z.shape == (10, 11)

    zz = z.compute()
    assert isinstance(zz, sparse.COO)
    assert zz.nnz == y.compute().nnz
