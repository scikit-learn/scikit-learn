from __future__ import absolute_import, division, print_function

import pytest
np = pytest.importorskip('numpy')

import dask.array as da
from dask.array.utils import assert_eq as _assert_eq, same_keys
from dask.core import get_deps
from dask.context import set_options


def assert_eq(a, b):
    _assert_eq(a, b, equal_nan=True)


def reduction_1d_test(da_func, darr, np_func, narr, use_dtype=True, split_every=True):
    assert_eq(da_func(darr), np_func(narr))
    assert_eq(da_func(darr, keepdims=True), np_func(narr, keepdims=True))
    assert same_keys(da_func(darr), da_func(darr))
    assert same_keys(da_func(darr, keepdims=True), da_func(darr, keepdims=True))
    if use_dtype:
        assert_eq(da_func(darr, dtype='f8'), np_func(narr, dtype='f8'))
        assert_eq(da_func(darr, dtype='i8'), np_func(narr, dtype='i8'))
        assert same_keys(da_func(darr, dtype='i8'), da_func(darr, dtype='i8'))
    if split_every:
        a1 = da_func(darr, split_every=2)
        a2 = da_func(darr, split_every={0: 2})
        assert same_keys(a1, a2)
        assert_eq(a1, np_func(narr))
        assert_eq(a2, np_func(narr))
        assert_eq(da_func(darr, keepdims=True, split_every=2),
                  np_func(narr, keepdims=True))


@pytest.mark.parametrize('dtype', ['f4', 'i4'])
def test_reductions_1D(dtype):
    x = np.arange(5).astype(dtype)
    a = da.from_array(x, chunks=(2,))

    reduction_1d_test(da.sum, a, np.sum, x)
    reduction_1d_test(da.prod, a, np.prod, x)
    reduction_1d_test(da.mean, a, np.mean, x)
    reduction_1d_test(da.var, a, np.var, x)
    reduction_1d_test(da.std, a, np.std, x)
    reduction_1d_test(da.min, a, np.min, x, False)
    reduction_1d_test(da.max, a, np.max, x, False)
    reduction_1d_test(da.any, a, np.any, x, False)
    reduction_1d_test(da.all, a, np.all, x, False)

    reduction_1d_test(da.nansum, a, np.nansum, x)
    reduction_1d_test(da.nanprod, a, np.nanprod, x)
    reduction_1d_test(da.nanmean, a, np.mean, x)
    reduction_1d_test(da.nanvar, a, np.var, x)
    reduction_1d_test(da.nanstd, a, np.std, x)
    reduction_1d_test(da.nanmin, a, np.nanmin, x, False)
    reduction_1d_test(da.nanmax, a, np.nanmax, x, False)


def reduction_2d_test(da_func, darr, np_func, narr, use_dtype=True,
                      split_every=True):
    assert_eq(da_func(darr), np_func(narr))
    assert_eq(da_func(darr, keepdims=True), np_func(narr, keepdims=True))
    assert_eq(da_func(darr, axis=0), np_func(narr, axis=0))
    assert_eq(da_func(darr, axis=1), np_func(narr, axis=1))
    assert_eq(da_func(darr, axis=-1), np_func(narr, axis=-1))
    assert_eq(da_func(darr, axis=-2), np_func(narr, axis=-2))
    assert_eq(da_func(darr, axis=1, keepdims=True),
              np_func(narr, axis=1, keepdims=True))
    assert_eq(da_func(darr, axis=(1, 0)), np_func(narr, axis=(1, 0)))

    assert same_keys(da_func(darr, axis=1), da_func(darr, axis=1))
    assert same_keys(da_func(darr, axis=(1, 0)), da_func(darr, axis=(1, 0)))

    if use_dtype:
        assert_eq(da_func(darr, dtype='f8'), np_func(narr, dtype='f8'))
        assert_eq(da_func(darr, dtype='i8'), np_func(narr, dtype='i8'))

    if split_every:
        a1 = da_func(darr, split_every=4)
        a2 = da_func(darr, split_every={0: 2, 1: 2})
        assert same_keys(a1, a2)
        assert_eq(a1, np_func(narr))
        assert_eq(a2, np_func(narr))
        assert_eq(da_func(darr, keepdims=True, split_every=4),
                  np_func(narr, keepdims=True))
        assert_eq(da_func(darr, axis=0, split_every=2), np_func(narr, axis=0))
        assert_eq(da_func(darr, axis=0, keepdims=True, split_every=2),
                  np_func(narr, axis=0, keepdims=True))
        assert_eq(da_func(darr, axis=1, split_every=2), np_func(narr, axis=1))
        assert_eq(da_func(darr, axis=1, keepdims=True, split_every=2),
                  np_func(narr, axis=1, keepdims=True))


def test_reduction_errors():
    x = da.ones((5, 5), chunks=(3, 3))
    with pytest.raises(ValueError):
        x.sum(axis=2)
    with pytest.raises(ValueError):
        x.sum(axis=-3)


@pytest.mark.slow
@pytest.mark.parametrize('dtype', ['f4', 'i4'])
def test_reductions_2D(dtype):
    x = np.arange(1, 122).reshape((11, 11)).astype(dtype)
    a = da.from_array(x, chunks=(4, 4))

    b = a.sum(keepdims=True)
    assert b.__dask_keys__() == [[(b.name, 0, 0)]]

    reduction_2d_test(da.sum, a, np.sum, x)
    reduction_2d_test(da.prod, a, np.prod, x)
    reduction_2d_test(da.mean, a, np.mean, x)
    reduction_2d_test(da.var, a, np.var, x, False)  # Difference in dtype algo
    reduction_2d_test(da.std, a, np.std, x, False)  # Difference in dtype algo
    reduction_2d_test(da.min, a, np.min, x, False)
    reduction_2d_test(da.max, a, np.max, x, False)
    reduction_2d_test(da.any, a, np.any, x, False)
    reduction_2d_test(da.all, a, np.all, x, False)

    reduction_2d_test(da.nansum, a, np.nansum, x)
    reduction_2d_test(da.nanprod, a, np.nanprod, x)
    reduction_2d_test(da.nanmean, a, np.mean, x)
    reduction_2d_test(da.nanvar, a, np.nanvar, x, False)  # Difference in dtype algo
    reduction_2d_test(da.nanstd, a, np.nanstd, x, False)  # Difference in dtype algo
    reduction_2d_test(da.nanmin, a, np.nanmin, x, False)
    reduction_2d_test(da.nanmax, a, np.nanmax, x, False)


@pytest.mark.parametrize(['dfunc', 'func'],
                         [(da.argmin, np.argmin), (da.argmax, np.argmax),
                          (da.nanargmin, np.nanargmin),
                          (da.nanargmax, np.nanargmax)])
def test_arg_reductions(dfunc, func):
    x = np.random.random((10, 10, 10))
    a = da.from_array(x, chunks=(3, 4, 5))

    assert_eq(dfunc(a), func(x))
    assert_eq(dfunc(a, 0), func(x, 0))
    assert_eq(dfunc(a, 1), func(x, 1))
    assert_eq(dfunc(a, 2), func(x, 2))
    with set_options(split_every=2):
        assert_eq(dfunc(a), func(x))
        assert_eq(dfunc(a, 0), func(x, 0))
        assert_eq(dfunc(a, 1), func(x, 1))
        assert_eq(dfunc(a, 2), func(x, 2))

    pytest.raises(ValueError, lambda: dfunc(a, 3))
    pytest.raises(TypeError, lambda: dfunc(a, (0, 1)))

    x2 = np.arange(10)
    a2 = da.from_array(x2, chunks=3)
    assert_eq(dfunc(a2), func(x2))
    assert_eq(dfunc(a2, 0), func(x2, 0))
    assert_eq(dfunc(a2, 0, split_every=2), func(x2, 0))


@pytest.mark.parametrize(['dfunc', 'func'],
                         [(da.nanargmin, np.nanargmin),
                          (da.nanargmax, np.nanargmax)])
def test_nanarg_reductions(dfunc, func):

    x = np.random.random((10, 10, 10))
    x[5] = np.nan
    a = da.from_array(x, chunks=(3, 4, 5))
    assert_eq(dfunc(a), func(x))
    assert_eq(dfunc(a, 0), func(x, 0))
    with pytest.raises(ValueError):
        with pytest.warns(None):  # All NaN axis
            dfunc(a, 1).compute()

    with pytest.raises(ValueError):
        with pytest.warns(None):  # All NaN axis
            dfunc(a, 2).compute()

    x[:] = np.nan
    a = da.from_array(x, chunks=(3, 4, 5))
    with pytest.raises(ValueError):
        with pytest.warns(None):  # All NaN axis
            dfunc(a).compute()


def test_reductions_2D_nans():
    # chunks are a mix of some/all/no NaNs
    x = np.full((4, 4), np.nan)
    x[:2, :2] = np.array([[1, 2], [3, 4]])
    x[2, 2] = 5
    x[3, 3] = 6
    a = da.from_array(x, chunks=(2, 2))

    reduction_2d_test(da.sum, a, np.sum, x, False, False)
    reduction_2d_test(da.prod, a, np.prod, x, False, False)
    reduction_2d_test(da.mean, a, np.mean, x, False, False)
    reduction_2d_test(da.var, a, np.var, x, False, False)
    reduction_2d_test(da.std, a, np.std, x, False, False)
    reduction_2d_test(da.min, a, np.min, x, False, False)
    reduction_2d_test(da.max, a, np.max, x, False, False)
    reduction_2d_test(da.any, a, np.any, x, False, False)
    reduction_2d_test(da.all, a, np.all, x, False, False)

    reduction_2d_test(da.nansum, a, np.nansum, x, False, False)
    reduction_2d_test(da.nanprod, a, np.nanprod, x, False, False)
    reduction_2d_test(da.nanmean, a, np.nanmean, x, False, False)
    with pytest.warns(None):  # division by 0 warning
        reduction_2d_test(da.nanvar, a, np.nanvar, x, False, False)
    with pytest.warns(None):  # division by 0 warning
        reduction_2d_test(da.nanstd, a, np.nanstd, x, False, False)
    with pytest.warns(None):  # all NaN axis warning
        reduction_2d_test(da.nanmin, a, np.nanmin, x, False, False)
    with pytest.warns(None):  # all NaN axis warning
        reduction_2d_test(da.nanmax, a, np.nanmax, x, False, False)

    assert_eq(da.argmax(a), np.argmax(x))
    assert_eq(da.argmin(a), np.argmin(x))
    with pytest.warns(None):  # all NaN axis warning
        assert_eq(da.nanargmax(a), np.nanargmax(x))
    with pytest.warns(None):  # all NaN axis warning
        assert_eq(da.nanargmin(a), np.nanargmin(x))
    assert_eq(da.argmax(a, axis=0), np.argmax(x, axis=0))
    assert_eq(da.argmin(a, axis=0), np.argmin(x, axis=0))
    with pytest.warns(None):  # all NaN axis warning
        assert_eq(da.nanargmax(a, axis=0), np.nanargmax(x, axis=0))
    with pytest.warns(None):  # all NaN axis warning
        assert_eq(da.nanargmin(a, axis=0), np.nanargmin(x, axis=0))
    assert_eq(da.argmax(a, axis=1), np.argmax(x, axis=1))
    assert_eq(da.argmin(a, axis=1), np.argmin(x, axis=1))
    with pytest.warns(None):  # all NaN axis warning
        assert_eq(da.nanargmax(a, axis=1), np.nanargmax(x, axis=1))
    with pytest.warns(None):  # all NaN axis warning
        assert_eq(da.nanargmin(a, axis=1), np.nanargmin(x, axis=1))


def test_moment():
    def moment(x, n, axis=None):
        return (((x - x.mean(axis=axis, keepdims=True)) ** n).sum(axis=axis) /
                np.ones_like(x).sum(axis=axis))

    # Poorly conditioned
    x = np.array([1., 2., 3.] * 10).reshape((3, 10)) + 1e8
    a = da.from_array(x, chunks=5)
    assert_eq(a.moment(2), moment(x, 2))
    assert_eq(a.moment(3), moment(x, 3))
    assert_eq(a.moment(4), moment(x, 4))

    x = np.arange(1, 122).reshape((11, 11)).astype('f8')
    a = da.from_array(x, chunks=(4, 4))
    assert_eq(a.moment(4, axis=1), moment(x, 4, axis=1))
    assert_eq(a.moment(4, axis=(1, 0)), moment(x, 4, axis=(1, 0)))

    # Tree reduction
    assert_eq(a.moment(order=4, split_every=4), moment(x, 4))
    assert_eq(a.moment(order=4, axis=0, split_every=4), moment(x, 4, axis=0))
    assert_eq(a.moment(order=4, axis=1, split_every=4), moment(x, 4, axis=1))


def test_reductions_with_negative_axes():
    x = np.random.random((4, 4, 4))
    a = da.from_array(x, chunks=2)

    assert_eq(a.argmin(axis=-1), x.argmin(axis=-1))
    assert_eq(a.argmin(axis=-1, split_every=2), x.argmin(axis=-1))

    assert_eq(a.sum(axis=-1), x.sum(axis=-1))
    assert_eq(a.sum(axis=(0, -1)), x.sum(axis=(0, -1)))


def test_nan():
    x = np.array([[1, np.nan, 3, 4],
                  [5, 6, 7, np.nan],
                  [9, 10, 11, 12]])
    d = da.from_array(x, chunks=(2, 2))

    assert_eq(np.nansum(x), da.nansum(d))
    assert_eq(np.nansum(x, axis=0), da.nansum(d, axis=0))
    assert_eq(np.nanmean(x, axis=1), da.nanmean(d, axis=1))
    assert_eq(np.nanmin(x, axis=1), da.nanmin(d, axis=1))
    assert_eq(np.nanmax(x, axis=(0, 1)), da.nanmax(d, axis=(0, 1)))
    assert_eq(np.nanvar(x), da.nanvar(d))
    assert_eq(np.nanstd(x, axis=0), da.nanstd(d, axis=0))
    assert_eq(np.nanargmin(x, axis=0), da.nanargmin(d, axis=0))
    assert_eq(np.nanargmax(x, axis=0), da.nanargmax(d, axis=0))
    assert_eq(np.nanprod(x), da.nanprod(d))


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='nanmax/nanmin for object dtype')
@pytest.mark.parametrize('func', ['nansum', 'sum', 'nanmin', 'min',
                                  'nanmax', 'max'])
def test_nan_object(func):
    x = np.array([[1, np.nan, 3, 4],
                  [5, 6, 7, np.nan],
                  [9, 10, 11, 12]]).astype(object)
    d = da.from_array(x, chunks=(2, 2))

    assert_eq(getattr(np, func)(x, axis=0), getattr(da, func)(d, axis=0))
    assert_eq(getattr(np, func)(x, axis=1), getattr(da, func)(d, axis=1))
    assert_eq(getattr(np, func)(x), getattr(da, func)(d))


def test_0d_array():
    x = da.mean(da.ones(4, chunks=4), axis=0).compute()
    y = np.mean(np.ones(4))
    assert type(x) == type(y)

    x = da.sum(da.zeros(4, chunks=1)).compute()
    y = np.sum(np.zeros(4))
    assert type(x) == type(y)


def test_reduction_on_scalar():
    x = da.from_array(np.array(1.0), chunks=())
    assert (x == x).all()


def test_reductions_with_empty_array():
    dx1 = da.ones((10, 0, 5), chunks=4)
    x1 = dx1.compute()
    dx2 = da.ones((0, 0, 0), chunks=4)
    x2 = dx2.compute()

    for dx, x in [(dx1, x1), (dx2, x2)]:
        with pytest.warns(None):  # empty slice warning
            assert_eq(dx.mean(), x.mean())
            assert_eq(dx.mean(axis=0), x.mean(axis=0))
            assert_eq(dx.mean(axis=1), x.mean(axis=1))
            assert_eq(dx.mean(axis=2), x.mean(axis=2))


def assert_max_deps(x, n, eq=True):
    dependencies, dependents = get_deps(x.dask)
    if eq:
        assert max(map(len, dependencies.values())) == n
    else:
        assert max(map(len, dependencies.values())) <= n


def test_tree_reduce_depth():
    # 2D
    x = da.from_array(np.arange(242).reshape((11, 22)), chunks=(3, 4))
    thresh = {0: 2, 1: 3}
    assert_max_deps(x.sum(split_every=thresh), 2 * 3)
    assert_max_deps(x.sum(axis=0, split_every=thresh), 2)
    assert_max_deps(x.sum(axis=1, split_every=thresh), 3)
    assert_max_deps(x.sum(split_every=20), 20, False)
    assert_max_deps(x.sum(axis=0, split_every=20), 4)
    assert_max_deps(x.sum(axis=1, split_every=20), 6)

    # 3D
    x = da.from_array(np.arange(11 * 22 * 29).reshape((11, 22, 29)), chunks=(3, 4, 5))
    thresh = {0: 2, 1: 3, 2: 4}
    assert_max_deps(x.sum(split_every=thresh), 2 * 3 * 4)
    assert_max_deps(x.sum(axis=0, split_every=thresh), 2)
    assert_max_deps(x.sum(axis=1, split_every=thresh), 3)
    assert_max_deps(x.sum(axis=2, split_every=thresh), 4)
    assert_max_deps(x.sum(axis=(0, 1), split_every=thresh), 2 * 3)
    assert_max_deps(x.sum(axis=(0, 2), split_every=thresh), 2 * 4)
    assert_max_deps(x.sum(axis=(1, 2), split_every=thresh), 3 * 4)
    assert_max_deps(x.sum(split_every=20), 20, False)
    assert_max_deps(x.sum(axis=0, split_every=20), 4)
    assert_max_deps(x.sum(axis=1, split_every=20), 6)
    assert_max_deps(x.sum(axis=2, split_every=20), 6)
    assert_max_deps(x.sum(axis=(0, 1), split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(0, 2), split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(1, 2), split_every=20), 20, False)
    assert_max_deps(x.sum(axis=(0, 1), split_every=40), 4 * 6)
    assert_max_deps(x.sum(axis=(0, 2), split_every=40), 4 * 6)
    assert_max_deps(x.sum(axis=(1, 2), split_every=40), 6 * 6)


def test_tree_reduce_set_options():
    x = da.from_array(np.arange(242).reshape((11, 22)), chunks=(3, 4))
    with set_options(split_every={0: 2, 1: 3}):
        assert_max_deps(x.sum(), 2 * 3)
        assert_max_deps(x.sum(axis=0), 2)


def test_reduction_names():
    x = da.ones(5, chunks=(2,))
    assert x.sum().name.startswith('sum')
    assert 'max' in x.max().name.split('-')[0]
    assert x.var().name.startswith('var')
    assert x.all().name.startswith('all')
    assert any(k[0].startswith('nansum') for k in da.nansum(x).dask)
    assert x.mean().name.startswith('mean')


@pytest.mark.skipif(np.__version__ < '1.12.0', reason='argmax out parameter')
@pytest.mark.parametrize('func', [np.sum,
                                  np.argmax])
def test_array_reduction_out(func):
    x = da.arange(10, chunks=(5,))
    y = da.ones((10, 10), chunks=(4, 4))
    func(y, axis=0, out=x)
    assert_eq(x, func(np.ones((10, 10)), axis=0))


@pytest.mark.parametrize("func", ["cumsum", "cumprod"])
@pytest.mark.parametrize("axis", [None, 0, 1, -1])
def test_array_cumreduction_axis(func, axis):
    np_func = getattr(np, func)
    da_func = getattr(da, func)

    s = (10, 11, 12)
    a = np.arange(np.prod(s)).reshape(s)
    d = da.from_array(a, chunks=(4, 5, 6))

    a_r = np_func(a, axis=axis)
    d_r = da_func(d, axis=axis)

    assert_eq(a_r, d_r)


@pytest.mark.parametrize('func', [np.cumsum, np.cumprod])
def test_array_cumreduction_out(func):
    x = da.ones((10, 10), chunks=(4, 4))
    func(x, axis=0, out=x)
    assert_eq(x, func(np.ones((10, 10)), axis=0))


@pytest.mark.parametrize('npfunc,daskfunc', [
    (np.sort, da.topk),
    (np.argsort, da.argtopk),
])
@pytest.mark.parametrize('split_every', [None, 2])
def test_topk_argtopk1(npfunc, daskfunc, split_every):
    # Test data
    k = 5
    a = da.random.random(1000, chunks=250)
    b = da.random.random((10, 20, 30), chunks=(4, 8, 8))
    c = da.from_array(np.array([(1, 'Hello'), (2, 'World')], dtype=[('foo', int), ('bar', '<U5')]),
                      chunks=1)

    # 1-dimensional arrays
    # top 5 elements, sorted descending
    assert_eq(npfunc(a)[-k:][::-1],
              daskfunc(a, k, split_every=split_every))
    # bottom 5 elements, sorted ascending
    assert_eq(npfunc(a)[:k],
              daskfunc(a, -k, split_every=split_every))

    # n-dimensional arrays
    # also testing when k > chunk
    # top 5 elements, sorted descending
    assert_eq(npfunc(b, axis=0)[-k:, :, :][::-1, :, :],
              daskfunc(b, k, axis=0, split_every=split_every))
    assert_eq(npfunc(b, axis=1)[:, -k:, :][:, ::-1, :],
              daskfunc(b, k, axis=1, split_every=split_every))
    assert_eq(npfunc(b, axis=-1)[:, :, -k:][:, :, ::-1],
              daskfunc(b, k, axis=-1, split_every=split_every))
    with pytest.raises(ValueError):
        daskfunc(b, k, axis=3, split_every=split_every)

    # bottom 5 elements, sorted ascending
    assert_eq(npfunc(b, axis=0)[:k, :, :],
              daskfunc(b, -k, axis=0, split_every=split_every))
    assert_eq(npfunc(b, axis=1)[:, :k, :],
              daskfunc(b, -k, axis=1, split_every=split_every))
    assert_eq(npfunc(b, axis=-1)[:, :, :k],
              daskfunc(b, -k, axis=-1, split_every=split_every))
    with pytest.raises(ValueError):
        daskfunc(b, -k, axis=3, split_every=split_every)

    # structured arrays
    assert_eq(npfunc(c, axis=0)[-1:][::-1],
              daskfunc(c, 1, split_every=split_every))
    assert_eq(npfunc(c, axis=0)[:1],
              daskfunc(c, -1, split_every=split_every))


def test_topk_argtopk2():
    a = da.random.random((10, 20, 30), chunks=(4, 8, 8))

    # Support for deprecated API for topk
    with pytest.warns(UserWarning):
        assert_eq(da.topk(a, 5), da.topk(5, a))

    # As Array methods
    assert_eq(a.topk(5, axis=1, split_every=2), da.topk(a, 5, axis=1, split_every=2))
    assert_eq(a.argtopk(5, axis=1, split_every=2), da.argtopk(a, 5, axis=1, split_every=2))
