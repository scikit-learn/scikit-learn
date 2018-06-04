from __future__ import division, print_function, absolute_import

import itertools
from numbers import Number
import textwrap

import pytest
from distutils.version import LooseVersion

np = pytest.importorskip('numpy')

import dask.array as da
from dask.utils import ignoring
from dask.array.utils import assert_eq, same_keys
from dask.array.einsumfuncs import einsum_can_optimize


def test_array():
    x = np.ones(5, dtype='i4')
    d = da.ones(5, chunks=3, dtype='i4')
    assert_eq(da.array(d, ndmin=3, dtype='i8'),
              np.array(x, ndmin=3, dtype='i8'))

    # regression #1847 this shall not raise an exception.
    x = da.ones((100,3), chunks=10)
    y = da.array(x)
    assert isinstance(y, da.Array)


@pytest.mark.parametrize("funcname", [
    "atleast_1d",
    "atleast_2d",
    "atleast_3d",
])
def test_atleast_nd_no_args(funcname):
    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    np_r_n = np_func()
    da_r_n = da_func()

    assert np_r_n == da_r_n


@pytest.mark.parametrize("funcname", [
    "atleast_1d",
    "atleast_2d",
    "atleast_3d",
])
@pytest.mark.parametrize("shape, chunks", [
    (tuple(), tuple()),
    ((4,), (2,)),
    ((4, 6), (2, 3)),
    ((4, 6, 8), (2, 3, 4)),
    ((4, 6, 8, 10), (2, 3, 4, 5)),
])
def test_atleast_nd_one_arg(funcname, shape, chunks):
    np_a = np.random.random(shape)
    da_a = da.from_array(np_a, chunks=chunks)

    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    np_r = np_func(np_a)
    da_r = da_func(da_a)

    assert_eq(np_r, da_r)


@pytest.mark.parametrize("funcname", [
    "atleast_1d",
    "atleast_2d",
    "atleast_3d",
])
@pytest.mark.parametrize("shape1, shape2", list(
    itertools.combinations_with_replacement(
        [
            tuple(),
            (4,),
            (4, 6),
            (4, 6, 8),
            (4, 6, 8, 10),
        ],
        2
    )
))
def test_atleast_nd_two_args(funcname, shape1, shape2):
    np_a_1 = np.random.random(shape1)
    da_a_1 = da.from_array(np_a_1, chunks=tuple(c // 2 for c in shape1))

    np_a_2 = np.random.random(shape2)
    da_a_2 = da.from_array(np_a_2, chunks=tuple(c // 2 for c in shape2))

    np_a_n = [np_a_1, np_a_2]
    da_a_n = [da_a_1, da_a_2]

    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    np_r_n = np_func(*np_a_n)
    da_r_n = da_func(*da_a_n)

    assert type(np_r_n) is type(da_r_n)

    assert len(np_r_n) == len(da_r_n)

    for np_r, da_r in zip(np_r_n, da_r_n):
        assert_eq(np_r, da_r)


def test_transpose():
    x = np.arange(240).reshape((4, 6, 10))
    d = da.from_array(x, (2, 3, 4))

    assert_eq(d.transpose((2, 0, 1)),
              x.transpose((2, 0, 1)))
    assert same_keys(d.transpose((2, 0, 1)), d.transpose((2, 0, 1)))

    assert_eq(d.transpose(2, 0, 1),
              x.transpose(2, 0, 1))
    assert same_keys(d.transpose(2, 0, 1), d.transpose(2, 0, 1))

    with pytest.raises(ValueError):
        d.transpose(1, 2)

    with pytest.raises(ValueError):
        d.transpose((1, 2))


def test_transpose_negative_axes():
    x = np.ones((2, 3, 4, 5))
    y = da.ones((2, 3, 4, 5), chunks=3)

    assert_eq(x.transpose([-1, -2, 0, 1]),
              y.transpose([-1, -2, 0, 1]))


def test_swapaxes():
    x = np.random.normal(0, 10, size=(10, 12, 7))
    d = da.from_array(x, chunks=(4, 5, 2))

    assert_eq(np.swapaxes(x, 0, 1), da.swapaxes(d, 0, 1))
    assert_eq(np.swapaxes(x, 2, 1), da.swapaxes(d, 2, 1))
    assert_eq(x.swapaxes(2, 1), d.swapaxes(2, 1))
    assert_eq(x.swapaxes(0, 0), d.swapaxes(0, 0))
    assert_eq(x.swapaxes(1, 2), d.swapaxes(1, 2))
    assert_eq(x.swapaxes(0, -1), d.swapaxes(0, -1))
    assert_eq(x.swapaxes(-1, 1), d.swapaxes(-1, 1))

    assert d.swapaxes(0, 1).name == d.swapaxes(0, 1).name
    assert d.swapaxes(0, 1).name != d.swapaxes(1, 0).name


@pytest.mark.parametrize("funcname, kwargs", [
    ("flipud", {}),
    ("fliplr", {}),
    ("flip", {"axis": 0}),
    ("flip", {"axis": 1}),
    ("flip", {"axis": 2}),
    ("flip", {"axis": -1}),
])
@pytest.mark.parametrize("shape", [
    tuple(),
    (4,),
    (4, 6),
    (4, 6, 8),
    (4, 6, 8, 10),
])
def test_flip(funcname, kwargs, shape):
    if (funcname == "flip" and
            LooseVersion(np.__version__) < LooseVersion("1.12.0")):
        pytest.skip(
            "NumPy %s doesn't support `flip`."
            " Need NumPy 1.12.0 or greater." % np.__version__
        )

    axis = kwargs.get("axis")
    if axis is None:
        if funcname == "flipud":
            axis = 0
        elif funcname == "fliplr":
            axis = 1

    np_a = np.random.random(shape)
    da_a = da.from_array(np_a, chunks=1)

    np_func = getattr(np, funcname)
    da_func = getattr(da, funcname)

    try:
        range(np_a.ndim)[axis]
    except IndexError:
        with pytest.raises(ValueError):
            da_func(da_a, **kwargs)
    else:
        np_r = np_func(np_a, **kwargs)
        da_r = da_func(da_a, **kwargs)

        assert_eq(np_r, da_r)


@pytest.mark.parametrize("x_shape, y_shape", [
    [(), ()],
    [(), (7,)],
    [(), (7, 11)],
    [(), (7, 11, 15)],
    [(), (7, 11, 15, 19)],
    [(7,), ()],
    [(7,), (7,)],
    [(11,), (11, 7)],
    [(15,), (7, 15, 11)],
    [(19,), (7, 11, 19, 15)],
    [(7, 11), ()],
    [(7, 11), (11,)],
    [(7, 11), (11, 7)],
    [(11, 15), (7, 15, 11)],
    [(15, 19), (7, 11, 19, 15)],
    [(7, 11, 15), ()],
    [(7, 11, 15), (15,)],
    [(7, 11, 15), (15, 7)],
    [(7, 11, 15), (7, 15, 11)],
    [(11, 15, 19), (7, 11, 19, 15)],
    [(7, 11, 15, 19), ()],
    [(7, 11, 15, 19), (19,)],
    [(7, 11, 15, 19), (19, 7)],
    [(7, 11, 15, 19), (11, 19, 13)],
    [(7, 11, 15, 19), (7, 11, 19, 15)],
])
def test_matmul(x_shape, y_shape):
    np.random.seed(3732)

    x = np.random.random(x_shape)[()]
    y = np.random.random(y_shape)[()]

    a = da.from_array(x, chunks=tuple((i // 2) for i in x.shape))
    b = da.from_array(y, chunks=tuple((i // 2) for i in y.shape))

    expected = None
    try:
        expected = np.matmul(x, y)
    except ValueError:
        pass

    for d1, d2 in itertools.product([a, x], [b, y]):
        if x.ndim == 0 or y.ndim == 0:
            with pytest.raises(ValueError):
                da.matmul(d1, d2)
        else:
            assert_eq(expected, da.matmul(d1, d2))


def test_tensordot():
    x = np.arange(400).reshape((20, 20))
    a = da.from_array(x, chunks=(5, 4))
    y = np.arange(200).reshape((20, 10))
    b = da.from_array(y, chunks=(4, 5))

    for axes in [1, (1, 0)]:
        assert_eq(da.tensordot(a, b, axes=axes), np.tensordot(x, y, axes=axes))
        assert_eq(da.tensordot(x, b, axes=axes), np.tensordot(x, y, axes=axes))
        assert_eq(da.tensordot(a, y, axes=axes), np.tensordot(x, y, axes=axes))

    assert same_keys(da.tensordot(a, b, axes=(1, 0)),
                     da.tensordot(a, b, axes=(1, 0)))
    with pytest.warns(None):  # Increasing number of chunks warning
        assert not same_keys(da.tensordot(a, b, axes=0),
                             da.tensordot(a, b, axes=1))


@pytest.mark.parametrize('axes', [
    0,
    1,
    (0, 1),
    (1, 0),
    ((1, 0), (2, 1)),
    ((1, 2), (2, 0)),
    ((2, 0), (1, 2))
])
def test_tensordot_2(axes):
    x = np.arange(4 * 4 * 4).reshape((4, 4, 4))
    y = da.from_array(x, chunks=2)

    assert_eq(da.tensordot(y, y, axes=axes),
              np.tensordot(x, x, axes=axes))


def test_dot_method():
    x = np.arange(400).reshape((20, 20))
    a = da.from_array(x, chunks=(5, 5))
    y = np.arange(200).reshape((20, 10))
    b = da.from_array(y, chunks=(5, 5))

    assert_eq(a.dot(b), x.dot(y))


@pytest.mark.parametrize('shape, chunks', [
    ((20,), (6,)),
    ((4, 5,), (2, 3)),
])
def test_vdot(shape, chunks):
    np.random.random(1337)

    x = 2 * np.random.random((2,) + shape) - 1
    x = x[0] + 1j * x[1]

    y = 2 * np.random.random((2,) + shape) - 1
    y = y[0] + 1j * y[1]

    a = da.from_array(x, chunks=chunks)
    b = da.from_array(y, chunks=chunks)

    assert_eq(np.vdot(x, y), da.vdot(a, b))
    assert_eq(np.vdot(y, x), da.vdot(b, a))
    assert_eq(da.vdot(a, b), da.vdot(b, a).conj())


@pytest.mark.parametrize('func1d_name, func1d', [
    ["ndim", lambda x: x.ndim],
    ["sum", lambda x: x.sum()],
    ["range", lambda x: [x.min(), x.max()]],
    ["range2", lambda x: [[x.min(), x.max()], [x.max(), x.min()]]],
])
@pytest.mark.parametrize('shape, axis', [
    [(10, 15, 20), 0],
    [(10, 15, 20), 1],
    [(10, 15, 20), 2],
    [(10, 15, 20), -1],
])
def test_apply_along_axis(func1d_name, func1d, shape, axis):
    a = np.random.randint(0, 10, shape)
    d = da.from_array(a, chunks=(len(shape) * (5,)))

    if (func1d_name == "range2" and
            LooseVersion(np.__version__) < LooseVersion("1.13.0")):
        with pytest.raises(ValueError):
            da.apply_along_axis(func1d, axis, d)
    else:
        assert_eq(
            da.apply_along_axis(func1d, axis, d),
            np.apply_along_axis(func1d, axis, a)
        )


@pytest.mark.parametrize('func_name, func', [
    ["sum0", lambda x, axis: x.sum(axis=axis)],
    ["sum1", lambda x, axis: x.sum(axis=axis, keepdims=True)],
    [
        "range", lambda x, axis:
            np.concatenate(
                [
                    x.min(axis=axis, keepdims=True),
                    x.max(axis=axis, keepdims=True)
                ],
                axis=axis
            )
    ],
])
@pytest.mark.parametrize('shape, axes', [
    [(10, 15, 20), tuple()],
    [(10, 15, 20), 0],
    [(10, 15, 20), (1,)],
    [(10, 15, 20), (-1, 1)],
    [(10, 15, 20), (2, 0, 1)],
])
def test_apply_over_axes(func_name, func, shape, axes):
    a = np.random.randint(0, 10, shape)
    d = da.from_array(a, chunks=(len(shape) * (5,)))

    assert_eq(
        da.apply_over_axes(func, d, axes),
        np.apply_over_axes(func, a, axes)
    )


@pytest.mark.parametrize('shape, axis', [
    [(10, 15, 20), None],
    [(10, 15, 20), 0],
    [(10, 15, 20), 1],
    [(10, 15, 20), 2],
    [(10, 15, 20), -1],
])
def test_ptp(shape, axis):
    a = np.random.randint(0, 10, shape)
    d = da.from_array(a, chunks=(len(shape) * (5,)))

    assert_eq(da.ptp(d, axis), np.ptp(a, axis))


@pytest.mark.parametrize('shape, axis', [
    [(10, 15, 20), 0],
    [(10, 15, 20), 1],
    [(10, 15, 20), 2],
    [(10, 15, 20), -1],
])
@pytest.mark.parametrize('n', [
    0,
    1,
    2,
])
def test_diff(shape, n, axis):
    x = np.random.randint(0, 10, shape)
    a = da.from_array(x, chunks=(len(shape) * (5,)))

    assert_eq(da.diff(a, n, axis), np.diff(x, n, axis))


@pytest.mark.parametrize('shape', [
    (10,),
    (10, 15),
])
@pytest.mark.parametrize('to_end, to_begin', [
    [None, None],
    [0, 0],
    [[1, 2], [3, 4]],
])
def test_ediff1d(shape, to_end, to_begin):
    x = np.random.randint(0, 10, shape)
    a = da.from_array(x, chunks=(len(shape) * (5,)))

    assert_eq(da.ediff1d(a, to_end, to_begin), np.ediff1d(x, to_end, to_begin))


@pytest.mark.parametrize('shape, varargs, axis', [
    [(10, 15, 20), (), None],
    [(10, 15, 20), (2,), None],
    [(10, 15, 20), (1.0, 1.5, 2.0), None],
    [(10, 15, 20), (), 0],
    [(10, 15, 20), (), 1],
    [(10, 15, 20), (), 2],
    [(10, 15, 20), (), -1],
    [(10, 15, 20), (), (0, 2)],
])
@pytest.mark.parametrize('edge_order', [
    1,
    2
])
def test_gradient(shape, varargs, axis, edge_order):
    a = np.random.randint(0, 10, shape)
    d_a = da.from_array(a, chunks=(len(shape) * (5,)))

    r = np.gradient(a, *varargs, axis=axis, edge_order=edge_order)
    r_a = da.gradient(d_a, *varargs, axis=axis, edge_order=edge_order)

    if isinstance(axis, Number):
        assert_eq(r, r_a)
    else:
        assert len(r) == len(r_a)

        for e_r, e_r_a in zip(r, r_a):
            assert_eq(e_r, e_r_a)


def test_bincount():
    x = np.array([2, 1, 5, 2, 1])
    d = da.from_array(x, chunks=2)
    e = da.bincount(d, minlength=6)
    assert_eq(e, np.bincount(x, minlength=6))
    assert same_keys(da.bincount(d, minlength=6), e)


def test_bincount_with_weights():
    x = np.array([2, 1, 5, 2, 1])
    d = da.from_array(x, chunks=2)
    weights = np.array([1, 2, 1, 0.5, 1])

    dweights = da.from_array(weights, chunks=2)
    e = da.bincount(d, weights=dweights, minlength=6)
    assert_eq(e, np.bincount(x, weights=dweights, minlength=6))
    assert same_keys(da.bincount(d, weights=dweights, minlength=6), e)


def test_bincount_raises_informative_error_on_missing_minlength_kwarg():
    x = np.array([2, 1, 5, 2, 1])
    d = da.from_array(x, chunks=2)
    try:
        da.bincount(d)
    except Exception as e:
        assert 'minlength' in str(e)
    else:
        assert False


def test_digitize():
    x = np.array([2, 4, 5, 6, 1])
    bins = np.array([1, 2, 3, 4, 5])
    for chunks in [2, 4]:
        for right in [False, True]:
            d = da.from_array(x, chunks=chunks)
            assert_eq(da.digitize(d, bins, right=right),
                      np.digitize(x, bins, right=right))

    x = np.random.random(size=(100, 100))
    bins = np.random.random(size=13)
    bins.sort()
    for chunks in [(10, 10), (10, 20), (13, 17), (87, 54)]:
        for right in [False, True]:
            d = da.from_array(x, chunks=chunks)
            assert_eq(da.digitize(d, bins, right=right),
                      np.digitize(x, bins, right=right))


def test_histogram():
    # Test for normal, flattened input
    n = 100
    v = da.random.random(n, chunks=10)
    bins = np.arange(0, 1.01, 0.01)
    (a1, b1) = da.histogram(v, bins=bins)
    (a2, b2) = np.histogram(v, bins=bins)

    # Check if the sum of the bins equals the number of samples
    assert a2.sum(axis=0) == n
    assert a1.sum(axis=0) == n
    assert_eq(a1, a2)
    assert same_keys(da.histogram(v, bins=bins)[0], a1)


def test_histogram_alternative_bins_range():
    v = da.random.random(100, chunks=10)
    (a1, b1) = da.histogram(v, bins=10, range=(0, 1))
    (a2, b2) = np.histogram(v, bins=10, range=(0, 1))
    assert_eq(a1, a2)
    assert_eq(b1, b2)


def test_histogram_return_type():
    v = da.random.random(100, chunks=10)
    bins = np.arange(0, 1.01, 0.01)
    # Check if return type is same as hist
    bins = np.arange(0, 11, 1, dtype='i4')
    assert_eq(da.histogram(v * 10, bins=bins)[0],
              np.histogram(v * 10, bins=bins)[0])


def test_histogram_extra_args_and_shapes():
    # Check for extra args and shapes
    bins = np.arange(0, 1.01, 0.01)
    v = da.random.random(100, chunks=10)
    data = [(v, bins, da.ones(100, chunks=v.chunks) * 5),
            (da.random.random((50, 50), chunks=10), bins, da.ones((50, 50), chunks=10) * 5)]

    for v, bins, w in data:
        # density
        assert_eq(da.histogram(v, bins=bins, normed=True)[0],
                  np.histogram(v, bins=bins, normed=True)[0])

        # normed
        assert_eq(da.histogram(v, bins=bins, density=True)[0],
                  np.histogram(v, bins=bins, density=True)[0])

        # weights
        assert_eq(da.histogram(v, bins=bins, weights=w)[0],
                  np.histogram(v, bins=bins, weights=w)[0])

        assert_eq(da.histogram(v, bins=bins, weights=w, density=True)[0],
                  da.histogram(v, bins=bins, weights=w, density=True)[0])


def test_cov():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(4, 4))

    assert_eq(da.cov(d), np.cov(x))
    assert_eq(da.cov(d, rowvar=0), np.cov(x, rowvar=0))
    with pytest.warns(None):  # warning dof <= 0 for slice
        assert_eq(da.cov(d, ddof=10), np.cov(x, ddof=10))
    assert_eq(da.cov(d, bias=1), np.cov(x, bias=1))
    assert_eq(da.cov(d, d), np.cov(x, x))

    y = np.arange(8)
    e = da.from_array(y, chunks=(4,))

    assert_eq(da.cov(d, e), np.cov(x, y))
    assert_eq(da.cov(e, d), np.cov(y, x))

    with pytest.raises(ValueError):
        da.cov(d, ddof=1.5)


def test_corrcoef():
    x = np.arange(56).reshape((7, 8))
    d = da.from_array(x, chunks=(4, 4))

    assert_eq(da.corrcoef(d), np.corrcoef(x))
    assert_eq(da.corrcoef(d, rowvar=0), np.corrcoef(x, rowvar=0))
    assert_eq(da.corrcoef(d, d), np.corrcoef(x, x))

    y = np.arange(8)
    e = da.from_array(y, chunks=(4,))

    assert_eq(da.corrcoef(d, e), np.corrcoef(x, y))
    assert_eq(da.corrcoef(e, d), np.corrcoef(y, x))


def test_round():
    x = np.random.random(10)
    d = da.from_array(x, chunks=4)

    for i in (0, 1, 4, 5):
        assert_eq(x.round(i), d.round(i))

    assert_eq(d.round(2), da.round(d, 2))


@pytest.mark.parametrize("return_index", [False, True])
@pytest.mark.parametrize("return_inverse", [False, True])
@pytest.mark.parametrize("return_counts", [False, True])
def test_unique_kwargs(return_index, return_inverse, return_counts):
    kwargs = dict(
        return_index=return_index,
        return_inverse=return_inverse,
        return_counts=return_counts
    )

    a = np.array([1, 2, 4, 4, 5, 2])
    d = da.from_array(a, chunks=(3,))

    r_a = np.unique(a, **kwargs)
    r_d = da.unique(d, **kwargs)

    if not any([return_index, return_inverse, return_counts]):
        assert isinstance(r_a, np.ndarray)
        assert isinstance(r_d, da.Array)

        r_a = (r_a,)
        r_d = (r_d,)

    assert len(r_a) == len(r_d)

    if return_inverse:
        i = 1 + int(return_index)
        assert (d.size,) == r_d[i].shape

    for e_r_a, e_r_d in zip(r_a, r_d):
        assert_eq(e_r_d, e_r_a)


@pytest.mark.parametrize("seed", [23, 796])
@pytest.mark.parametrize("low, high", [
    [0, 10]
])
@pytest.mark.parametrize("shape, chunks", [
    [(10,), (5,)],
    [(10,), (3,)],
    [(4, 5), (3, 2)],
    [(20, 20), (4, 5)],
])
def test_unique_rand(seed, low, high, shape, chunks):
    np.random.seed(seed)

    a = np.random.randint(low, high, size=shape)
    d = da.from_array(a, chunks=chunks)

    kwargs = dict(
        return_index=True,
        return_inverse=True,
        return_counts=True
    )

    r_a = np.unique(a, **kwargs)
    r_d = da.unique(d, **kwargs)

    assert len(r_a) == len(r_d)

    assert (d.size,) == r_d[2].shape

    for e_r_a, e_r_d in zip(r_a, r_d):
        assert_eq(e_r_d, e_r_a)


@pytest.mark.parametrize("seed", [23, 796])
@pytest.mark.parametrize("low, high", [
    [0, 10]
])
@pytest.mark.parametrize("elements_shape, elements_chunks", [
    [(10,), (5,)],
    [(10,), (3,)],
    [(4, 5), (3, 2)],
    [(20, 20), (4, 5)],
])
@pytest.mark.parametrize("test_shape, test_chunks", [
    [(10,), (5,)],
    [(10,), (3,)],
    [(4, 5), (3, 2)],
    [(20, 20), (4, 5)],
])
@pytest.mark.parametrize("invert", [True, False])
@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="np.isin is new in numpy 1.13")
def test_isin_rand(seed, low, high, elements_shape, elements_chunks,
                   test_shape, test_chunks, invert):
    rng = np.random.RandomState(seed)

    a1 = rng.randint(low, high, size=elements_shape)
    d1 = da.from_array(a1, chunks=elements_chunks)

    a2 = rng.randint(low, high, size=test_shape) - 5
    d2 = da.from_array(a2, chunks=test_chunks)

    r_a = np.isin(a1, a2, invert=invert)
    r_d = da.isin(d1, d2, invert=invert)
    assert_eq(r_a, r_d)


@pytest.mark.parametrize("assume_unique", [True, False])
@pytest.mark.skipif(LooseVersion(np.__version__) < '1.13.0',
                    reason="np.isin is new in numpy 1.13")
def test_isin_assume_unique(assume_unique):
    a1 = np.arange(10)
    d1 = da.from_array(a1, chunks=(5,))

    test_elements = np.arange(0, 10, 2)
    r_a = np.isin(a1, test_elements, assume_unique=assume_unique)
    r_d = da.isin(d1, test_elements, assume_unique=assume_unique)
    assert_eq(r_a, r_d)


def _maybe_len(l):
    try:
        return len(l)
    except TypeError:
        return 0


@pytest.mark.parametrize('chunks', [(4, 6), (2, 6)])
@pytest.mark.parametrize('shift', [3, 7, 9, (3, 9), (7, 2)])
@pytest.mark.parametrize('axis', [None, 0, 1, -1, (0, 1), (1, 0)])
def test_roll(chunks, shift, axis):
    x = np.random.randint(10, size=(4, 6))
    a = da.from_array(x, chunks=chunks)

    if _maybe_len(shift) != _maybe_len(axis):
        with pytest.raises(TypeError if axis is None else ValueError):
            da.roll(a, shift, axis)
    else:
        if (_maybe_len(shift) > 1 and
                LooseVersion(np.__version__) < LooseVersion("1.12.0")):
            pytest.skip(
                "NumPy %s doesn't support multiple axes with `roll`."
                " Need NumPy 1.12.0 or greater." % np.__version__
            )
        assert_eq(np.roll(x, shift, axis), da.roll(a, shift, axis))


def test_ravel():
    x = np.random.randint(10, size=(4, 6))

    # 2d
    for chunks in [(4, 6), (2, 6)]:
        a = da.from_array(x, chunks=chunks)
        assert_eq(x.ravel(), a.ravel())
        assert len(a.ravel().dask) == len(a.dask) + len(a.chunks[0])

    # 0d
    assert_eq(x[0, 0].ravel(), a[0, 0].ravel())

    # 1d
    a_flat = a.ravel()
    assert_eq(a_flat.ravel(), a_flat)

    # 3d
    x = np.random.randint(10, size=(2, 3, 4))
    for chunks in [4, (1, 3, 4)]:
        a = da.from_array(x, chunks=chunks)
        assert_eq(x.ravel(), a.ravel())

    assert_eq(x.flatten(), a.flatten())
    assert_eq(np.ravel(x), da.ravel(a))


@pytest.mark.parametrize('is_func', [True, False])
@pytest.mark.parametrize('axis', [None, 0, -1, (0, -1)])
def test_squeeze(is_func, axis):
    a = np.arange(10)[None, :, None, None]
    d = da.from_array(a, chunks=(1, 3, 1, 1))

    if is_func:
        a_s = np.squeeze(a, axis=axis)
        d_s = da.squeeze(d, axis=axis)
    else:
        a_s = a.squeeze(axis=axis)
        d_s = d.squeeze(axis=axis)

    assert_eq(d_s, a_s)
    assert same_keys(d_s, da.squeeze(d, axis=axis))

    if axis is None:
        axis = tuple(range(a.ndim))
    else:
        axis = axis if isinstance(axis, tuple) else (axis,)
        axis = tuple(i % a.ndim for i in axis)
    axis = tuple(
        i for i, c in enumerate(d.chunks) if i in axis and len(c) == 1
    )

    exp_d_s_chunks = tuple(
        c for i, c in enumerate(d.chunks) if i not in axis
    )
    assert d_s.chunks == exp_d_s_chunks


def test_vstack():
    x = np.arange(5)
    y = np.ones(5)
    a = da.arange(5, chunks=2)
    b = da.ones(5, chunks=2)

    assert_eq(np.vstack((x, y)), da.vstack((a, b)))
    assert_eq(np.vstack((x, y[None, :])), da.vstack((a, b[None, :])))


def test_hstack():
    x = np.arange(5)
    y = np.ones(5)
    a = da.arange(5, chunks=2)
    b = da.ones(5, chunks=2)

    assert_eq(np.hstack((x[None, :], y[None, :])),
              da.hstack((a[None, :], b[None, :])))
    assert_eq(np.hstack((x, y)), da.hstack((a, b)))


def test_dstack():
    x = np.arange(5)
    y = np.ones(5)
    a = da.arange(5, chunks=2)
    b = da.ones(5, chunks=2)

    assert_eq(np.dstack((x[None, None, :], y[None, None, :])),
              da.dstack((a[None, None, :], b[None, None, :])))
    assert_eq(np.dstack((x[None, :], y[None, :])),
              da.dstack((a[None, :], b[None, :])))
    assert_eq(np.dstack((x, y)), da.dstack((a, b)))


def test_take():
    x = np.arange(400).reshape((20, 20))
    a = da.from_array(x, chunks=(5, 5))

    assert_eq(np.take(x, 3, axis=0), da.take(a, 3, axis=0))
    assert_eq(np.take(x, [3, 4, 5], axis=-1), da.take(a, [3, 4, 5], axis=-1))

    with pytest.raises(ValueError):
        da.take(a, 3, axis=2)

    assert same_keys(da.take(a, [3, 4, 5], axis=-1),
                     da.take(a, [3, 4, 5], axis=-1))


def test_take_dask_from_numpy():
    x = np.arange(5).astype('f8')
    y = da.from_array(np.array([1, 2, 3, 3, 2 ,1]), chunks=3)

    z = da.take(x * 2, y)

    assert z.chunks == y.chunks
    assert_eq(z, np.array([2., 4., 6., 6., 4., 2.]))


def test_compress():
    x = np.arange(25).reshape((5, 5))
    a = da.from_array(x, chunks=(2, 2))

    c1 = np.array([True, False, True, False, True])
    c2 = np.array([True, False])
    c3 = [True, False]
    dc1 = da.from_array(c1, chunks=3)
    dc2 = da.from_array(c2, chunks=2)

    for c, dc in [(c1, c1), (c2, c2), (c3, c3),
                  (c1, dc1), (c2, dc2), (c3, dc2)]:
        for axis in [None, 0, 1]:
            res = da.compress(dc, a, axis=axis)
            assert_eq(np.compress(c, x, axis=axis), res)
            if isinstance(dc, da.Array):
                axis = axis or 0
                assert np.isnan(res.chunks[axis]).all()

    with pytest.raises(ValueError):
        da.compress([True, False], a, axis=100)

    with pytest.raises(ValueError):
        da.compress([[True], [False]], a, axis=100)


def test_extract():
    x = np.arange(25).reshape((5, 5))
    a = da.from_array(x, chunks=(2, 2))

    c1 = np.array([True, False, True, False, True])
    c2 = np.array([[True, False], [True, False]])
    c3 = np.array([True, False])
    dc1 = da.from_array(c1, chunks=3)
    dc2 = da.from_array(c2, chunks=(2, 1))
    dc3 = da.from_array(c3, chunks=2)

    for c, dc in [(c1, c1), (c2, c2), (c3, c3),
                  (c1, dc1), (c2, dc2), (c3, dc3)]:
        res = da.extract(dc, a)
        assert_eq(np.extract(c, x), res)
        if isinstance(dc, da.Array):
            assert np.isnan(res.chunks[0]).all()


def test_isnull():
    x = np.array([1, np.nan])
    a = da.from_array(x, chunks=(2,))
    with ignoring(ImportError):
        assert_eq(da.isnull(a), np.isnan(x))
        assert_eq(da.notnull(a), ~np.isnan(x))


def test_isclose():
    x = np.array([0, np.nan, 1, 1.5])
    y = np.array([1e-9, np.nan, 1, 2])
    a = da.from_array(x, chunks=(2,))
    b = da.from_array(y, chunks=(2,))
    assert_eq(da.isclose(a, b, equal_nan=True),
              np.isclose(x, y, equal_nan=True))


def test_allclose():
    n_a = np.array([0, np.nan, 1, 1.5])
    n_b = np.array([1e-9, np.nan, 1, 2])

    d_a = da.from_array(n_a, chunks=(2,))
    d_b = da.from_array(n_b, chunks=(2,))

    n_r = np.allclose(n_a, n_b, equal_nan=True)
    d_r = da.allclose(d_a, d_b, equal_nan=True)

    assert_eq(np.array(n_r)[()], d_r)


def test_choose():
    # test choose function
    x = np.random.randint(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))

    assert_eq(da.choose(d > 5, [0, d]), np.choose(x > 5, [0, x]))
    assert_eq(da.choose(d > 5, [-d, d]), np.choose(x > 5, [-x, x]))

    # test choose method
    index_dask = d > 5
    index_numpy = x > 5
    assert_eq(index_dask.choose([0, d]), index_numpy.choose([0, x]))
    assert_eq(index_dask.choose([-d, d]), index_numpy.choose([-x, x]))


def test_piecewise():
    np.random.seed(1337)

    x = np.random.randint(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))

    assert_eq(
        np.piecewise(x, [x < 5, x >= 5], [lambda e, v, k: e + 1, 5], 1, k=2),
        da.piecewise(d, [d < 5, d >= 5], [lambda e, v, k: e + 1, 5], 1, k=2)
    )


@pytest.mark.skipif(
    LooseVersion(np.__version__) < '1.12.0',
    reason=textwrap.dedent(
        """\
            NumPy piecewise mishandles the otherwise condition pre-1.12.0.

            xref: https://github.com/numpy/numpy/issues/5737
        """
    )
)
def test_piecewise_otherwise():
    np.random.seed(1337)

    x = np.random.randint(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))

    assert_eq(
        np.piecewise(
            x,
            [x > 5, x <= 2],
            [lambda e, v, k: e + 1, lambda e, v, k: v * e, lambda e, v, k: 0],
            1, k=2
        ),
        da.piecewise(
            d,
            [d > 5, d <= 2],
            [lambda e, v, k: e + 1, lambda e, v, k: v * e, lambda e, v, k: 0],
            1, k=2
        )
    )


def test_argwhere():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_nz = np.argwhere(x)
        d_nz = da.argwhere(d)

        assert_eq(d_nz, x_nz)


def test_argwhere_obj():
    x = np.random.randint(10, size=(15, 16)).astype(object)
    d = da.from_array(x, chunks=(4, 5))

    x_nz = np.argwhere(x)
    d_nz = da.argwhere(d)

    assert_eq(d_nz, x_nz)


def test_argwhere_str():
    x = np.array(list("Hello world"))
    d = da.from_array(x, chunks=(4,))

    x_nz = np.argwhere(x)
    d_nz = da.argwhere(d)

    assert_eq(d_nz, x_nz)


def test_where():
    x = np.random.randint(10, size=(15, 14))
    x[5, 5] = x[4, 4] = 0 # Ensure some false elements
    d = da.from_array(x, chunks=(4, 5))
    y = np.random.randint(10, size=15).astype(np.uint8)
    e = da.from_array(y, chunks=(4,))

    for c1, c2 in [(d > 5, x > 5),
                   (d, x),
                   (1, 1),
                   (0, 0),
                   (5, 5),
                   (True, True),
                   (np.True_, np.True_),
                   (False, False),
                   (np.False_, np.False_)]:
        for b1, b2 in [(0, 0), (-e[:, None], -y[:, None]), (e[:14], y[:14])]:
            w1 = da.where(c1, d, b1)
            w2 = np.where(c2, x, b2)
            assert_eq(w1, w2)


def test_where_scalar_dtype():
    x = np.int32(3)
    y1 = np.array([4, 5, 6], dtype=np.int16)
    c1 = np.array([1, 0, 1])
    y2 = da.from_array(y1, chunks=2)
    c2 = da.from_array(c1, chunks=2)
    w1 = np.where(c1, x, y1)
    w2 = da.where(c2, x, y2)
    assert_eq(w1, w2)
    # Test again for the bool optimization
    w3 = np.where(True, x, y1)
    w4 = da.where(True, x, y1)
    assert_eq(w3, w4)


def test_where_bool_optimization():
    x = np.random.randint(10, size=(15, 16))
    d = da.from_array(x, chunks=(4, 5))
    y = np.random.randint(10, size=(15, 16))
    e = da.from_array(y, chunks=(4, 5))

    for c in [True, False, np.True_, np.False_, 1, 0]:
        w1 = da.where(c, d, e)
        w2 = np.where(c, x, y)

        assert_eq(w1, w2)

        ex_w1 = d if c else e

        assert w1 is ex_w1


def test_where_nonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_w = np.where(x)
        d_w = da.where(d)

        assert isinstance(d_w, type(x_w))
        assert len(d_w) == len(x_w)

        for i in range(len(x_w)):
            assert_eq(d_w[i], x_w[i])


def test_where_incorrect_args():
    a = da.ones(5, chunks=3)

    for kwd in ["x", "y"]:
        kwargs = {kwd: a}
        try:
            da.where(a > 0, **kwargs)
        except ValueError as e:
            assert 'either both or neither of x and y should be given' in str(e)


def test_count_nonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_c = np.count_nonzero(x)
        d_c = da.count_nonzero(d)

        if d_c.shape == tuple():
            assert x_c == d_c.compute()
        else:
            assert_eq(x_c, d_c)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.12.0',
                    reason="NumPy's count_nonzero doesn't yet support axis")
@pytest.mark.parametrize('axis', [None, 0, (1,), (0, 1)])
def test_count_nonzero_axis(axis):
    for shape, chunks in [((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_c = np.count_nonzero(x, axis)
        d_c = da.count_nonzero(d, axis)

        if d_c.shape == tuple():
            assert x_c == d_c.compute()
        else:
            assert_eq(x_c, d_c)


def test_count_nonzero_obj():
    x = np.random.randint(10, size=(15, 16)).astype(object)
    d = da.from_array(x, chunks=(4, 5))

    x_c = np.count_nonzero(x)
    d_c = da.count_nonzero(d)

    if d_c.shape == tuple():
        assert x_c == d_c.compute()
    else:
        assert_eq(x_c, d_c)


@pytest.mark.skipif(LooseVersion(np.__version__) < '1.12.0',
                    reason="NumPy's count_nonzero doesn't yet support axis")
@pytest.mark.parametrize('axis', [None, 0, (1,), (0, 1)])
def test_count_nonzero_obj_axis(axis):
    x = np.random.randint(10, size=(15, 16)).astype(object)
    d = da.from_array(x, chunks=(4, 5))

    x_c = np.count_nonzero(x, axis)
    d_c = da.count_nonzero(d, axis)

    if d_c.shape == tuple():
        assert x_c == d_c.compute()
    else:
        #######################################################
        # Workaround oddness with Windows and object arrays.  #
        #                                                     #
        # xref: https://github.com/numpy/numpy/issues/9468    #
        #######################################################
        assert_eq(x_c.astype(np.intp), d_c)


def test_count_nonzero_str():
    x = np.array(list("Hello world"))
    d = da.from_array(x, chunks=(4,))

    x_c = np.count_nonzero(x)
    d_c = da.count_nonzero(d)

    assert x_c == d_c.compute()


def test_flatnonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_fnz = np.flatnonzero(x)
        d_fnz = da.flatnonzero(d)

        assert_eq(d_fnz, x_fnz)


def test_nonzero():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_nz = np.nonzero(x)
        d_nz = da.nonzero(d)

        assert isinstance(d_nz, type(x_nz))
        assert len(d_nz) == len(x_nz)

        for i in range(len(x_nz)):
            assert_eq(d_nz[i], x_nz[i])


def test_nonzero_method():
    for shape, chunks in [(0, ()), ((0, 0), (0, 0)), ((15, 16), (4, 5))]:
        x = np.random.randint(10, size=shape)
        d = da.from_array(x, chunks=chunks)

        x_nz = x.nonzero()
        d_nz = d.nonzero()

        assert isinstance(d_nz, type(x_nz))
        assert len(d_nz) == len(x_nz)

        for i in range(len(x_nz)):
            assert_eq(d_nz[i], x_nz[i])


def test_coarsen():
    x = np.random.randint(10, size=(24, 24))
    d = da.from_array(x, chunks=(4, 8))

    assert_eq(da.chunk.coarsen(np.sum, x, {0: 2, 1: 4}),
              da.coarsen(np.sum, d, {0: 2, 1: 4}))
    assert_eq(da.chunk.coarsen(np.sum, x, {0: 2, 1: 4}),
              da.coarsen(da.sum, d, {0: 2, 1: 4}))


def test_coarsen_with_excess():
    x = da.arange(10, chunks=5)
    assert_eq(da.coarsen(np.min, x, {0: 3}, trim_excess=True),
              np.array([0, 5]))
    assert_eq(da.coarsen(np.sum, x, {0: 3}, trim_excess=True),
              np.array([0 + 1 + 2, 5 + 6 + 7]))


def test_insert():
    x = np.random.randint(10, size=(10, 10))
    a = da.from_array(x, chunks=(5, 5))
    y = np.random.randint(10, size=(5, 10))
    b = da.from_array(y, chunks=(4, 4))

    assert_eq(np.insert(x, 0, -1, axis=0), da.insert(a, 0, -1, axis=0))
    assert_eq(np.insert(x, 3, -1, axis=-1), da.insert(a, 3, -1, axis=-1))
    assert_eq(np.insert(x, 5, -1, axis=1), da.insert(a, 5, -1, axis=1))
    assert_eq(np.insert(x, -1, -1, axis=-2), da.insert(a, -1, -1, axis=-2))
    assert_eq(np.insert(x, [2, 3, 3], -1, axis=1),
              da.insert(a, [2, 3, 3], -1, axis=1))
    assert_eq(np.insert(x, [2, 3, 8, 8, -2, -2], -1, axis=0),
              da.insert(a, [2, 3, 8, 8, -2, -2], -1, axis=0))
    assert_eq(np.insert(x, slice(1, 4), -1, axis=1),
              da.insert(a, slice(1, 4), -1, axis=1))
    assert_eq(np.insert(x, [2] * 3 + [5] * 2, y, axis=0),
              da.insert(a, [2] * 3 + [5] * 2, b, axis=0))
    assert_eq(np.insert(x, 0, y[0], axis=1),
              da.insert(a, 0, b[0], axis=1))

    assert same_keys(da.insert(a, [2, 3, 8, 8, -2, -2], -1, axis=0),
                     da.insert(a, [2, 3, 8, 8, -2, -2], -1, axis=0))

    with pytest.raises(NotImplementedError):
        da.insert(a, [4, 2], -1, axis=0)

    with pytest.raises(IndexError):
        da.insert(a, [3], -1, axis=2)

    with pytest.raises(IndexError):
        da.insert(a, [3], -1, axis=-3)


def test_multi_insert():
    z = np.random.randint(10, size=(1, 2))
    c = da.from_array(z, chunks=(1, 2))
    assert_eq(np.insert(np.insert(z, [0, 1], -1, axis=0), [1], -1, axis=1),
              da.insert(da.insert(c, [0, 1], -1, axis=0), [1], -1, axis=1))


def test_result_type():
    a = da.from_array(np.ones(5, np.float32), chunks=(3,))
    b = da.from_array(np.ones(5, np.int16), chunks=(3,))
    c = da.from_array(np.ones(5, np.int64), chunks=(3,))
    x = np.ones(5, np.float32)
    assert da.result_type(b, c) == np.int64
    assert da.result_type(a, b, c) == np.float64
    assert da.result_type(b, np.float32) == np.float32
    assert da.result_type(b, np.dtype(np.float32)) == np.float32
    assert da.result_type(b, x) == np.float32
    # Effect of scalars depends on their value
    assert da.result_type(1, b) == np.int16
    assert da.result_type(1.0, a) == np.float32
    assert da.result_type(np.int64(1), b) == np.int16
    assert da.result_type(np.ones((), np.int64), b) == np.int16  # 0d array
    assert da.result_type(1e200, a) == np.float64   # 1e200 is too big for float32
    # dask 0d-arrays are NOT treated like scalars
    c = da.from_array(np.ones((), np.float64), chunks=())
    assert da.result_type(a, c) == np.float64


def _numpy_and_dask_inputs(input_sigs):
    # einsum label dimensions
    _dimensions = {'a': 5, 'b': 6, 'c': 7,
                   'd': 5, 'e': 6, 'f': 10,
                   'g': 1, 'h': 2, '*': 11}

    # dimension chunks sizes
    _chunks = {'a': (2, 3), 'b': (2, 3, 1), 'c': (2, 3, 2),
               'd': (4, 1), 'e': (2, 4),    'f': (1, 2, 3, 4),
               'g': 1,      'h': (1, 1),    '*': 11}

    def _shape_from_string(s):
        return tuple(_dimensions[c] for c in s)

    def _chunks_from_string(s):
        return tuple(_chunks[c] for c in s)

    shapes = [_shape_from_string(s) for s in input_sigs]
    chunks = [_chunks_from_string(s) for s in input_sigs]

    np_inputs = [np.random.random(s) for s in shapes]
    da_inputs = [da.from_array(i, chunks=c) for i, c in zip(np_inputs, chunks)]

    return np_inputs, da_inputs


@pytest.mark.parametrize('einsum_signature', [
    'abc,bad->abcd',
    'abcdef,bcdfg->abcdeg',
    'ea,fb,abcd,gc,hd->efgh',
    'ab,b',
    'aa',
    'a,a->',
    'a,a->a',
    'a,a',
    'a,b',
    'a,b,c',
    'a',
    'ba,b',
    'ba,b->',
    'defab,fedbc->defac',
    'ab...,bc...->ac...',
    'a...a',
    'abc...->cba...',
    '...ab->...a',
    'a...a->a...',
    # Following 2 from # https://stackoverflow.com/a/19203475/1611416
    '...abc,...abcd->...d',
    'ab...,b->ab...',
    # https://github.com/dask/dask/pull/3412#discussion_r182413444
    'aa->a',
    'ab,ab,c->c',
    'aab,bc->ac',
    'aab,bcc->ac',
    'fdf,cdd,ccd,afe->ae',
    'fff,fae,bef,def->abd',
])
def test_einsum(einsum_signature):
    input_sigs = (einsum_signature.split('->')[0]
                                  .replace("...", "*")
                                  .split(','))

    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    assert_eq(np.einsum(einsum_signature, *np_inputs),
              da.einsum(einsum_signature, *da_inputs))


@pytest.mark.skipif(not einsum_can_optimize,
                    reason="np.einsum(optimize) unavailable")
@pytest.mark.parametrize('optimize_opts', [
    (True, False),
    ('greedy', False),
    ('optimal', False)
])
def test_einsum_optimize(optimize_opts):
    sig = 'ea,fb,abcd,gc,hd->efgh'
    input_sigs = sig.split('->')[0].split(',')
    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    opt1, opt2 = optimize_opts

    assert_eq(np.einsum(sig, *np_inputs, optimize=opt1),
              da.einsum(sig, *np_inputs, optimize=opt2))

    assert_eq(np.einsum(sig, *np_inputs, optimize=opt2),
              da.einsum(sig, *np_inputs, optimize=opt1))


@pytest.mark.parametrize('order', ['C', 'F', 'A', 'K'])
def test_einsum_order(order):
    sig = 'ea,fb,abcd,gc,hd->efgh'
    input_sigs = sig.split('->')[0].split(',')
    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    assert_eq(np.einsum(sig, *np_inputs, order=order),
              da.einsum(sig, *np_inputs, order=order))


@pytest.mark.parametrize('casting', [
    'no', 'equiv', 'safe', 'same_kind', 'unsafe'])
def test_einsum_casting(casting):
    sig = 'ea,fb,abcd,gc,hd->efgh'
    input_sigs = sig.split('->')[0].split(',')
    np_inputs, da_inputs = _numpy_and_dask_inputs(input_sigs)

    assert_eq(np.einsum(sig, *np_inputs, casting=casting),
              da.einsum(sig, *np_inputs, casting=casting))


@pytest.mark.parametrize('split_every', [None, 2])
def test_einsum_split_every(split_every):
    np_inputs, da_inputs = _numpy_and_dask_inputs('a')
    assert_eq(np.einsum('a', *np_inputs),
              da.einsum('a', *da_inputs, split_every=split_every))


def test_einsum_invalid_args():
    _, da_inputs = _numpy_and_dask_inputs('a')
    with pytest.raises(TypeError):
        da.einsum('a', *da_inputs, foo=1, bar=2)


def test_einsum_broadcasting_contraction():
    a = np.random.rand(1, 5, 4)
    b = np.random.rand(4, 6)
    c = np.random.rand(5, 6)
    d = np.random.rand(10)

    d_a = da.from_array(a, chunks=(1, (2, 3), (2, 2)))
    d_b = da.from_array(b, chunks=((2, 2), (4, 2)))
    d_c = da.from_array(c, chunks=((2, 3), (4, 2)))
    d_d = da.from_array(d, chunks=((7, 3)))

    np_res = np.einsum('ijk,kl,jl', a, b, c)
    da_res = da.einsum('ijk,kl,jl', d_a, d_b, d_c)
    assert_eq(np_res, da_res)

    mul_res = da_res * d

    np_res = np.einsum('ijk,kl,jl,i->i', a, b, c, d)
    da_res = da.einsum('ijk,kl,jl,i->i', d_a, d_b, d_c, d_d)
    assert_eq(np_res, da_res)
    assert_eq(np_res, mul_res)


def test_einsum_broadcasting_contraction2():
    a = np.random.rand(1, 1, 5, 4)
    b = np.random.rand(4, 6)
    c = np.random.rand(5, 6)
    d = np.random.rand(7, 7)

    d_a = da.from_array(a, chunks=(1, 1, (2, 3), (2, 2)))
    d_b = da.from_array(b, chunks=((2, 2), (4, 2)))
    d_c = da.from_array(c, chunks=((2, 3), (4, 2)))
    d_d = da.from_array(d, chunks=((7, 3)))

    np_res = np.einsum('abjk,kl,jl', a, b, c)
    da_res = da.einsum('abjk,kl,jl', d_a, d_b, d_c)
    assert_eq(np_res, da_res)

    mul_res = da_res * d

    np_res = np.einsum('abjk,kl,jl,ab->ab', a, b, c, d)
    da_res = da.einsum('abjk,kl,jl,ab->ab', d_a, d_b, d_c, d_d)
    assert_eq(np_res, da_res)
    assert_eq(np_res, mul_res)


def test_einsum_broadcasting_contraction3():
    a = np.random.rand(1, 5, 4)
    b = np.random.rand(4, 1, 6)
    c = np.random.rand(5, 6)
    d = np.random.rand(7, 7)

    d_a = da.from_array(a, chunks=(1, (2, 3), (2, 2)))
    d_b = da.from_array(b, chunks=((2, 2), 1,  (4, 2)))
    d_c = da.from_array(c, chunks=((2, 3), (4, 2)))
    d_d = da.from_array(d, chunks=((7, 3)))

    np_res = np.einsum('ajk,kbl,jl,ab->ab', a, b, c, d)
    da_res = da.einsum('ajk,kbl,jl,ab->ab', d_a, d_b, d_c, d_d)
    assert_eq(np_res, da_res)
