from __future__ import absolute_import, division, print_function

import pickle
from functools import partial
from operator import add

import pytest
np = pytest.importorskip('numpy')

import dask.array as da
from dask.array.ufunc import da_frompyfunc
from dask.array.utils import assert_eq
from dask.base import tokenize


def test_ufunc_meta():
    assert da.log.__name__ == 'log'
    assert da.log.__doc__.replace('  # doctest: +SKIP', '') == np.log.__doc__

    assert da.modf.__name__ == 'modf'
    assert da.modf.__doc__.replace('  # doctest: +SKIP', '') == np.modf.__doc__

    assert da.frexp.__name__ == 'frexp'
    assert da.frexp.__doc__.replace('  # doctest: +SKIP', '') == np.frexp.__doc__


def test_ufunc():
    for attr in ['nin', 'nargs', 'nout', 'ntypes', 'identity',
                 'signature', 'types']:
        assert getattr(da.log, attr) == getattr(np.log, attr)

    with pytest.raises(AttributeError):
        da.log.not_an_attribute

    assert repr(da.log) == repr(np.log)
    assert 'nin' in dir(da.log)
    assert 'outer' in dir(da.log)


binary_ufuncs = ['add', 'arctan2', 'copysign', 'divide', 'equal',
                 'bitwise_and', 'bitwise_or', 'bitwise_xor',
                 'floor_divide', 'fmax', 'fmin', 'fmod', 'greater',
                 'greater_equal', 'hypot', 'ldexp', 'less', 'less_equal',
                 'logaddexp', 'logaddexp2', 'logical_and', 'logical_or',
                 'logical_xor', 'maximum', 'minimum', 'mod', 'multiply',
                 'nextafter', 'not_equal', 'power', 'remainder', 'subtract',
                 'true_divide']

try:
    da.float_power
    binary_ufuncs += ['float_power']
except AttributeError:
    # Absent for NumPy versions prior to 1.12.
    pass

unary_ufuncs = ['absolute', 'arccos', 'arccosh', 'arcsin', 'arcsinh', 'arctan',
                'arctanh', 'bitwise_not', 'cbrt', 'ceil', 'conj', 'cos',
                'cosh', 'deg2rad', 'degrees', 'exp', 'exp2', 'expm1', 'fabs',
                'fix', 'floor', 'i0', 'isfinite', 'isinf', 'isnan', 'log',
                'log10', 'log1p', 'log2', 'logical_not', 'nan_to_num',
                'negative', 'rad2deg', 'radians', 'reciprocal', 'rint', 'sign',
                'signbit', 'sin', 'sinc', 'sinh', 'spacing', 'sqrt', 'square',
                'tan', 'tanh', 'trunc']


@pytest.mark.parametrize('ufunc', unary_ufuncs)
def test_unary_ufunc(ufunc):
    if ufunc == 'fix' and np.__version__ >= '1.13.0':
        pytest.skip('fix calls floor in a way that we do not yet support')
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    arr = np.random.randint(1, 100, size=(20, 20))
    darr = da.from_array(arr, 3)

    with pytest.warns(None):  # some invalid values (arccos, arcsin, etc.)
        # applying Dask ufunc doesn't trigger computation
        assert isinstance(dafunc(darr), da.Array)
        assert_eq(dafunc(darr), npfunc(arr), equal_nan=True)

    with pytest.warns(None):  # some invalid values (arccos, arcsin, etc.)
        # applying NumPy ufunc is lazy
        if isinstance(npfunc, np.ufunc) and np.__version__ >= '1.13.0':
            assert isinstance(npfunc(darr), da.Array)
        else:
            assert isinstance(npfunc(darr), np.ndarray)
        assert_eq(npfunc(darr), npfunc(arr), equal_nan=True)

    with pytest.warns(None):  # some invalid values (arccos, arcsin, etc.)
        # applying Dask ufunc to normal ndarray triggers computation
        assert isinstance(dafunc(arr), np.ndarray)
        assert_eq(dafunc(arr), npfunc(arr), equal_nan=True)


@pytest.mark.parametrize('ufunc', binary_ufuncs)
def test_binary_ufunc(ufunc):
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    arr1 = np.random.randint(1, 100, size=(20, 20))
    darr1 = da.from_array(arr1, 3)

    arr2 = np.random.randint(1, 100, size=(20, 20))
    darr2 = da.from_array(arr2, 3)

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(darr1, darr2), da.Array)
    assert_eq(dafunc(darr1, darr2), npfunc(arr1, arr2))

    # applying NumPy ufunc triggers computation or is lazy (np >= 1.13.0)
    if np.__version__ >= '1.13.0':
        assert isinstance(npfunc(darr1, darr2), da.Array)
    else:
        assert isinstance(npfunc(darr1, darr2), np.ndarray)
    assert_eq(npfunc(darr1, darr2), npfunc(arr1, arr2))

    # applying Dask ufunc to normal ndarray triggers computation
    assert isinstance(dafunc(arr1, arr2), np.ndarray)
    assert_eq(dafunc(arr1, arr2), npfunc(arr1, arr2))

    # with scalar
    assert isinstance(dafunc(darr1, 10), da.Array)
    assert_eq(dafunc(darr1, 10), npfunc(arr1, 10))

    with pytest.warns(None):  # overflow in ldexp
        assert isinstance(dafunc(10, darr1), da.Array)
        assert_eq(dafunc(10, darr1), npfunc(10, arr1))

    assert isinstance(dafunc(arr1, 10), np.ndarray)
    assert_eq(dafunc(arr1, 10), npfunc(arr1, 10))

    with pytest.warns(None):  # overflow in ldexp
        assert isinstance(dafunc(10, arr1), np.ndarray)
        assert_eq(dafunc(10, arr1), npfunc(10, arr1))


def test_ufunc_outer():
    arr1 = np.random.randint(1, 100, size=20)
    darr1 = da.from_array(arr1, 3)

    arr2 = np.random.randint(1, 100, size=(10, 3))
    darr2 = da.from_array(arr2, 3)

    # Check output types
    assert isinstance(da.add.outer(darr1, darr2), da.Array)
    assert isinstance(da.add.outer(arr1, darr2), da.Array)
    assert isinstance(da.add.outer(darr1, arr2), da.Array)
    assert isinstance(da.add.outer(arr1, arr2), np.ndarray)

    # Check mix of dimensions, dtypes, and numpy/dask/object
    cases = [((darr1, darr2), (arr1, arr2)),
             ((darr2, darr1), (arr2, arr1)),
             ((darr2, darr1.astype('f8')), (arr2, arr1.astype('f8'))),
             ((darr1, arr2), (arr1, arr2)),
             ((darr1, 1), (arr1, 1)),
             ((1, darr2), (1, arr2)),
             ((1.5, darr2), (1.5, arr2)),
             (([1, 2, 3], darr2), ([1, 2, 3], arr2)),
             ((darr1.sum(), darr2), (arr1.sum(), arr2)),
             ((np.array(1), darr2), (np.array(1), arr2))]

    for (dA, dB), (A, B) in cases:
        assert_eq(da.add.outer(dA, dB), np.add.outer(A, B))

    # Check dtype kwarg works
    assert_eq(da.add.outer(darr1, darr2, dtype='f8'),
              np.add.outer(arr1, arr2, dtype='f8'))

    with pytest.raises(ValueError):
        da.add.outer(darr1, darr2, out=arr1)

    with pytest.raises(ValueError):
        da.sin.outer(darr1, darr2)


@pytest.mark.parametrize('ufunc', ['isreal', 'iscomplex', 'real', 'imag'])
def test_complex(ufunc):

    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    real = np.random.randint(1, 100, size=(20, 20))
    imag = np.random.randint(1, 100, size=(20, 20)) * 1j
    comp = real + imag

    dareal = da.from_array(real, 3)
    daimag = da.from_array(imag, 3)
    dacomp = da.from_array(comp, 3)

    assert_eq(dacomp.real, comp.real)
    assert_eq(dacomp.imag, comp.imag)
    assert_eq(dacomp.conj(), comp.conj())

    for darr, arr in [(dacomp, comp), (dareal, real), (daimag, imag)]:
        # applying Dask ufunc doesn't trigger computation
        assert isinstance(dafunc(darr), da.Array)
        assert_eq(dafunc(darr), npfunc(arr))

        # applying NumPy ufunc triggers computation
        if np.__version__ < '1.13.0':
            assert isinstance(npfunc(darr), np.ndarray)
        assert_eq(npfunc(darr), npfunc(arr))

        # applying Dask ufunc to normal ndarray triggers computation
        assert isinstance(dafunc(arr), np.ndarray)
        assert_eq(dafunc(arr), npfunc(arr))


@pytest.mark.parametrize('ufunc', ['frexp', 'modf'])
def test_ufunc_2results(ufunc):

    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    arr = np.random.randint(1, 100, size=(20, 20))
    darr = da.from_array(arr, 3)

    # applying Dask ufunc doesn't trigger computation
    res1, res2 = dafunc(darr)
    assert isinstance(res1, da.Array)
    assert isinstance(res2, da.Array)
    exp1, exp2 = npfunc(arr)
    assert_eq(res1, exp1)
    assert_eq(res2, exp2)

    # applying NumPy ufunc is now lazy
    res1, res2 = npfunc(darr)
    if np.__version__ >= '1.13.0':
        assert isinstance(res1, da.Array)
        assert isinstance(res2, da.Array)
    else:
        assert isinstance(res1, np.ndarray)
        assert isinstance(res2, np.ndarray)
    exp1, exp2 = npfunc(arr)
    assert_eq(res1, exp1)
    assert_eq(res2, exp2)

    # applying Dask ufunc to normal ndarray triggers computation
    res1, res2 = dafunc(arr)
    assert isinstance(res1, da.Array)
    assert isinstance(res2, da.Array)
    exp1, exp2 = npfunc(arr)
    assert_eq(res1, exp1)
    assert_eq(res2, exp2)


def test_clip():
    x = np.random.normal(0, 10, size=(10, 10))
    d = da.from_array(x, chunks=(3, 4))

    assert_eq(x.clip(5), d.clip(5))
    assert_eq(x.clip(1, 5), d.clip(1, 5))
    assert_eq(x.clip(min=5), d.clip(min=5))
    assert_eq(x.clip(max=5), d.clip(max=5))
    assert_eq(x.clip(max=1, min=5), d.clip(max=1, min=5))
    assert_eq(x.clip(min=1, max=5), d.clip(min=1, max=5))


def test_angle():
    real = np.random.randint(1, 100, size=(20, 20))
    imag = np.random.randint(1, 100, size=(20, 20)) * 1j
    comp = real + imag
    dacomp = da.from_array(comp, 3)

    assert_eq(da.angle(dacomp), np.angle(comp))
    assert_eq(da.angle(dacomp, deg=True), np.angle(comp, deg=True))
    assert isinstance(da.angle(comp), np.ndarray)
    assert_eq(da.angle(comp), np.angle(comp))


def test_frompyfunc():
    myadd = da.frompyfunc(add, 2, 1)
    np_myadd = np.frompyfunc(add, 2, 1)

    x = np.random.normal(0, 10, size=(10, 10))
    dx = da.from_array(x, chunks=(3, 4))
    y = np.random.normal(0, 10, size=10)
    dy = da.from_array(y, chunks=2)

    assert_eq(myadd(dx, dy), np_myadd(x, y))
    assert_eq(myadd.outer(dx, dy), np_myadd.outer(x, y))

    with pytest.raises(NotImplementedError):
        da.frompyfunc(lambda x, y: (x + y, x - y), 2, 2)


def test_frompyfunc_wrapper():
    f = da_frompyfunc(add, 2, 1)
    np_f = np.frompyfunc(add, 2, 1)
    x = np.array([1, 2, 3])

    # Callable
    np.testing.assert_equal(f(x, 1), np_f(x, 1))

    # picklable
    f2 = pickle.loads(pickle.dumps(f))
    np.testing.assert_equal(f2(x, 1), np_f(x, 1))

    # Attributes
    assert f.ntypes == np_f.ntypes
    with pytest.raises(AttributeError):
        f.not_an_attribute

    # Tab completion
    assert 'ntypes' in dir(f)

    # Methods
    np.testing.assert_equal(f.outer(x, x), np_f.outer(x, x))

    # funcname
    assert f.__name__ == 'frompyfunc-add'

    # repr
    assert repr(f) == "da.frompyfunc<add, 2, 1>"

    # tokenize
    assert (tokenize(da_frompyfunc(add, 2, 1)) ==
            tokenize(da_frompyfunc(add, 2, 1)))


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_array_ufunc():
    x = np.arange(24).reshape((4, 6))
    d = da.from_array(x, chunks=(2, 3))

    for func in [np.sin, np.isreal, np.sum, np.negative, partial(np.prod, axis=0)]:
        assert isinstance(func(d), da.Array)
        assert_eq(func(d), func(x))


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_array_ufunc_binop():
    x = np.arange(25).reshape((5, 5))
    d = da.from_array(x, chunks=(2, 2))

    for func in [np.add, np.multiply]:
        assert isinstance(func(d, d), da.Array)
        assert_eq(func(d, d), func(x, x))

        assert isinstance(func.outer(d, d), da.Array)
        assert_eq(func.outer(d, d), func.outer(x, x))


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_array_ufunc_out():
    x = da.arange(10, chunks=(5,))
    np.sin(x, out=x)
    np.add(x, 10, out=x)
    assert_eq(x, np.sin(np.arange(10)) + 10)


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_unsupported_ufunc_methods():
    x = da.arange(10, chunks=(5,))
    with pytest.raises(TypeError):
        assert np.add.reduce(x)


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_out_numpy():
    x = da.arange(10, chunks=(5,))
    empty = np.empty(10, dtype=x.dtype)
    with pytest.raises((TypeError, NotImplementedError)) as info:
        np.add(x, 1, out=empty)

    assert 'ndarray' in str(info.value)
    assert 'Array' in str(info.value)


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_out_shape_mismatch():
    x = da.arange(10, chunks=(5,))
    y = da.arange(15, chunks=(5,))
    with pytest.raises(ValueError):
        assert np.log(x, out=y)
