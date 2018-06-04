from __future__ import absolute_import, division, print_function

import pytest

pd = pytest.importorskip('pandas')
import pandas.util.testing as tm

import numpy as np

import dask.array as da
import dask.dataframe as dd
from dask.dataframe.utils import assert_eq


_BASE_UFUNCS = ['conj', 'exp', 'log', 'log2', 'log10', 'log1p',
                'expm1', 'sqrt', 'square', 'sin', 'cos', 'tan',
                'arcsin','arccos', 'arctan', 'sinh', 'cosh', 'tanh',
                'arcsinh', 'arccosh', 'arctanh', 'deg2rad', 'rad2deg',
                'isfinite', 'isinf', 'isnan', 'signbit',
                'degrees', 'radians', 'rint', 'fabs', 'sign', 'absolute',
                'floor', 'ceil', 'trunc', 'logical_not', 'cbrt', 'exp2',
                'negative', 'reciprocal', 'spacing']


@pytest.mark.parametrize('pandas_input', [
                         pd.Series(np.random.randint(1, 100, size=20)),
                         pd.Series(np.abs(np.random.randn(100))),
                         pd.DataFrame({'A': np.random.randint(1, 100, size=20),
                                       'B': np.random.randint(1, 100, size=20),
                                       'C': np.abs(np.random.randn(20))}),
                         pd.Series(np.random.randint(1, 100, size=20),
                                   index=list('abcdefghijklmnopqrst')),
                         pd.Series(np.abs(np.random.randn(20)),
                                   index=list('abcdefghijklmnopqrst')),
                         pd.DataFrame({'A': np.random.randint(1, 100, size=20),
                                       'B': np.random.randint(1, 100, size=20),
                                       'C': np.abs(np.random.randn(20))},
                                      index=list('abcdefghijklmnopqrst'))])
@pytest.mark.parametrize('ufunc', _BASE_UFUNCS)
def test_ufunc(pandas_input, ufunc):

    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    dask_input = dd.from_pandas(pandas_input, 3)
    pandas_type = pandas_input.__class__
    dask_type = dask_input.__class__

    # applying Dask ufunc doesn't trigger computation
    with pytest.warns(None):
        # Some cause warnings (arcsine)
        assert isinstance(dafunc(dask_input), dask_type)
        assert_eq(dafunc(dask_input), npfunc(pandas_input))

        # applying NumPy ufunc is lazy
        if isinstance(npfunc, np.ufunc) and np.__version__ >= '1.13.0':
            assert isinstance(npfunc(dask_input), dask_type)
        else:
            assert isinstance(npfunc(dask_input), pandas_type)
        assert_eq(npfunc(dask_input), npfunc(pandas_input))

        # applying Dask ufunc to normal Series triggers computation
        assert isinstance(dafunc(pandas_input), pandas_type)
        assert_eq(dafunc(dask_input), npfunc(pandas_input))

    # Index
    if pandas_input.index.dtype in [object, str]:
        return
    if ufunc in ('logical_not', 'signbit', 'isnan', 'isinf', 'isfinite'):
        return

    with pytest.warns(None):
        assert isinstance(dafunc(dask_input.index), dd.Index)
        assert_eq(dafunc(dask_input.index), npfunc(pandas_input.index))

        # applying NumPy ufunc is lazy
        if isinstance(npfunc, np.ufunc) and np.__version__ >= '1.13.0':
            assert isinstance(npfunc(dask_input.index), dd.Index)
        else:
            assert isinstance(npfunc(dask_input.index), pd.Index)

        assert_eq(npfunc(dask_input.index), npfunc(dask_input.index))

    # applying Dask ufunc to normal Series triggers computation
    with pytest.warns(None):
        # some (da.log) cause warnings
        assert isinstance(dafunc(pandas_input.index), pd.Index)
        assert_eq(dafunc(pandas_input), npfunc(pandas_input))


@pytest.mark.parametrize('ufunc', ['isreal', 'iscomplex', 'real', 'imag',
                                   'angle', 'fix', 'i0', 'sinc', 'nan_to_num'])
def test_ufunc_array_wrap(ufunc):
    """
    some np.ufuncs doesn't call __array_wrap__
    (or __array_ufunc__ starting from numpy v.1.13.0), it should work as below

    - da.ufunc(dd.Series) => dd.Series
    - da.ufunc(pd.Series) => np.ndarray
    - np.ufunc(dd.Series) => np.ndarray
    - np.ufunc(pd.Series) => np.ndarray
    """
    if ufunc == 'fix' and np.__version__ >= '1.13.0':
        pytest.skip('fix calls floor in a way that we do not yet support')

    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    s = pd.Series(np.random.randint(1, 100, size=20),
                  index=list('abcdefghijklmnopqrst'))
    ds = dd.from_pandas(s, 3)

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(ds), dd.Series)
    assert_eq(dafunc(ds), pd.Series(npfunc(s), index=s.index))

    assert isinstance(npfunc(ds), np.ndarray)
    tm.assert_numpy_array_equal(npfunc(ds), npfunc(s))

    assert isinstance(dafunc(s), np.ndarray)
    tm.assert_numpy_array_equal(dafunc(s), npfunc(s))

    df = pd.DataFrame({'A': np.random.randint(1, 100, size=20),
                       'B': np.random.randint(1, 100, size=20),
                       'C': np.abs(np.random.randn(20))},
                      index=list('abcdefghijklmnopqrst'))
    ddf = dd.from_pandas(df, 3)

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(ddf), dd.DataFrame)
    # result may be read-only ndarray
    exp = pd.DataFrame(npfunc(df).copy(), columns=df.columns, index=df.index)
    assert_eq(dafunc(ddf), exp)

    assert isinstance(npfunc(ddf), np.ndarray)
    tm.assert_numpy_array_equal(npfunc(ddf), npfunc(df))

    assert isinstance(dafunc(df), np.ndarray)
    tm.assert_numpy_array_equal(dafunc(df), npfunc(df))


_UFUNCS_2ARG = ['logaddexp', 'logaddexp2', 'arctan2',
                'hypot', 'copysign', 'nextafter', 'ldexp',
                'fmod', 'logical_and', 'logical_or',
                'logical_xor', 'maximum', 'minimum',
                'fmax', 'fmin', 'greater',
                'greater_equal', 'less', 'less_equal',
                'not_equal', 'equal', 'logical_or',
                'logical_and', 'logical_xor']


@pytest.mark.parametrize('ufunc', _UFUNCS_2ARG)
@pytest.mark.parametrize('make_pandas_input', [
    lambda: pd.Series(np.random.randint(1, 100, size=20)),
    lambda: pd.DataFrame(np.random.randint(1, 100, size=(20, 2)),
                         columns=['A', 'B'])
])
def test_ufunc_with_2args(ufunc, make_pandas_input):

    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    pandas1 = make_pandas_input()
    pandas2 = make_pandas_input()

    dask1 = dd.from_pandas(pandas1, 3)
    dask2 = dd.from_pandas(pandas2, 4)

    pandas_type = pandas1.__class__
    dask_type = dask1.__class__

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(dask1, dask2), dask_type)
    assert_eq(dafunc(dask1, dask2), npfunc(pandas1, pandas2))

    # should be fine with pandas as a second arg, too
    assert isinstance(dafunc(dask1, pandas2), dask_type)
    assert_eq(dafunc(dask1, pandas2), npfunc(pandas1, pandas2))

    # applying NumPy ufunc is lazy
    if isinstance(npfunc, np.ufunc) and np.__version__ >= '1.13.0':
        assert isinstance(npfunc(dask1, dask2), dask_type)
        assert isinstance(npfunc(dask1, pandas2), dask_type)
    else:
        assert isinstance(npfunc(dask1, dask2), pandas_type)
        assert isinstance(npfunc(dask1, pandas2), pandas_type)

    assert_eq(npfunc(dask1, dask2), npfunc(pandas1, pandas2))
    assert_eq(npfunc(dask1, pandas2), npfunc(pandas1, pandas2))

    # applying Dask ufunc to normal Series triggers computation
    assert isinstance(dafunc(pandas1, pandas2), pandas_type)
    assert_eq(dafunc(pandas1, pandas2), npfunc(pandas1, pandas2))


@pytest.mark.parametrize('pandas,min,max', [
    (pd.Series(np.random.randint(1, 100, size=20)), 5, 50),
    (pd.DataFrame(np.random.randint(1, 100, size=(20, 2)),
                  columns=['A', 'B']), 5.5, 40.5)
])
def test_clip(pandas, min, max):

    dask = dd.from_pandas(pandas, 3)
    pandas_type = pandas.__class__
    dask_type = dask.__class__

    # clip internally calls dd.Series.clip

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(da.clip(dask, min, max), dask_type)
    assert_eq(da.clip(dask, min, max), np.clip(pandas, min, max))

    # applying Numpy ufunc doesn't trigger computation
    assert isinstance(np.clip(dask, min, max), dask_type)
    assert_eq(np.clip(dask, min, max), np.clip(pandas, min, max))

    # applying Dask ufunc to normal pandas objects triggers computation
    assert isinstance(da.clip(pandas, min, max), pandas_type)
    assert_eq(da.clip(pandas, min, max), np.clip(pandas, min, max))


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
@pytest.mark.parametrize('ufunc', _BASE_UFUNCS)
def test_frame_ufunc_out(ufunc):
    npfunc = getattr(np, ufunc)
    dafunc = getattr(da, ufunc)

    input_matrix = np.random.randint(1, 100, size=(20, 2))

    df = pd.DataFrame(input_matrix, columns=['A', 'B'])
    ddf = dd.from_pandas(df, 3)
    df_out = pd.DataFrame(np.random.randint(1, 100, size=(20, 2)),
                          columns=['Y', 'Z'])
    ddf_out_np = dd.from_pandas(df_out, 3)
    ddf_out_da = dd.from_pandas(df_out, 3)

    with pytest.warns(None):
        npfunc(ddf, out=ddf_out_np)
        dafunc(ddf, out=ddf_out_da)
        assert_eq(ddf_out_np, ddf_out_da)

    with pytest.warns(None):
        expected = pd.DataFrame(npfunc(input_matrix), columns=['A', 'B'])
        assert_eq(ddf_out_np, expected)


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
def test_frame_2ufunc_out():
    input_matrix = np.random.randint(1, 100, size=(20, 2))

    df = pd.DataFrame(input_matrix, columns=['A', 'B'])
    ddf = dd.from_pandas(df, 3)

    # column number mismatch
    df_out = pd.DataFrame(np.random.randint(1, 100, size=(20, 3)),
                          columns=['X', 'Y', 'Z'])
    ddf_out = dd.from_pandas(df_out, 3)

    with pytest.raises(ValueError):
        np.sin(ddf, out=ddf_out)

    # types mismatch
    ddf_out = dd.from_pandas(pd.Series([0]),1)
    with pytest.raises(TypeError):
        np.sin(ddf, out=ddf_out)

    df_out = pd.DataFrame(np.random.randint(1, 100, size=(20, 2)),
                          columns=['X', 'Y'])
    ddf_out = dd.from_pandas(df_out, 3)

    np.sin(ddf, out=ddf_out)
    np.add(ddf_out, 10, out=ddf_out)

    expected = pd.DataFrame(np.sin(input_matrix) + 10, columns=['A', 'B'])

    assert_eq(ddf_out, expected)


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
@pytest.mark.parametrize('arg1', [
    pd.Series(np.abs(np.random.randn(100))),
    pd.DataFrame({'A': np.random.randint(1, 100, size=20),
                  'B': np.random.randint(1, 100, size=20),
                  'C': np.abs(np.random.randn(20))})])
@pytest.mark.parametrize('arg2', [2, dd.from_pandas(pd.Series([0]), 1).sum()])
@pytest.mark.parametrize('ufunc', _UFUNCS_2ARG)
def test_mixed_types(ufunc, arg1, arg2):
    npfunc = getattr(np, ufunc)
    dafunc = getattr(da, ufunc)

    dask = dd.from_pandas(arg1, 3)

    pandas_type = arg1.__class__
    dask_type = dask.__class__

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(dask, arg2), dask_type)
    assert_eq(dafunc(dask, arg2), npfunc(dask, arg2))

    # applying NumPy ufunc is lazy
    assert isinstance(npfunc(dask, arg2), dask_type)
    assert_eq(npfunc(dask, arg2), npfunc(arg1, arg2))

    # applying Dask ufunc to normal Series triggers computation
    assert isinstance(dafunc(arg1, arg2), pandas_type)
    assert_eq(dafunc(arg1, arg2), npfunc(arg1, arg2))

    # swapping arguments

    # first parameter of ldexp should be array-like
    if ufunc == 'ldexp':
        return

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(arg2, dask), dask_type)
    assert_eq(dafunc(arg2, dask), npfunc(arg2, dask))

    # applying NumPy ufunc is lazy
    assert isinstance(npfunc(arg2, dask), dask_type)
    assert_eq(npfunc(arg2, dask), npfunc(arg2, dask))

    # applying Dask ufunc to normal Series triggers computation
    assert isinstance(dafunc(arg2, arg1), pandas_type)
    assert_eq(dafunc(arg2, arg1), npfunc(arg2, arg1))


@pytest.mark.skipif(np.__version__ < '1.13.0', reason='array_ufunc not present')
@pytest.mark.parametrize('ufunc', _UFUNCS_2ARG)
@pytest.mark.parametrize('pandas,darray',
                         [(pd.Series(np.random.randint(1, 100, size=(100,))),
                           da.from_array(np.random.randint(1, 100, size=(100,)),
                                         chunks=(50,))),
                          (pd.DataFrame(np.random.randint(1, 100, size=(20, 2)),
                                        columns=['A', 'B']),
                           da.from_array(np.random.randint(1, 100, size=(20, 2)),
                                         chunks=(10, 2)))])
def test_2args_with_array(ufunc, pandas, darray):
    dafunc = getattr(da, ufunc)
    npfunc = getattr(np, ufunc)

    dask = dd.from_pandas(pandas, 2)
    dask_type = dask.__class__

    # applying Dask ufunc doesn't trigger computation
    assert isinstance(dafunc(dask, darray), dask_type)
    assert isinstance(dafunc(darray, dask), dask_type)

    tm.assert_numpy_array_equal(dafunc(dask, darray).compute().values,
                                npfunc(pandas.values, darray).compute())

    # applying NumPy ufunc is lazy
    assert isinstance(npfunc(dask, darray), dask_type)
    assert isinstance(npfunc(darray, dask), dask_type)

    tm.assert_numpy_array_equal(npfunc(dask, darray).compute().values,
                                npfunc(pandas.values, darray.compute()))
    tm.assert_numpy_array_equal(npfunc(darray, dask).compute().values,
                                npfunc(darray.compute(), pandas.values))


@pytest.mark.parametrize('redfunc', ['sum', 'prod', 'min', 'max', 'mean'])
@pytest.mark.parametrize('ufunc', _BASE_UFUNCS)
@pytest.mark.parametrize('pandas',
                         [pd.Series(np.abs(np.random.randn(100))),
                          pd.DataFrame({'A': np.random.randint(1, 100, size=20),
                                        'B': np.random.randint(1, 100, size=20),
                                        'C': np.abs(np.random.randn(20))})])
def test_ufunc_with_reduction(redfunc, ufunc, pandas):
    dask = dd.from_pandas(pandas, 3)

    np_redfunc = getattr(np, redfunc)
    np_ufunc = getattr(np, ufunc)

    with pytest.warns(None):
        assert isinstance(np_redfunc(dask), (dd.DataFrame, dd.Series, dd.core.Scalar))
        assert_eq(np_redfunc(np_ufunc(dask)), np_redfunc(np_ufunc(pandas)))


@pytest.mark.parametrize('pandas',
                         [pd.Series(np.random.randint(1, 100, size=100)),
                          pd.DataFrame({'A': np.random.randint(1, 100, size=20),
                                        'B': np.random.randint(1, 100, size=20),
                                        'C': np.abs(np.random.randn(20))})])
@pytest.mark.parametrize('scalar', [15, 16.4, np.int64(15), np.float64(16.4)])
def test_ufunc_numpy_scalar_comparison(pandas, scalar):
    # Regression test for issue #3392

    dask_compare = scalar >= dd.from_pandas(pandas, npartitions=3)
    pandas_compare = scalar >= pandas

    assert_eq(dask_compare, pandas_compare)
