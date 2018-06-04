import numpy as np
import pandas as pd
import pandas.util.testing as tm

import pytest
from threading import Lock
from multiprocessing.pool import ThreadPool

import dask.array as da
import dask.dataframe as dd
from dask.dataframe.io.io import _meta_from_array
from dask.delayed import Delayed, delayed

from dask.utils import tmpfile
from dask.local import get_sync

from dask.dataframe.utils import assert_eq, is_categorical_dtype


####################
# Arrays and BColz #
####################


def test_meta_from_array():
    x = np.array([[1, 2], [3, 4]], dtype=np.int64)
    res = _meta_from_array(x)
    assert isinstance(res, pd.DataFrame)
    assert res[0].dtype == np.int64
    assert res[1].dtype == np.int64
    tm.assert_index_equal(res.columns, pd.Index([0, 1]))

    x = np.array([[1., 2.], [3., 4.]], dtype=np.float64)
    res = _meta_from_array(x, columns=['a', 'b'])
    assert isinstance(res, pd.DataFrame)
    assert res['a'].dtype == np.float64
    assert res['b'].dtype == np.float64
    tm.assert_index_equal(res.columns, pd.Index(['a', 'b']))

    with pytest.raises(ValueError):
        _meta_from_array(x, columns=['a', 'b', 'c'])

    np.random.seed(42)
    x = np.random.rand(201, 2)
    x = dd.from_array(x, chunksize=50, columns=['a', 'b'])
    assert len(x.divisions) == 6   # Should be 5 partitions and the end


def test_meta_from_1darray():
    x = np.array([1., 2., 3.], dtype=np.float64)
    res = _meta_from_array(x)
    assert isinstance(res, pd.Series)
    assert res.dtype == np.float64

    x = np.array([1, 2, 3], dtype=np.object_)
    res = _meta_from_array(x, columns='x')
    assert isinstance(res, pd.Series)
    assert res.name == 'x'
    assert res.dtype == np.object_

    x = np.array([1, 2, 3], dtype=np.object_)
    res = _meta_from_array(x, columns=['x'])
    assert isinstance(res, pd.DataFrame)
    assert res['x'].dtype == np.object_
    tm.assert_index_equal(res.columns, pd.Index(['x']))

    with pytest.raises(ValueError):
        _meta_from_array(x, columns=['a', 'b'])


def test_meta_from_recarray():
    x = np.array([(i, i * 10) for i in range(10)],
                 dtype=[('a', np.float64), ('b', np.int64)])
    res = _meta_from_array(x)
    assert isinstance(res, pd.DataFrame)
    assert res['a'].dtype == np.float64
    assert res['b'].dtype == np.int64
    tm.assert_index_equal(res.columns, pd.Index(['a', 'b']))

    res = _meta_from_array(x, columns=['b', 'a'])
    assert isinstance(res, pd.DataFrame)
    assert res['a'].dtype == np.float64
    assert res['b'].dtype == np.int64
    tm.assert_index_equal(res.columns, pd.Index(['b', 'a']))

    with pytest.raises(ValueError):
        _meta_from_array(x, columns=['a', 'b', 'c'])


def test_from_array():
    x = np.arange(10 * 3).reshape(10, 3)
    d = dd.from_array(x, chunksize=4)
    assert isinstance(d, dd.DataFrame)
    tm.assert_index_equal(d.columns, pd.Index([0, 1, 2]))
    assert d.divisions == (0, 4, 8, 9)
    assert (d.compute().values == x).all()

    d = dd.from_array(x, chunksize=4, columns=list('abc'))
    assert isinstance(d, dd.DataFrame)
    tm.assert_index_equal(d.columns, pd.Index(['a', 'b', 'c']))
    assert d.divisions == (0, 4, 8, 9)
    assert (d.compute().values == x).all()

    with pytest.raises(ValueError):
        dd.from_array(np.ones(shape=(10, 10, 10)))


def test_from_array_with_record_dtype():
    x = np.array([(i, i * 10) for i in range(10)],
                 dtype=[('a', 'i4'), ('b', 'i4')])
    d = dd.from_array(x, chunksize=4)
    assert isinstance(d, dd.DataFrame)
    assert list(d.columns) == ['a', 'b']
    assert d.divisions == (0, 4, 8, 9)

    assert (d.compute().to_records(index=False) == x).all()


def test_from_bcolz_multiple_threads():
    bcolz = pytest.importorskip('bcolz')
    pool = ThreadPool(processes=5)

    def check(i):
        t = bcolz.ctable([[1, 2, 3], [1., 2., 3.], ['a', 'b', 'a']],
                         names=['x', 'y', 'a'])
        d = dd.from_bcolz(t, chunksize=2)
        assert d.npartitions == 2
        assert is_categorical_dtype(d.dtypes['a'])
        assert list(d.x.compute(get=get_sync)) == [1, 2, 3]
        assert list(d.a.compute(get=get_sync)) == ['a', 'b', 'a']

        d = dd.from_bcolz(t, chunksize=2, index='x')
        L = list(d.index.compute(get=get_sync))
        assert L == [1, 2, 3] or L == [1, 3, 2]

        # Names
        assert (sorted(dd.from_bcolz(t, chunksize=2).dask) ==
                sorted(dd.from_bcolz(t, chunksize=2).dask))
        assert (sorted(dd.from_bcolz(t, chunksize=2).dask) !=
                sorted(dd.from_bcolz(t, chunksize=3).dask))

    pool.map(check, range(5))


def test_from_bcolz():
    bcolz = pytest.importorskip('bcolz')

    t = bcolz.ctable([[1, 2, 3], [1., 2., 3.], ['a', 'b', 'a']],
                     names=['x', 'y', 'a'])
    d = dd.from_bcolz(t, chunksize=2)
    assert d.npartitions == 2
    assert is_categorical_dtype(d.dtypes['a'])
    assert list(d.x.compute(get=get_sync)) == [1, 2, 3]
    assert list(d.a.compute(get=get_sync)) == ['a', 'b', 'a']
    L = list(d.index.compute(get=get_sync))
    assert L == [0, 1, 2]

    d = dd.from_bcolz(t, chunksize=2, index='x')
    L = list(d.index.compute(get=get_sync))
    assert L == [1, 2, 3] or L == [1, 3, 2]

    # Names
    assert (sorted(dd.from_bcolz(t, chunksize=2).dask) ==
            sorted(dd.from_bcolz(t, chunksize=2).dask))
    assert (sorted(dd.from_bcolz(t, chunksize=2).dask) !=
            sorted(dd.from_bcolz(t, chunksize=3).dask))

    dsk = dd.from_bcolz(t, chunksize=3).dask

    t.append((4, 4., 'b'))
    t.flush()

    assert (sorted(dd.from_bcolz(t, chunksize=2).dask) !=
            sorted(dsk))


def test_from_bcolz_no_lock():
    bcolz = pytest.importorskip('bcolz')
    locktype = type(Lock())

    t = bcolz.ctable([[1, 2, 3], [1., 2., 3.], ['a', 'b', 'a']],
                     names=['x', 'y', 'a'], chunklen=2)
    a = dd.from_bcolz(t, chunksize=2)
    b = dd.from_bcolz(t, chunksize=2, lock=True)
    c = dd.from_bcolz(t, chunksize=2, lock=False)
    assert_eq(a, b)
    assert_eq(a, c)

    assert not any(isinstance(item, locktype)
                   for v in c.dask.values()
                   for item in v)


def test_from_bcolz_filename():
    bcolz = pytest.importorskip('bcolz')

    with tmpfile('.bcolz') as fn:
        t = bcolz.ctable([[1, 2, 3], [1., 2., 3.], ['a', 'b', 'a']],
                         names=['x', 'y', 'a'],
                         rootdir=fn)
        t.flush()

        d = dd.from_bcolz(fn, chunksize=2)
        assert list(d.x.compute()) == [1, 2, 3]


def test_from_bcolz_column_order():
    bcolz = pytest.importorskip('bcolz')

    t = bcolz.ctable([[1, 2, 3], [1., 2., 3.], ['a', 'b', 'a']],
                     names=['x', 'y', 'a'])
    df = dd.from_bcolz(t, chunksize=2)
    assert list(df.loc[0].compute().columns) == ['x', 'y', 'a']


def test_from_pandas_dataframe():
    a = list('aaaaaaabbbbbbbbccccccc')
    df = pd.DataFrame(dict(a=a, b=np.random.randn(len(a))),
                      index=pd.date_range(start='20120101', periods=len(a)))
    ddf = dd.from_pandas(df, 3)
    assert len(ddf.dask) == 3
    assert len(ddf.divisions) == len(ddf.dask) + 1
    assert isinstance(ddf.divisions[0], type(df.index[0]))
    tm.assert_frame_equal(df, ddf.compute())
    ddf = dd.from_pandas(df, chunksize=8)
    msg = 'Exactly one of npartitions and chunksize must be specified.'
    with pytest.raises(ValueError) as err:
        dd.from_pandas(df, npartitions=2, chunksize=2)
    assert msg in str(err.value)
    with pytest.raises((ValueError, AssertionError)) as err:
        dd.from_pandas(df)
    assert msg in str(err.value)
    assert len(ddf.dask) == 3
    assert len(ddf.divisions) == len(ddf.dask) + 1
    assert isinstance(ddf.divisions[0], type(df.index[0]))
    tm.assert_frame_equal(df, ddf.compute())


def test_from_pandas_small():
    df = pd.DataFrame({'x': [1, 2, 3]})
    for i in [1, 2, 30]:
        a = dd.from_pandas(df, i)
        assert len(a.compute()) == 3
        assert a.divisions[0] == 0
        assert a.divisions[-1] == 2

        a = dd.from_pandas(df, chunksize=i)
        assert len(a.compute()) == 3
        assert a.divisions[0] == 0
        assert a.divisions[-1] == 2

    for sort in [True, False]:
        for i in [0, 2]:
            df = pd.DataFrame({'x': [0] * i})
            ddf = dd.from_pandas(df, npartitions=5, sort=sort)
            assert_eq(df, ddf)

            s = pd.Series([0] * i, name='x')
            ds = dd.from_pandas(s, npartitions=5, sort=sort)
            assert_eq(s, ds)


@pytest.mark.xfail(reason="")
def test_from_pandas_npartitions_is_accurate():
    df = pd.DataFrame({'x': [1, 2, 3, 4, 5, 6], 'y': list('abdabd')},
                      index=[10, 20, 30, 40, 50, 60])
    for n in [1, 2, 4, 5]:
        assert dd.from_pandas(df, npartitions=n).npartitions == n


def test_from_pandas_series():
    n = 20
    s = pd.Series(np.random.randn(n),
                  index=pd.date_range(start='20120101', periods=n))
    ds = dd.from_pandas(s, 3)
    assert len(ds.dask) == 3
    assert len(ds.divisions) == len(ds.dask) + 1
    assert isinstance(ds.divisions[0], type(s.index[0]))
    tm.assert_series_equal(s, ds.compute())

    ds = dd.from_pandas(s, chunksize=8)
    assert len(ds.dask) == 3
    assert len(ds.divisions) == len(ds.dask) + 1
    assert isinstance(ds.divisions[0], type(s.index[0]))
    tm.assert_series_equal(s, ds.compute())


def test_from_pandas_non_sorted():
    df = pd.DataFrame({'x': [1, 2, 3]}, index=[3, 1, 2])
    ddf = dd.from_pandas(df, npartitions=2, sort=False)
    assert not ddf.known_divisions
    assert_eq(df, ddf)

    ddf = dd.from_pandas(df, chunksize=2, sort=False)
    assert not ddf.known_divisions
    assert_eq(df, ddf)


def test_from_pandas_single_row():
    df = pd.DataFrame({'x': [1]}, index=[1])
    ddf = dd.from_pandas(df, npartitions=1)
    assert ddf.divisions == (1, 1)
    assert_eq(ddf, df)


@pytest.mark.skipif(np.__version__ < '1.11',
                    reason='datetime unit unsupported in NumPy < 1.11')
def test_from_pandas_with_datetime_index():
    df = pd.DataFrame({"Date": ["2015-08-28", "2015-08-27", "2015-08-26",
                                "2015-08-25", "2015-08-24", "2015-08-21",
                                "2015-08-20", "2015-08-19", "2015-08-18"],
                       "Val": list(range(9))})
    df.Date = df.Date.astype('datetime64[ns]')
    ddf = dd.from_pandas(df, 2)
    assert_eq(df, ddf)
    ddf = dd.from_pandas(df, chunksize=2)
    assert_eq(df, ddf)


def test_DataFrame_from_dask_array():
    x = da.ones((10, 3), chunks=(4, 2))

    df = dd.from_dask_array(x, ['a', 'b', 'c'])
    assert isinstance(df, dd.DataFrame)
    tm.assert_index_equal(df.columns, pd.Index(['a', 'b', 'c']))
    assert list(df.divisions) == [0, 4, 8, 9]
    assert (df.compute(get=get_sync).values == x.compute(get=get_sync)).all()

    # dd.from_array should re-route to from_dask_array
    df2 = dd.from_array(x, columns=['a', 'b', 'c'])
    assert isinstance(df, dd.DataFrame)
    tm.assert_index_equal(df2.columns, df.columns)
    assert df2.divisions == df.divisions


def test_Series_from_dask_array():
    x = da.ones(10, chunks=4)

    ser = dd.from_dask_array(x, 'a')
    assert isinstance(ser, dd.Series)
    assert ser.name == 'a'
    assert list(ser.divisions) == [0, 4, 8, 9]
    assert (ser.compute(get=get_sync).values == x.compute(get=get_sync)).all()

    ser = dd.from_dask_array(x)
    assert isinstance(ser, dd.Series)
    assert ser.name is None

    # dd.from_array should re-route to from_dask_array
    ser2 = dd.from_array(x)
    assert isinstance(ser2, dd.Series)
    assert_eq(ser, ser2)


def test_from_dask_array_compat_numpy_array():
    x = da.ones((3, 3, 3), chunks=2)

    with pytest.raises(ValueError):
        dd.from_dask_array(x)       # dask

    with pytest.raises(ValueError):
        dd.from_array(x.compute())  # numpy

    x = da.ones((10, 3), chunks=(3, 3))
    d1 = dd.from_dask_array(x)       # dask
    assert isinstance(d1, dd.DataFrame)
    assert (d1.compute().values == x.compute()).all()
    tm.assert_index_equal(d1.columns, pd.Index([0, 1, 2]))

    d2 = dd.from_array(x.compute())  # numpy
    assert isinstance(d1, dd.DataFrame)
    assert (d2.compute().values == x.compute()).all()
    tm.assert_index_equal(d2.columns, pd.Index([0, 1, 2]))

    with pytest.raises(ValueError):
        dd.from_dask_array(x, columns=['a'])       # dask

    with pytest.raises(ValueError):
        dd.from_array(x.compute(), columns=['a'])  # numpy

    d1 = dd.from_dask_array(x, columns=['a', 'b', 'c'])       # dask
    assert isinstance(d1, dd.DataFrame)
    assert (d1.compute().values == x.compute()).all()
    tm.assert_index_equal(d1.columns, pd.Index(['a', 'b', 'c']))

    d2 = dd.from_array(x.compute(), columns=['a', 'b', 'c'])  # numpy
    assert isinstance(d1, dd.DataFrame)
    assert (d2.compute().values == x.compute()).all()
    tm.assert_index_equal(d2.columns, pd.Index(['a', 'b', 'c']))


def test_from_dask_array_compat_numpy_array_1d():

    x = da.ones(10, chunks=3)
    d1 = dd.from_dask_array(x)       # dask
    assert isinstance(d1, dd.Series)
    assert (d1.compute().values == x.compute()).all()
    assert d1.name is None

    d2 = dd.from_array(x.compute())  # numpy
    assert isinstance(d1, dd.Series)
    assert (d2.compute().values == x.compute()).all()
    assert d2.name is None

    d1 = dd.from_dask_array(x, columns='name')       # dask
    assert isinstance(d1, dd.Series)
    assert (d1.compute().values == x.compute()).all()
    assert d1.name == 'name'

    d2 = dd.from_array(x.compute(), columns='name')  # numpy
    assert isinstance(d1, dd.Series)
    assert (d2.compute().values == x.compute()).all()
    assert d2.name == 'name'

    # passing list via columns results in DataFrame
    d1 = dd.from_dask_array(x, columns=['name'])       # dask
    assert isinstance(d1, dd.DataFrame)
    assert (d1.compute().values == x.compute()).all()
    tm.assert_index_equal(d1.columns, pd.Index(['name']))

    d2 = dd.from_array(x.compute(), columns=['name'])  # numpy
    assert isinstance(d1, dd.DataFrame)
    assert (d2.compute().values == x.compute()).all()
    tm.assert_index_equal(d2.columns, pd.Index(['name']))


def test_from_dask_array_struct_dtype():
    x = np.array([(1, 'a'), (2, 'b')], dtype=[('a', 'i4'), ('b', 'object')])
    y = da.from_array(x, chunks=(1,))
    df = dd.from_dask_array(y)
    tm.assert_index_equal(df.columns, pd.Index(['a', 'b']))
    assert_eq(df, pd.DataFrame(x))

    assert_eq(dd.from_dask_array(y, columns=['b', 'a']),
              pd.DataFrame(x, columns=['b', 'a']))


def test_from_dask_array_unknown_chunks():
    # Series
    dx = da.Array({('x', 0): np.arange(5), ('x', 1): np.arange(5, 11)}, 'x',
                  ((np.nan, np.nan,),), np.arange(1).dtype)
    df = dd.from_dask_array(dx)
    assert isinstance(df, dd.Series)
    assert not df.known_divisions
    assert_eq(df, pd.Series(np.arange(11)), check_index=False)

    # DataFrame
    dsk = {('x', 0, 0): np.random.random((2, 3)),
           ('x', 1, 0): np.random.random((5, 3))}
    dx = da.Array(dsk, 'x', ((np.nan, np.nan,), (3,)), np.float64)
    df = dd.from_dask_array(dx)
    assert isinstance(df, dd.DataFrame)
    assert not df.known_divisions
    assert_eq(df, pd.DataFrame(dx.compute()), check_index=False)

    # Unknown width
    dx = da.Array(dsk, 'x', ((np.nan, np.nan,), (np.nan,)), np.float64)
    with pytest.raises(ValueError):
        df = dd.from_dask_array(dx)


def test_to_bag():
    pytest.importorskip('dask.bag')
    a = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                      'y': [2, 3, 4, 5]},
                     index=pd.Index([1., 2., 3., 4.], name='ind'))
    ddf = dd.from_pandas(a, 2)

    assert ddf.to_bag().compute() == list(a.itertuples(False))
    assert ddf.to_bag(True).compute() == list(a.itertuples(True))
    assert ddf.x.to_bag(True).compute() == list(a.x.iteritems())
    assert ddf.x.to_bag().compute() == list(a.x)


def test_to_records():
    pytest.importorskip('dask.array')
    from dask.array.utils import assert_eq
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [2, 3, 4, 5]},
                      index=pd.Index([1., 2., 3., 4.], name='ind'))
    ddf = dd.from_pandas(df, 2)

    assert_eq(df.to_records(), ddf.to_records())


def test_from_delayed():
    df = pd.DataFrame(data=np.random.normal(size=(10, 4)), columns=list('abcd'))
    parts = [df.iloc[:1], df.iloc[1:3], df.iloc[3:6], df.iloc[6:10]]
    dfs = [delayed(parts.__getitem__)(i) for i in range(4)]
    meta = dfs[0].compute()

    my_len = lambda x: pd.Series([len(x)])

    for divisions in [None, [0, 1, 3, 6, 10]]:
        ddf = dd.from_delayed(dfs, meta=meta, divisions=divisions)
        assert_eq(ddf, df)
        assert list(ddf.map_partitions(my_len).compute()) == [1, 2, 3, 4]
        assert ddf.known_divisions == (divisions is not None)

        s = dd.from_delayed([d.a for d in dfs], meta=meta.a,
                            divisions=divisions)
        assert_eq(s, df.a)
        assert list(s.map_partitions(my_len).compute()) == [1, 2, 3, 4]
        assert ddf.known_divisions == (divisions is not None)

    meta2 = [(c, 'f8') for c in df.columns]
    assert_eq(dd.from_delayed(dfs, meta=meta2), df)
    assert_eq(dd.from_delayed([d.a for d in dfs], meta=('a', 'f8')), df.a)

    with pytest.raises(ValueError):
        dd.from_delayed(dfs, meta=meta, divisions=[0, 1, 3, 6])

    with pytest.raises(ValueError) as e:
        dd.from_delayed(dfs, meta=meta.a).compute()
    assert str(e.value).startswith('Metadata mismatch found in `from_delayed`')


def test_from_delayed_sorted():
    a = pd.DataFrame({'x': [1, 2]}, index=[1, 10])
    b = pd.DataFrame({'x': [4, 1]}, index=[100, 200])

    A = dd.from_delayed([delayed(a), delayed(b)], divisions='sorted')
    assert A.known_divisions

    assert A.divisions == (1, 100, 200)


def test_to_delayed():
    df = pd.DataFrame({'x': [1, 2, 3, 4], 'y': [10, 20, 30, 40]})
    ddf = dd.from_pandas(df, npartitions=2)

    # Frame
    a, b = ddf.to_delayed()
    assert isinstance(a, Delayed)
    assert isinstance(b, Delayed)
    assert_eq(a.compute(), df.iloc[:2])

    # Scalar
    x = ddf.x.sum()
    dx = x.to_delayed()
    assert isinstance(dx, Delayed)
    assert_eq(dx.compute(), x)


def test_to_delayed_optimize_graph():
    df = pd.DataFrame({'x': list(range(20))})
    ddf = dd.from_pandas(df, npartitions=20)
    ddf2 = (ddf + 1).loc[:2]

    # Frame
    d = ddf2.to_delayed()[0]
    assert len(d.dask) < 20
    d2 = ddf2.to_delayed(optimize_graph=False)[0]
    assert sorted(d2.dask) == sorted(ddf2.dask)
    assert_eq(ddf2.get_partition(0), d.compute())
    assert_eq(ddf2.get_partition(0), d2.compute())

    # Scalar
    x = ddf2.x.sum()
    dx = x.to_delayed()
    dx2 = x.to_delayed(optimize_graph=False)
    assert len(dx.dask) < len(dx2.dask)
    assert_eq(dx.compute(), dx2.compute())
