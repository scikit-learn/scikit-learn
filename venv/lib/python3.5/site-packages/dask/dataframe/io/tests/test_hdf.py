import numpy as np
import pandas as pd
import pandas.util.testing as tm

import sys
import os
import dask
import pytest

from time import sleep

import dask.dataframe as dd

from dask.utils import tmpfile, tmpdir, dependency_depth

from dask.dataframe.utils import assert_eq


def test_to_hdf():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    a = dd.from_pandas(df, 2)

    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data')
        out = pd.read_hdf(fn, '/data')
        tm.assert_frame_equal(df, out[:])

    with tmpfile('h5') as fn:
        a.x.to_hdf(fn, '/data')
        out = pd.read_hdf(fn, '/data')
        tm.assert_series_equal(df.x, out[:])

    a = dd.from_pandas(df, 1)
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data')
        out = pd.read_hdf(fn, '/data')
        tm.assert_frame_equal(df, out[:])

    # test compute = False
    with tmpfile('h5') as fn:
        r = a.to_hdf(fn, '/data', compute=False)
        r.compute()
        out = pd.read_hdf(fn, '/data')
        tm.assert_frame_equal(df, out[:])


def test_to_hdf_multiple_nodes():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    a = dd.from_pandas(df, 2)
    df16 = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                               'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                         'y': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
                               14, 15, 16]},
                        index=[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
                               12., 13., 14., 15., 16.])
    b = dd.from_pandas(df16, 16)

    # saving to multiple nodes
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data*')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df, out)

    # saving to multiple nodes making sure order is kept
    with tmpfile('h5') as fn:
        b.to_hdf(fn, '/data*')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df16, out)

    # saving to multiple datasets with custom name_function
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data_*', name_function=lambda i: 'a' * (i + 1))
        out = dd.read_hdf(fn, '/data_*')
        assert_eq(df, out)

        out = pd.read_hdf(fn, '/data_a')
        tm.assert_frame_equal(out, df.iloc[:2])
        out = pd.read_hdf(fn, '/data_aa')
        tm.assert_frame_equal(out, df.iloc[2:])

    # test multiple nodes with hdf object
    with tmpfile('h5') as fn:
        with pd.HDFStore(fn) as hdf:
            b.to_hdf(hdf, '/data*')
            out = dd.read_hdf(fn, '/data*')
            assert_eq(df16, out)


def test_to_hdf_multiple_files():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    a = dd.from_pandas(df, 2)
    df16 = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                               'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                         'y': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                               15, 16]},
                        index=[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
                               12., 13., 14., 15., 16.])
    b = dd.from_pandas(df16, 16)

    # saving to multiple files
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data_*.h5')
        a.to_hdf(fn, '/data')
        out = dd.read_hdf(fn, '/data')
        assert_eq(df, out)

    # saving to multiple files making sure order is kept
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data_*.h5')
        b.to_hdf(fn, '/data')
        out = dd.read_hdf(fn, '/data')
        assert_eq(df16, out)

    # saving to multiple files with custom name_function
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data_*.h5')
        a.to_hdf(fn, '/data', name_function=lambda i: 'a' * (i + 1))
        out = dd.read_hdf(fn, '/data')
        assert_eq(df, out)

        out = pd.read_hdf(os.path.join(dn, 'data_a.h5'), '/data')
        tm.assert_frame_equal(out, df.iloc[:2])
        out = pd.read_hdf(os.path.join(dn, 'data_aa.h5'), '/data')
        tm.assert_frame_equal(out, df.iloc[2:])

    # test hdf object
    with tmpfile('h5') as fn:
        with pd.HDFStore(fn) as hdf:
            a.to_hdf(hdf, '/data*')
            out = dd.read_hdf(fn, '/data*')
            assert_eq(df, out)


def test_to_hdf_modes_multiple_nodes():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])

    # appending a single partition to existing data
    a = dd.from_pandas(df, 1)
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data2')
        a.to_hdf(fn, '/data*', mode='a')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df.append(df), out)

    # overwriting a file with a single partition
    a = dd.from_pandas(df, 1)
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data2')
        a.to_hdf(fn, '/data*', mode='w')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df, out)

    # appending two partitions to existing data
    a = dd.from_pandas(df, 2)
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data2')
        a.to_hdf(fn, '/data*', mode='a')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df.append(df), out)

    # overwriting a file with two partitions
    a = dd.from_pandas(df, 2)
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data2')
        a.to_hdf(fn, '/data*', mode='w')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df, out)

    # overwriting a single partition, keeping other partitions
    a = dd.from_pandas(df, 2)
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data1')
        a.to_hdf(fn, '/data2')
        a.to_hdf(fn, '/data*', mode='a', append=False)
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df.append(df), out)


def test_to_hdf_modes_multiple_files():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])

    # appending a single partition to existing data
    a = dd.from_pandas(df, 1)
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data*')
        a.to_hdf(os.path.join(dn, 'data2'), '/data')
        a.to_hdf(fn, '/data', mode='a')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df.append(df), out)

    # appending two partitions to existing data
    a = dd.from_pandas(df, 2)
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data*')
        a.to_hdf(os.path.join(dn, 'data2'), '/data')
        a.to_hdf(fn, '/data', mode='a')
        out = dd.read_hdf(fn, '/data')
        assert_eq(df.append(df), out)

    # overwriting a file with two partitions
    a = dd.from_pandas(df, 2)
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data*')
        a.to_hdf(os.path.join(dn, 'data1'), '/data')
        a.to_hdf(fn, '/data', mode='w')
        out = dd.read_hdf(fn, '/data')
        assert_eq(df, out)

    # overwriting a single partition, keeping other partitions
    a = dd.from_pandas(df, 2)
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data*')
        a.to_hdf(os.path.join(dn, 'data1'), '/data')
        a.to_hdf(fn, '/data', mode='a', append=False)
        out = dd.read_hdf(fn, '/data')
        assert_eq(df.append(df), out)


def test_to_hdf_link_optimizations():
    """testing dask link levels is correct by calculating the depth of the dask graph"""
    pytest.importorskip('tables')
    df16 = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                               'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                         'y': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                               15, 16]},
                        index=[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
                               12., 13., 14., 15., 16.])
    a = dd.from_pandas(df16, 16)

    # saving to multiple hdf files, no links are needed
    # expected layers: from_pandas, to_hdf, list = depth of 3
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data*')
        d = a.to_hdf(fn, '/data', compute=False)
        assert dependency_depth(d.dask) == 3

    # saving to a single hdf file with multiple nodes
    # all subsequent nodes depend on the first
    # expected layers: from_pandas, first to_hdf(creates file+node), subsequent to_hdfs, list = 4
    with tmpfile() as fn:
        d = a.to_hdf(fn, '/data*', compute=False)
        assert dependency_depth(d.dask) == 4

    # saving to a single hdf file with a single node
    # every node depends on the previous node
    # expected layers: from_pandas, to_hdf times npartitions(15), list = 2 + npartitions = 17
    with tmpfile() as fn:
        d = a.to_hdf(fn, '/data', compute=False)
        assert dependency_depth(d.dask) == 2 + a.npartitions


@pytest.mark.slow
def test_to_hdf_lock_delays():
    pytest.importorskip('tables')
    df16 = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i',
                               'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                         'y': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
                               15, 16]},
                        index=[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11.,
                               12., 13., 14., 15., 16.])
    a = dd.from_pandas(df16, 16)

    # adding artifichial delays to make sure last tasks finish first
    # that's a way to simulate last tasks finishing last
    def delayed_nop(i):
        if i[1] < 10:
            sleep(0.1 * (10 - i[1]))
        return i

    # saving to multiple hdf nodes
    with tmpfile() as fn:
        a = a.apply(delayed_nop, axis=1, meta=a)
        a.to_hdf(fn, '/data*')
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df16, out)

    # saving to multiple hdf files
    # adding artifichial delays to make sure last tasks finish first
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data*')
        a = a.apply(delayed_nop, axis=1, meta=a)
        a.to_hdf(fn, '/data')
        out = dd.read_hdf(fn, '/data')
        assert_eq(df16, out)


def test_to_hdf_exceptions():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    a = dd.from_pandas(df, 1)

    # triggering too many asterisks error
    with tmpdir() as dn:
        with pytest.raises(ValueError):
            fn = os.path.join(dn, 'data_*.h5')
            a.to_hdf(fn, '/data_*')

    # triggering too many asterisks error
    with tmpfile() as fn:
        with pd.HDFStore(fn) as hdf:
            with pytest.raises(ValueError):
                a.to_hdf(hdf, '/data_*_*')


@pytest.mark.parametrize('get', [dask.get,
                                 dask.threaded.get,
                                 dask.multiprocessing.get])
@pytest.mark.parametrize('npartitions', [1, 4, 10])
def test_to_hdf_schedulers(get, npartitions):
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                       'y': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]},
                      index=[1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16.])
    a = dd.from_pandas(df, npartitions=npartitions)

    # test single file single node
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data', get=get)
        out = pd.read_hdf(fn, '/data')
        assert_eq(df, out)

    # test multiple files single node
    with tmpdir() as dn:
        fn = os.path.join(dn, 'data_*.h5')
        a.to_hdf(fn, '/data', get=get)
        out = dd.read_hdf(fn, '/data')
        assert_eq(df, out)

    # test single file multiple nodes
    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data*', get=get)
        out = dd.read_hdf(fn, '/data*')
        assert_eq(df, out)


def test_to_hdf_kwargs():
    pytest.importorskip('tables')
    df = pd.DataFrame({'A': ['a', 'aaaa']})
    ddf = dd.from_pandas(df, npartitions=2)
    with tmpfile('h5') as fn:
        ddf.to_hdf(fn, 'foo4', format='table', min_itemsize=4)
        df2 = pd.read_hdf(fn, 'foo4')
        tm.assert_frame_equal(df, df2)

    # test shorthand 't' for table
    with tmpfile('h5') as fn:
        ddf.to_hdf(fn, 'foo4', format='t', min_itemsize=4)
        df2 = pd.read_hdf(fn, 'foo4')
        tm.assert_frame_equal(df, df2)


@pytest.mark.skipif(sys.version_info[:2] == (3, 3),
                    reason="Python3.3 uses pytest2.7.2, w/o warns method")
def test_to_fmt_warns():
    pytest.importorskip('tables')
    df16 = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
                               'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                         'y': [1, 2, 3, 4, 5, 6, 7, 8, 9,
                               10, 11, 12, 13, 14, 15, 16]},
                        index=[1., 2., 3., 4., 5., 6., 7., 8., 9.,
                               10., 11., 12., 13., 14., 15., 16.])
    a = dd.from_pandas(df16, 16)

    # testing warning when breaking order
    with tmpfile('h5') as fn:
        with pytest.warns(None):
            a.to_hdf(fn, '/data*', name_function=str)

    # testing warning when breaking order
    with tmpdir() as dn:
        with pytest.warns(None):
            fn = os.path.join(dn, "data_*.csv")
            a.to_csv(fn, name_function=str)


@pytest.mark.parametrize('data, compare', [
    (pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                   'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.]),
     tm.assert_frame_equal),
    (pd.Series([1, 2, 3, 4], name='a'),
     tm.assert_series_equal),
])
def test_read_hdf(data, compare):
    pytest.importorskip('tables')
    with tmpfile('h5') as fn:
        data.to_hdf(fn, '/data')
        try:
            dd.read_hdf(fn, 'data', chunksize=2, mode='r')
            assert False
        except TypeError as e:
            assert "format='table'" in str(e)

    with tmpfile('h5') as fn:
        data.to_hdf(fn, '/data', format='table')
        a = dd.read_hdf(fn, '/data', chunksize=2, mode='r')
        assert a.npartitions == 2

        compare(a.compute(), data)

        compare(dd.read_hdf(fn, '/data', chunksize=2, start=1, stop=3,
                            mode='r').compute(),
                pd.read_hdf(fn, '/data', start=1, stop=3))

        assert (sorted(dd.read_hdf(fn, '/data', mode='r').dask) ==
                sorted(dd.read_hdf(fn, '/data', mode='r').dask))

    with tmpfile('h5') as fn:
        sorted_data = data.sort_index()
        sorted_data.to_hdf(fn, '/data', format='table')
        a = dd.read_hdf(fn, '/data', chunksize=2, sorted_index=True, mode='r')
        assert a.npartitions == 2

        compare(a.compute(), sorted_data)


def test_read_hdf_multiply_open():
    """Test that we can read from a file that's already opened elsewhere in
    read-only mode."""
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    with tmpfile('h5') as fn:
        df.to_hdf(fn, '/data', format='table')
        with pd.HDFStore(fn, mode='r'):
            dd.read_hdf(fn, '/data', chunksize=2, mode='r')


def test_read_hdf_multiple():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h',
                             'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p'],
                       'y': [1, 2, 3, 4, 5, 6, 7, 8, 9,
                             10, 11, 12, 13, 14, 15, 16]},
                      index=[1., 2., 3., 4., 5., 6., 7., 8., 9.,
                             10., 11., 12., 13., 14., 15., 16.])
    a = dd.from_pandas(df, 16)

    with tmpfile('h5') as fn:
        a.to_hdf(fn, '/data*')
        r = dd.read_hdf(fn, '/data*', sorted_index=True)
        assert a.npartitions == r.npartitions
        assert a.divisions == r.divisions
        assert_eq(a, r)


def test_read_hdf_start_stop_values():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    with tmpfile('h5') as fn:
        df.to_hdf(fn, '/data', format='table')

        with pytest.raises(ValueError) as e:
            dd.read_hdf(fn, '/data', stop=10)
        assert 'number of rows' in str(e)

        with pytest.raises(ValueError) as e:
            dd.read_hdf(fn, '/data', start=10)
        assert 'is above or equal to' in str(e)

        with pytest.raises(ValueError) as e:
            dd.read_hdf(fn, '/data', chunksize=-1)
        assert 'positive integer' in str(e)


def test_hdf_globbing():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])

    with tmpdir() as tdir:
        df.to_hdf(os.path.join(tdir, 'one.h5'), '/foo/data', format='table')
        df.to_hdf(os.path.join(tdir, 'two.h5'), '/bar/data', format='table')
        df.to_hdf(os.path.join(tdir, 'two.h5'), '/foo/data', format='table')

        with dask.set_options(get=dask.get):
            res = dd.read_hdf(os.path.join(tdir, 'one.h5'), '/*/data',
                              chunksize=2)
            assert res.npartitions == 2
            tm.assert_frame_equal(res.compute(), df)

            res = dd.read_hdf(os.path.join(tdir, 'one.h5'), '/*/data',
                              chunksize=2, start=1, stop=3)
            expected = pd.read_hdf(os.path.join(tdir, 'one.h5'), '/foo/data',
                                   start=1, stop=3)
            tm.assert_frame_equal(res.compute(), expected)

            res = dd.read_hdf(os.path.join(tdir, 'two.h5'), '/*/data', chunksize=2)
            assert res.npartitions == 2 + 2
            tm.assert_frame_equal(res.compute(), pd.concat([df] * 2))

            res = dd.read_hdf(os.path.join(tdir, '*.h5'), '/foo/data', chunksize=2)
            assert res.npartitions == 2 + 2
            tm.assert_frame_equal(res.compute(), pd.concat([df] * 2))

            res = dd.read_hdf(os.path.join(tdir, '*.h5'), '/*/data', chunksize=2)
            assert res.npartitions == 2 + 2 + 2
            tm.assert_frame_equal(res.compute(), pd.concat([df] * 3))


def test_hdf_file_list():
    pytest.importorskip('tables')
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])

    with tmpdir() as tdir:
        df.iloc[:2].to_hdf(os.path.join(tdir, 'one.h5'), 'dataframe', format='table')
        df.iloc[2:].to_hdf(os.path.join(tdir, 'two.h5'), 'dataframe', format='table')

        with dask.set_options(get=dask.get):
            input_files = [os.path.join(tdir, 'one.h5'), os.path.join(tdir, 'two.h5')]
            res = dd.read_hdf(input_files, 'dataframe')
            tm.assert_frame_equal(res.compute(), df)


def test_read_hdf_doesnt_segfault():
    pytest.importorskip('tables')
    with tmpfile('h5') as fn:
        N = 40
        df = pd.DataFrame(np.random.randn(N, 3))
        with pd.HDFStore(fn, mode='w') as store:
            store.append('/x', df)

        ddf = dd.read_hdf(fn, '/x', chunksize=2)
        assert len(ddf) == N


def test_hdf_filenames():
    df = pd.DataFrame({'x': ['a', 'b', 'c', 'd'],
                       'y': [1, 2, 3, 4]}, index=[1., 2., 3., 4.])
    ddf = dd.from_pandas(df, npartitions=2)
    assert ddf.to_hdf("foo*.hdf5", "key") == ["foo0.hdf5", "foo1.hdf5"]
    os.remove("foo0.hdf5")
    os.remove("foo1.hdf5")
