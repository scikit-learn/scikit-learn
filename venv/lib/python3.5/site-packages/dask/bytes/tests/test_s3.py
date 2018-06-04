from __future__ import print_function, division, absolute_import

import sys
from contextlib import contextmanager

import pytest

s3fs = pytest.importorskip('s3fs')
boto3 = pytest.importorskip('boto3')
moto = pytest.importorskip('moto')
httpretty = pytest.importorskip('httpretty')

from toolz import concat, valmap, partial

from dask import compute, get
from dask.bytes.s3 import DaskS3FileSystem
from dask.bytes.core import read_bytes, open_files, get_pyarrow_filesystem
from dask.bytes.compression import compress, files as compress_files, seekable_files


compute = partial(compute, get=get)


test_bucket_name = 'test'
files = {'test/accounts.1.json':  (b'{"amount": 100, "name": "Alice"}\n'
                                   b'{"amount": 200, "name": "Bob"}\n'
                                   b'{"amount": 300, "name": "Charlie"}\n'
                                   b'{"amount": 400, "name": "Dennis"}\n'),
         'test/accounts.2.json':  (b'{"amount": 500, "name": "Alice"}\n'
                                   b'{"amount": 600, "name": "Bob"}\n'
                                   b'{"amount": 700, "name": "Charlie"}\n'
                                   b'{"amount": 800, "name": "Dennis"}\n')}


@pytest.yield_fixture
def s3():
    # writable local S3 system
    import moto
    try:
        with moto.mock_s3():
            client = boto3.client('s3')
            client.create_bucket(Bucket=test_bucket_name, ACL='public-read-write')
            for f, data in files.items():
                client.put_object(Bucket=test_bucket_name, Key=f, Body=data)
            yield s3fs.S3FileSystem(anon=True)
    finally:
        httpretty.HTTPretty.disable()
        httpretty.HTTPretty.reset()


@contextmanager
def s3_context(bucket, files):
    m = moto.mock_s3()
    m.start()
    client = boto3.client('s3')
    client.create_bucket(Bucket=bucket, ACL='public-read-write')
    for f, data in files.items():
        client.put_object(Bucket=bucket, Key=f, Body=data)

    yield DaskS3FileSystem(anon=True)

    for f, data in files.items():
        try:
            client.delete_object(Bucket=bucket, Key=f, Body=data)
        except Exception:
            pass
        finally:
            m.stop()
            httpretty = pytest.importorskip('httpretty')
            httpretty.HTTPretty.disable()
            httpretty.HTTPretty.reset()


def test_get_s3():
    s3 = DaskS3FileSystem(key='key', secret='secret')
    assert s3.key == 'key'
    assert s3.secret == 'secret'

    s3 = DaskS3FileSystem(username='key', password='secret')
    assert s3.key == 'key'
    assert s3.secret == 'secret'

    with pytest.raises(KeyError):
        DaskS3FileSystem(key='key', username='key')
    with pytest.raises(KeyError):
        DaskS3FileSystem(secret='key', password='key')


def test_open_files_write(s3):
    paths = ['s3://' + test_bucket_name + '/more/' + f for f in files]
    fils = open_files(paths, mode='wb')
    for fil, data in zip(fils, files.values()):
        with fil as f:
            f.write(data)
    sample, values = read_bytes('s3://' + test_bucket_name + '/more/test/accounts.*')
    results = compute(*concat(values))
    assert set(list(files.values())) == set(results)


def test_read_bytes(s3):
    sample, values = read_bytes('s3://' + test_bucket_name + '/test/accounts.*')
    assert isinstance(sample, bytes)
    assert sample[:5] == files[sorted(files)[0]][:5]
    assert sample.endswith(b'\n')

    assert isinstance(values, (list, tuple))
    assert isinstance(values[0], (list, tuple))
    assert hasattr(values[0][0], 'dask')

    assert sum(map(len, values)) >= len(files)
    results = compute(*concat(values))
    assert set(results) == set(files.values())


def test_read_bytes_sample_delimiter(s3):
    sample, values = read_bytes('s3://' + test_bucket_name + '/test/accounts.*',
                                sample=80, delimiter=b'\n')
    assert sample.endswith(b'\n')
    sample, values = read_bytes('s3://' + test_bucket_name + '/test/accounts.1.json',
                                sample=80, delimiter=b'\n')
    assert sample.endswith(b'\n')
    sample, values = read_bytes('s3://' + test_bucket_name + '/test/accounts.1.json',
                                sample=2, delimiter=b'\n')
    assert sample.endswith(b'\n')


def test_read_bytes_non_existing_glob(s3):
    with pytest.raises(IOError):
        read_bytes('s3://' + test_bucket_name + '/non-existing/*')


def test_read_bytes_blocksize_none(s3):
    _, values = read_bytes('s3://' + test_bucket_name + '/test/accounts.*',
                           blocksize=None)
    assert sum(map(len, values)) == len(files)


@pytest.mark.slow
@pytest.mark.network
def test_read_bytes_blocksize_on_large_data():
    _, L = read_bytes('s3://dask-data/nyc-taxi/2015/yellow_tripdata_2015-01.csv',
                      blocksize=None, anon=True)
    assert len(L) == 1

    _, L = read_bytes('s3://dask-data/nyc-taxi/2014/*.csv', blocksize=None, anon=True)
    assert len(L) == 12


@pytest.mark.parametrize('blocksize', [5, 15, 45, 1500])
def test_read_bytes_block(s3, blocksize):
    _, vals = read_bytes('s3://' + test_bucket_name + '/test/account*',
                         blocksize=blocksize)
    assert (list(map(len, vals)) ==
            [(len(v) // blocksize + 1) for v in files.values()])

    results = compute(*concat(vals))
    assert (sum(len(r) for r in results) ==
            sum(len(v) for v in files.values()))

    ourlines = b"".join(results).split(b'\n')
    testlines = b"".join(files.values()).split(b'\n')
    assert set(ourlines) == set(testlines)


@pytest.mark.parametrize('blocksize', [5, 15, 45, 1500])
def test_read_bytes_delimited(s3, blocksize):
    _, values = read_bytes('s3://' + test_bucket_name + '/test/accounts*',
                           blocksize=blocksize, delimiter=b'\n')
    _, values2 = read_bytes('s3://' + test_bucket_name + '/test/accounts*',
                            blocksize=blocksize, delimiter=b'foo')
    assert ([a.key for a in concat(values)] !=
            [b.key for b in concat(values2)])

    results = compute(*concat(values))
    res = [r for r in results if r]
    assert all(r.endswith(b'\n') for r in res)
    ourlines = b''.join(res).split(b'\n')
    testlines = b"".join(files[k] for k in sorted(files)).split(b'\n')
    assert ourlines == testlines

    # delimiter not at the end
    d = b'}'
    _, values = read_bytes('s3://' + test_bucket_name + '/test/accounts*',
                           blocksize=blocksize, delimiter=d)
    results = compute(*concat(values))
    res = [r for r in results if r]
    # All should end in } except EOF
    assert sum(r.endswith(b'}') for r in res) == len(res) - 2
    ours = b"".join(res)
    test = b"".join(files[v] for v in sorted(files))
    assert ours == test


@pytest.mark.parametrize('fmt,blocksize',
                         [(fmt, None) for fmt in compress_files] +
                         [(fmt, 10) for fmt in seekable_files])
def test_compression(s3, fmt, blocksize):
    with s3_context('compress', valmap(compress[fmt], files)):
        sample, values = read_bytes('s3://compress/test/accounts.*',
                                    compression=fmt, blocksize=blocksize)
        assert sample.startswith(files[sorted(files)[0]][:10])
        assert sample.endswith(b'\n')

        results = compute(*concat(values))
        assert b''.join(results) == b''.join([files[k] for k in sorted(files)])


@pytest.mark.parametrize('mode', ['rt', 'rb'])
def test_open_files(s3, mode):
    myfiles = open_files('s3://' + test_bucket_name + '/test/accounts.*',
                         mode=mode)
    assert len(myfiles) == len(files)
    for lazy_file, path in zip(myfiles, sorted(files)):
        with lazy_file as f:
            data = f.read()
            sol = files[path]
            assert data == sol if mode == 'rb' else sol.decode()


double = lambda x: x * 2


def test_modification_time_read_bytes():
    with s3_context('compress', files):
        _, a = read_bytes('s3://compress/test/accounts.*')
        _, b = read_bytes('s3://compress/test/accounts.*')

        assert [aa._key for aa in concat(a)] == [bb._key for bb in concat(b)]

    with s3_context('compress', valmap(double, files)):
        _, c = read_bytes('s3://compress/test/accounts.*')

    assert [aa._key for aa in concat(a)] != [cc._key for cc in concat(c)]


def test_read_csv_passes_through_options():
    dd = pytest.importorskip('dask.dataframe')
    with s3_context('csv', {'a.csv': b'a,b\n1,2\n3,4'}) as s3:
        df = dd.read_csv('s3://csv/*.csv', storage_options={'s3': s3})
        assert df.a.sum().compute() == 1 + 3


def test_read_text_passes_through_options():
    db = pytest.importorskip('dask.bag')
    with s3_context('csv', {'a.csv': b'a,b\n1,2\n3,4'}) as s3:
        df = db.read_text('s3://csv/*.csv', storage_options={'s3': s3})
        assert df.count().compute(get=get) == 3


@pytest.mark.parametrize("engine", ['pyarrow', 'fastparquet'])
def test_parquet(s3, engine):
    dd = pytest.importorskip('dask.dataframe')
    pytest.importorskip(engine)
    import pandas as pd
    import numpy as np

    url = 's3://%s/test.parquet' % test_bucket_name

    data = pd.DataFrame({'i32': np.arange(1000, dtype=np.int32),
                         'i64': np.arange(1000, dtype=np.int64),
                         'f': np.arange(1000, dtype=np.float64),
                         'bhello': np.random.choice(
                             [u'hello', u'you', u'people'],
                             size=1000).astype("O")},
                        index=pd.Index(np.arange(1000), name='foo'))
    df = dd.from_pandas(data, chunksize=500)
    df.to_parquet(url, engine=engine)

    files = [f.split('/')[-1] for f in s3.ls(url)]
    assert '_common_metadata' in files
    assert 'part.0.parquet' in files

    df2 = dd.read_parquet(url, index='foo', engine=engine)
    assert len(df2.divisions) > 1

    pd.util.testing.assert_frame_equal(data, df2.compute())


def test_parquet_wstoragepars(s3):
    dd = pytest.importorskip('dask.dataframe')
    pytest.importorskip('fastparquet')

    import pandas as pd
    import numpy as np

    url = 's3://%s/test.parquet' % test_bucket_name

    data = pd.DataFrame({'i32': np.array([0, 5, 2, 5])})
    df = dd.from_pandas(data, chunksize=500)
    df.to_parquet(url, write_index=False)

    dd.read_parquet(url, storage_options={'default_fill_cache': False})
    assert s3.current().default_fill_cache is False
    dd.read_parquet(url, storage_options={'default_fill_cache': True})
    assert s3.current().default_fill_cache is True

    dd.read_parquet(url, storage_options={'default_block_size': 2**20})
    assert s3.current().default_block_size == 2**20
    with s3.current().open(url + '/_metadata') as f:
        assert f.blocksize == 2**20


@pytest.mark.skipif(sys.platform == 'win32',
                    reason="pathlib and moto clash on windows")
def test_pathlib_s3(s3):
    pathlib = pytest.importorskip("pathlib")
    with pytest.raises(ValueError):
        url = pathlib.Path('s3://bucket/test.accounts.*')
        sample, values = read_bytes(url, blocksize=None)


def test_get_pyarrow_fs_s3(s3):
    pa = pytest.importorskip('pyarrow')
    fs = DaskS3FileSystem(anon=True)
    assert isinstance(get_pyarrow_filesystem(fs), pa.filesystem.S3FSWrapper)
