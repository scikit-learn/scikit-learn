from __future__ import print_function, division, absolute_import

import gzip
import os
from time import sleep
import sys

import pytest
from toolz import concat, valmap, partial

from dask import compute, get
from dask.compatibility import FileNotFoundError, unicode
from dask.utils import filetexts
from dask.bytes import compression
from dask.bytes.local import LocalFileSystem
from dask.bytes.core import (read_bytes, open_files, FileSystem,
                             get_pyarrow_filesystem, logical_size,
                             get_fs_token_paths)

compute = partial(compute, get=get)

files = {'.test.accounts.1.json': (b'{"amount": 100, "name": "Alice"}\n'
                                   b'{"amount": 200, "name": "Bob"}\n'
                                   b'{"amount": 300, "name": "Charlie"}\n'
                                   b'{"amount": 400, "name": "Dennis"}\n'),
         '.test.accounts.2.json': (b'{"amount": 500, "name": "Alice"}\n'
                                   b'{"amount": 600, "name": "Bob"}\n'
                                   b'{"amount": 700, "name": "Charlie"}\n'
                                   b'{"amount": 800, "name": "Dennis"}\n')}


try:
    # used only in test_with_urls - may be more generally useful
    import pathlib

    def to_uri(path):
        return pathlib.Path(os.path.abspath(path)).as_uri()

except (ImportError, NameError):
    import urlparse, urllib

    def to_uri(path):
        return urlparse.urljoin(
            'file:', urllib.pathname2url(os.path.abspath(path)))


def test_urlpath_inference_strips_protocol(tmpdir):
    tmpdir = str(tmpdir)
    paths = [os.path.join(tmpdir, 'test.%02d.csv' % i) for i in range(20)]

    for path in paths:
        with open(path, 'wb') as f:
            f.write(b'1,2,3\n' * 10)

    # globstring
    protocol = 'file:///' if sys.platform == 'win32' else 'file://'
    urlpath = protocol + os.path.join(tmpdir, 'test.*.csv')
    _, _, paths2 = get_fs_token_paths(urlpath)
    assert paths2 == paths

    # list of paths
    _, _, paths2 = get_fs_token_paths([protocol + p for p in paths])
    assert paths2 == paths


def test_urlpath_inference_errors():
    # Empty list
    with pytest.raises(ValueError) as err:
        get_fs_token_paths([])
    assert 'empty' in str(err)

    # Protocols differ
    with pytest.raises(ValueError) as err:
        get_fs_token_paths(['s3://test/path.csv', '/other/path.csv'])
    assert 'same protocol and options' in str(err)

    # Options differ
    with pytest.raises(ValueError) as err:
        get_fs_token_paths(['hdfs://myuser@node.com/test/path.csv',
                            'hdfs://otheruser@node.com/other/path.csv'])
    assert 'same protocol and options' in str(err)

    # Unknown type
    with pytest.raises(TypeError):
        get_fs_token_paths({'sets/are.csv', 'unordered/so/they.csv',
                            'should/not/be.csv' 'allowed.csv'})


def test_read_bytes():
    with filetexts(files, mode='b'):
        sample, values = read_bytes('.test.accounts.*')
        assert isinstance(sample, bytes)
        assert sample[:5] == files[sorted(files)[0]][:5]
        assert sample.endswith(b'\n')

        assert isinstance(values, (list, tuple))
        assert isinstance(values[0], (list, tuple))
        assert hasattr(values[0][0], 'dask')

        assert sum(map(len, values)) >= len(files)
        results = compute(*concat(values))
        assert set(results) == set(files.values())


def test_read_bytes_sample_delimiter():
    with filetexts(files, mode='b'):
        sample, values = read_bytes('.test.accounts.*',
                                    sample=80, delimiter=b'\n')
        assert sample.endswith(b'\n')
        sample, values = read_bytes('.test.accounts.1.json',
                                    sample=80, delimiter=b'\n')
        assert sample.endswith(b'\n')
        sample, values = read_bytes('.test.accounts.1.json',
                                    sample=2, delimiter=b'\n')
        assert sample.endswith(b'\n')


def test_read_bytes_blocksize_none():
    with filetexts(files, mode='b'):
        sample, values = read_bytes('.test.accounts.*', blocksize=None)
        assert sum(map(len, values)) == len(files)


def test_read_bytes_blocksize_float():
    with filetexts(files, mode='b'):
        sample, vals = read_bytes('.test.account*', blocksize=5.0)
        results = compute(*concat(vals))
        ourlines = b"".join(results).split(b'\n')
        testlines = b"".join(files.values()).split(b'\n')
        assert set(ourlines) == set(testlines)

        with pytest.raises(TypeError):
            read_bytes('.test.account*', blocksize=5.5)


def test_with_urls():
    with filetexts(files, mode='b'):
        # OS-independent file:// URI with glob *
        url = to_uri('.test.accounts.') + '*'
        sample, values = read_bytes(url, blocksize=None)
        assert sum(map(len, values)) == len(files)


@pytest.mark.skipif(sys.platform == 'win32',
                    reason="pathlib and moto clash on windows")
def test_with_paths():
    pathlib = pytest.importorskip('pathlib')
    with filetexts(files, mode='b'):
        url = pathlib.Path('./.test.accounts.*')
        sample, values = read_bytes(url, blocksize=None)
        assert sum(map(len, values)) == len(files)
    with pytest.raises(OSError):
        # relative path doesn't work
        url = pathlib.Path('file://.test.accounts.*')
        read_bytes(url, blocksize=None)


def test_read_bytes_block():
    with filetexts(files, mode='b'):
        for bs in [5, 15, 45, 1500]:
            sample, vals = read_bytes('.test.account*', blocksize=bs)
            assert (list(map(len, vals)) ==
                    [(len(v) // bs + 1) for v in files.values()])

            results = compute(*concat(vals))
            assert (sum(len(r) for r in results) ==
                    sum(len(v) for v in files.values()))

            ourlines = b"".join(results).split(b'\n')
            testlines = b"".join(files.values()).split(b'\n')
            assert set(ourlines) == set(testlines)


def test_read_bytes_delimited():
    with filetexts(files, mode='b'):
        for bs in [5, 15, 45, 1500]:
            _, values = read_bytes('.test.accounts*',
                                   blocksize=bs, delimiter=b'\n')
            _, values2 = read_bytes('.test.accounts*',
                                    blocksize=bs, delimiter=b'foo')
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
            _, values = read_bytes('.test.accounts*', blocksize=bs, delimiter=d)
            results = compute(*concat(values))
            res = [r for r in results if r]
            # All should end in } except EOF
            assert sum(r.endswith(b'}') for r in res) == len(res) - 2
            ours = b"".join(res)
            test = b"".join(files[v] for v in sorted(files))
            assert ours == test


fmt_bs = ([(fmt, None) for fmt in compression.files] +
          [(fmt, 10) for fmt in compression.seekable_files])


@pytest.mark.parametrize('fmt,blocksize', fmt_bs)
def test_compression(fmt, blocksize):
    compress = compression.compress[fmt]
    files2 = valmap(compress, files)
    with filetexts(files2, mode='b'):
        sample, values = read_bytes('.test.accounts.*.json',
                                    blocksize=blocksize, delimiter=b'\n',
                                    compression=fmt)
        assert sample[:5] == files[sorted(files)[0]][:5]
        assert sample.endswith(b'\n')

        results = compute(*concat(values))
        assert (b''.join(results) ==
                b''.join([files[k] for k in sorted(files)]))


def test_open_files():
    with filetexts(files, mode='b'):
        myfiles = open_files('.test.accounts.*')
        assert len(myfiles) == len(files)
        for lazy_file, data_file in zip(myfiles, sorted(files)):
            with lazy_file as f:
                x = f.read()
                assert x == files[data_file]


@pytest.mark.parametrize('encoding', ['utf-8', 'ascii'])
def test_open_files_text_mode(encoding):
    with filetexts(files, mode='b'):
        myfiles = open_files('.test.accounts.*', mode='rt', encoding=encoding)
        assert len(myfiles) == len(files)
        data = []
        for file in myfiles:
            with file as f:
                data.append(f.read())
        assert list(data) == [files[k].decode(encoding)
                              for k in sorted(files)]


@pytest.mark.parametrize('mode', ['rt', 'rb'])
@pytest.mark.parametrize('fmt', list(compression.files))
def test_open_files_compression(mode, fmt):
    files2 = valmap(compression.compress[fmt], files)
    with filetexts(files2, mode='b'):
        myfiles = open_files('.test.accounts.*', mode=mode, compression=fmt)
        data = []
        for file in myfiles:
            with file as f:
                data.append(f.read())
        sol = [files[k] for k in sorted(files)]
        if mode == 'rt':
            sol = [b.decode() for b in sol]
        assert list(data) == sol


@pytest.mark.parametrize('fmt', list(compression.seekable_files))
def test_getsize(fmt):
    compress = compression.compress[fmt]
    with filetexts({'.tmp.getsize': compress(b'1234567890')}, mode='b'):
        fs = LocalFileSystem()
        assert logical_size(fs, '.tmp.getsize', fmt) == 10


def test_bad_compression():
    with filetexts(files, mode='b'):
        for func in [read_bytes, open_files]:
            with pytest.raises(ValueError):
                sample, values = func('.test.accounts.*',
                                      compression='not-found')


def test_not_found():
    fn = 'not-a-file'
    with pytest.raises((FileNotFoundError, OSError)) as e:
        read_bytes(fn)
    assert fn in str(e)


@pytest.mark.slow
def test_names():
    with filetexts(files, mode='b'):
        _, a = read_bytes('.test.accounts.*')
        _, b = read_bytes('.test.accounts.*')
        a = list(concat(a))
        b = list(concat(b))

        assert [aa._key for aa in a] == [bb._key for bb in b]

        sleep(1)
        for fn in files:
            with open(fn, 'ab') as f:
                f.write(b'x')

        _, c = read_bytes('.test.accounts.*')
        c = list(concat(c))
        assert [aa._key for aa in a] != [cc._key for cc in c]


@pytest.mark.parametrize('compression_opener',
                         [(None, open), ('gzip', gzip.open)])
def test_open_files_write(tmpdir, compression_opener):
    compression, opener = compression_opener
    tmpdir = str(tmpdir)
    files = open_files(tmpdir, num=2, mode='wb', compression=compression)
    assert len(files) == 2
    assert {f.mode for f in files} == {'wb'}
    for fil in files:
        with fil as f:
            f.write(b'000')
    files = sorted(os.listdir(tmpdir))
    assert files == ['0.part', '1.part']

    with opener(os.path.join(tmpdir, files[0]), 'rb') as f:
        d = f.read()
    assert d == b'000'


def test_pickability_of_lazy_files(tmpdir):
    tmpdir = str(tmpdir)
    cloudpickle = pytest.importorskip('cloudpickle')

    with filetexts(files, mode='b'):
        myfiles = open_files('.test.accounts.*')
        myfiles2 = cloudpickle.loads(cloudpickle.dumps(myfiles))

        for f, f2 in zip(myfiles, myfiles2):
            assert f.path == f2.path
            assert type(f.fs) == type(f2.fs)
            with f as f_open, f2 as f2_open:
                assert f_open.read() == f2_open.read()


def test_py2_local_bytes(tmpdir):
    fn = str(tmpdir / 'myfile.txt.gz')
    with gzip.open(fn, mode='wb') as f:
        f.write(b'hello\nworld')

    files = open_files(fn, compression='gzip', mode='rt')

    with files[0] as f:
        assert all(isinstance(line, unicode) for line in f)


def test_abs_paths(tmpdir):
    tmpdir = str(tmpdir)
    here = os.getcwd()
    os.chdir(tmpdir)
    with open('tmp', 'w') as f:
        f.write('hi')
    out = LocalFileSystem().glob('*')
    assert len(out) == 1
    assert os.sep in out[0]
    assert tmpdir in out[0] and 'tmp' in out[0]

    fs = LocalFileSystem()
    os.chdir(here)
    with fs.open('tmp', 'r') as f:
        res = f.read()
    assert res == 'hi'


class UnknownFileSystem(FileSystem):
    pass


def test_get_pyarrow_filesystem():
    pa = pytest.importorskip('pyarrow')

    fs = LocalFileSystem()
    assert isinstance(get_pyarrow_filesystem(fs), pa.filesystem.LocalFileSystem)

    with pytest.raises(NotImplementedError):
        get_pyarrow_filesystem(UnknownFileSystem())
