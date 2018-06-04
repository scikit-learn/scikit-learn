import os
import pytest
import requests
import subprocess
import sys
import time

from dask.bytes.core import open_files
from dask.compatibility import PY2
from dask.utils import tmpdir

files = ['a', 'b']


@pytest.fixture(scope='module')
def dir_server():
    with tmpdir() as d:
        for fn in files:
            with open(os.path.join(d, fn), 'wb') as f:
                f.write(b'a' * 10000)

        if PY2:
            cmd = [sys.executable, '-m', 'SimpleHTTPServer', '8999']
        else:
            cmd = [sys.executable, '-m', 'http.server', '8999']
        p = subprocess.Popen(cmd, cwd=d)
        timeout = 10
        while True:
            try:
                requests.get('http://localhost:8999')
                break
            except requests.exceptions.ConnectionError:
                time.sleep(0.1)
                timeout -= 0.1
                if timeout < 0:
                    raise RuntimeError('Server did not appear')
        yield d
        p.terminate()


def test_simple(dir_server):
    root = 'http://localhost:8999/'
    fn = files[0]
    f = open_files(root + fn)[0]
    with f as f:
        data = f.read()
    assert data == open(os.path.join(dir_server, fn), 'rb').read()


@pytest.mark.parametrize('block_size', [None, 99999])
def test_ops(dir_server, block_size):
    root = 'http://localhost:8999/'
    fn = files[0]
    f = open_files(root + fn)[0]
    data = open(os.path.join(dir_server, fn), 'rb').read()
    with f as f:
        # these pass because the default
        assert f.read(10) == data[:10]
        f.seek(0)
        assert f.read(10) == data[:10]
        assert f.read(10) == data[10:20]
        f.seek(-10, 2)
        assert f.read() == data[-10:]


def test_ops_blocksize(dir_server):
    root = 'http://localhost:8999/'
    fn = files[0]
    f = open_files(root + fn, block_size=2)[0]
    data = open(os.path.join(dir_server, fn), 'rb').read()
    with f as f:
        # it's OK to read the whole file
        assert f.read() == data

    # note that if we reuse f from above, because it is tokenized, we get
    # the same open file - where is this cached?
    fn = files[1]
    f = open_files(root + fn, block_size=2)[0]
    with f as f:
        # fails becasue we want only 12 bytes
        with pytest.raises(ValueError):
            assert f.read(10) == data[:10]


def test_errors(dir_server):
    f = open_files('http://localhost:8999/doesnotexist')[0]
    with pytest.raises(requests.exceptions.RequestException):
        with f:
            pass
    f = open_files('http://nohost/')[0]
    with pytest.raises(requests.exceptions.RequestException):
        with f:
            pass
    root = 'http://localhost:8999/'
    fn = files[0]
    f = open_files(root + fn, mode='wb')[0]
    with pytest.raises(NotImplementedError):
        with f:
            pass
    f = open_files(root + fn)[0]
    with f as f:
        with pytest.raises(ValueError):
            f.seek(-1)


def test_files(dir_server):
    root = 'http://localhost:8999/'
    fs = open_files([root + f for f in files])
    for f, f2 in zip(fs, files):
        with f as f:
            assert f.read() == open(os.path.join(dir_server, f2), 'rb').read()


@pytest.mark.network
def test_parquet():
    dd = pytest.importorskip('dask.dataframe')
    pytest.importorskip('fastparquet')  # no pyarrow compatability FS yet
    df = dd.read_parquet([
        'https://github.com/Parquet/parquet-compatibility/raw/'
        'master/parquet-testdata/impala/1.1.1-NONE/'
        'nation.impala.parquet']).compute()
    assert df.n_nationkey.tolist() == list(range(25))
    assert df.columns.tolist() == ['n_nationkey', 'n_name', 'n_regionkey',
                                   'n_comment']


@pytest.mark.network
def test_bag():
    # This test pulls from different hosts
    db = pytest.importorskip('dask.bag')
    urls = ['https://raw.githubusercontent.com/weierophinney/pastebin/'
            'master/public/js-src/dojox/data/tests/stores/patterns.csv',
            'https://en.wikipedia.org']
    b = db.read_text(urls)
    assert b.npartitions == 2
    b.compute()
