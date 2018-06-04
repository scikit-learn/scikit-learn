import os
import pytest
import shutil
import tempfile

from dask.dataframe import read_orc
from dask.dataframe.utils import assert_eq
import dask.dataframe as dd

pytest.importorskip('pyarrow.orc')
url = ('https://www.googleapis.com/download/storage/v1/b/anaconda-public-data/o'
       '/orc%2FTestOrcFile.testDate1900.orc?generation=1522611448751555&alt='
       'media')
columns = ['time', 'date']


@pytest.mark.network
def test_orc_with_backend():
    d = read_orc(url)
    assert set(d.columns) == {'time', 'date'}  # order is not guranteed
    assert len(d) == 70000


@pytest.fixture(scope='module')
def orc_files():
    requests = pytest.importorskip('requests')
    data = requests.get(url).content
    d = tempfile.mkdtemp()
    files = [os.path.join(d, fn) for fn in ['test1.orc', 'test2.orc']]
    for fn in files:
        with open(fn, 'wb') as f:
            f.write(data)
    try:
        yield files
    finally:
        shutil.rmtree(d, ignore_errors=True)


def test_orc_single(orc_files):
    fn = orc_files[0]
    d = read_orc(fn)
    assert len(d) == 70000
    assert d.npartitions == 8
    d2 = read_orc(fn, columns=['time', 'date'])
    assert_eq(d[columns], d2[columns])
    with pytest.raises(ValueError) as e:
        read_orc(fn, columns=['time', 'nonexist'])
    assert 'nonexist' in str(e)


def test_orc_multiple(orc_files):
    d = read_orc(orc_files[0])
    d2 = read_orc(orc_files)
    assert_eq(d2[columns], dd.concat([d, d])[columns], check_index=False)
    d2 = read_orc(os.path.dirname(orc_files[0]) + '/*.orc')
    assert_eq(d2[columns], dd.concat([d, d])[columns], check_index=False)
