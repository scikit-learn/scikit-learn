"""Test functionality of mldata fetching utilities."""

from scikits.learn import datasets
from scikits.learn.datasets import mldata_filename, fetch_mldata
from nose.tools import assert_equal, assert_raises
from nose import SkipTest, with_setup
from numpy.testing import assert_array_equal
import os
import shutil
import tempfile
import urllib2
import scipy as sp
from scipy.io import savemat


def setup_tmpdata():
    # create temporary dir
    global tmpdir
    tmpdir = tempfile.mkdtemp()
    os.makedirs(os.path.join(tmpdir, 'mldata'))


def teardown_tmpdata():
    # remove temporary dir
    shutil.rmtree(tmpdir)


def fake_mldata_cache(columns_dict, dataname, cachedir, ordering=None):
    """Create a fake mldata data set in the cache_path.
    
    Parameters
    ----------
    columns_dict: contains data as
                  columns_dict[column_name] = array of data
    dataname: name of data set
    cachedir: cache path
    ordering: list of column_names, determines the ordering in the data set
    
    Note: this function transposes all arrays, while fetch_mldata only
    transposes 'data', keep that into account in the tests.
    """
    dset = dict(columns_dict)

    # transpose all variables
    for name in dset:
        dset[name] = dset[name].T

    if ordering is None:
        ordering = dset.keys()
    # NOTE: setting up this array is tricky, because of the way Matlab
    # re-packages 1D arrays
    dset['mldata_descr_ordering'] = sp.empty((1, len(ordering)),
                                             dtype='object')
    for i, name in enumerate(ordering):
        dset['mldata_descr_ordering'][0, i] = name

    savemat(os.path.join(cachedir, 'mldata', dataname), dset)


def test_mldata_filename():
    cases = [('datasets-UCI iris', 'datasets-uci-iris'),
             ('news20.binary', 'news20binary'),
             ('book-crossing-ratings-1.0', 'book-crossing-ratings-10'),
             ('Nile Water Level', 'nile-water-level'),
             ('MNIST (original)', 'mnist-original')]
    for name, desired in cases:
        assert_equal(mldata_filename(name), desired)


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_download():
    """Test that fetch_mldata is able to download and cache a data set."""

    class mock_urllib2(object):
        @staticmethod
        def urlopen(urlname):
            print urlname
            if urlname == 'http://mldata.org/repository/data/download/matlab/mock':
                label = sp.ones((150,))
                data = sp.ones((150, 4))
                fake_mldata_cache({'label': label, 'data': data},
                                  '_mock', tmpdir)
                return open(os.path.join(tmpdir, 'mldata', '_mock.mat'), 'rb')
            else:
                raise mock_urllib2.URLError('%s not found.', urlname)

        class URLError(Exception):
            pass

    _urllib2_ref = datasets.mldata.urllib2
    datasets.mldata.urllib2 = mock_urllib2
    try:
        mock = fetch_mldata('mock', data_home=tmpdir)
        assert 'COL_NAMES' in mock
        assert 'DESCR' in mock
        assert 'target' in mock
        assert 'data' in mock
    
        assert_equal(mock.target.shape, (150,))
        assert_equal(mock.data.shape, (150, 4))

        assert_raises(IOError, fetch_mldata, 'not_existing_name')
    finally:
        datasets.mldata.urllib2 = _urllib2_ref


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_fetch_one_column():
    dataname = 'onecol'
    # create fake data set in cache
    x = sp.arange(6).reshape(2, 3)
    fake_mldata_cache({'x': x}, dataname, tmpdir)

    dset = fetch_mldata(dataname, data_home=tmpdir)
    assert 'COL_NAMES' in dset
    assert 'DESCR' in dset
    assert 'target' not in dset
    assert 'data' in dset

    assert_equal(dset.data.shape, (2, 3))
    assert_array_equal(dset.data, x)

    # transposing the data array
    dset = fetch_mldata(dataname, transpose_data=False, data_home=tmpdir)
    assert_equal(dset.data.shape, (3, 2))


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_fetch_multiple_column():
    dataname = 'threecol'
    # create fake data set in cache
    x = sp.arange(6).reshape(2, 3)
    y = sp.array([1, -1])
    z = sp.arange(12).reshape(4, 3)

    # by default
    fake_mldata_cache({'label': y, 'data': x, 'z': z}, dataname, tmpdir,
                      ordering=['z', 'data', 'label'])

    dset = fetch_mldata(dataname, data_home=tmpdir)
    assert 'target' in dset
    assert 'data' in dset
    assert 'x' not in dset
    assert 'y' not in dset
    assert 'z' in dset

    assert_array_equal(dset.data, x)
    assert_array_equal(dset.target, y)
    assert_array_equal(dset.z, z.T)

    # by order
    fake_mldata_cache({'y': y, 'x': x, 'z': z}, dataname, tmpdir,
                      ordering=['y', 'x', 'z'])

    dset = fetch_mldata(dataname, data_home=tmpdir)
    assert 'target' in dset
    assert 'data' in dset
    assert_array_equal(dset.data, x)
    assert_array_equal(dset.target, y)
    assert_array_equal(dset.z, z.T)

    # by number
    fake_mldata_cache({'y': y, 'x': x, 'z': z}, dataname, tmpdir,
                      ordering=['z', 'x', 'y'])

    dset = fetch_mldata(dataname, data_home=tmpdir,
                        target_name=2, data_name=0)
    assert 'target' in dset
    assert 'data' in dset
    assert_array_equal(dset.data, z)
    assert_array_equal(dset.target, y)

    # by name
    dset = fetch_mldata(dataname, data_home=tmpdir,
                        target_name='y', data_name='z')
    assert 'target' in dset
    assert 'data' in dset
    assert_array_equal(dset.data, z)
    assert_array_equal(dset.target, y)

