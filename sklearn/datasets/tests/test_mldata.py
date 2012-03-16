"""Test functionality of mldata fetching utilities."""

from sklearn import datasets
from sklearn.datasets import mldata_filename, fetch_mldata
from sklearn.utils.testing import assert_in
from sklearn.utils.testing import UrlopenMock
from nose.tools import assert_equal, assert_raises
from nose import with_setup
from numpy.testing import assert_array_equal
import os
import shutil
import tempfile
import scipy as sp

tmpdir = None


def setup_tmpdata():
    # create temporary dir
    global tmpdir
    tmpdir = tempfile.mkdtemp()
    os.makedirs(os.path.join(tmpdir, 'mldata'))


def teardown_tmpdata():
    # remove temporary dir
    if tmpdir is not None:
        shutil.rmtree(tmpdir)


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
    mock = UrlopenMock(datasets.mldata,
                      {'mock':
                       {'label': sp.ones((150,)),
                        'data': sp.ones((150, 4))}}).install()
    try:
        dset = fetch_mldata('mock', data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data'])

        assert_equal(dset.target.shape, (150,))
        assert_equal(dset.data.shape, (150, 4))

        assert_raises(UrlopenMock.HTTPError,
                      fetch_mldata, 'not_existing_name')
    finally:
        mock.restaure()


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_fetch_one_column():
    dataname = 'onecol'
    x = sp.arange(6).reshape(2, 3)

    # create fake data set in cache
    mock = UrlopenMock(datasets.mldata, {dataname: {'x': x}}).install()
    try:
        dset = fetch_mldata(dataname, data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'data'], out_=['target'])

        assert_equal(dset.data.shape, (2, 3))
        assert_array_equal(dset.data, x)

        # transposing the data array
        dset = fetch_mldata(dataname, transpose_data=False, data_home=tmpdir)
        assert_equal(dset.data.shape, (3, 2))
    finally:
        mock.restaure()


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_fetch_multiple_column():
    # create fake data set in cache
    x = sp.arange(6).reshape(2, 3)
    y = sp.array([1, -1])
    z = sp.arange(12).reshape(4, 3)

    # by default
    dataname = 'threecol-default'
    mock = UrlopenMock(datasets.mldata,
                      {dataname:
                       ({'label': y, 'data': x, 'z': z},
                        ['z', 'data', 'label'])}).install()
    try:
        dset = fetch_mldata(dataname, data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data', 'z'],
                  out_=['x', 'y'])

        assert_array_equal(dset.data, x)
        assert_array_equal(dset.target, y)
        assert_array_equal(dset.z, z.T)
    finally:
        mock.restaure()

    # by order
    dataname = 'threecol-order'
    mock = UrlopenMock(datasets.mldata,
                      {dataname:
                       ({'y': y, 'x': x, 'z': z},
                        ['y', 'x', 'z'])}).install()
    try:
        dset = fetch_mldata(dataname, data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data', 'z'],
                  out_=['x', 'y'])
        assert_array_equal(dset.data, x)
        assert_array_equal(dset.target, y)
        assert_array_equal(dset.z, z.T)
    finally:
        mock.restaure()

    # by number
    dataname = 'threecol-number'
    mock = UrlopenMock(datasets.mldata,
                      {dataname:
                       ({'y': y, 'x': x, 'z': z},
                        ['z', 'x', 'y'])}).install()
    try:
        dset = fetch_mldata(dataname, target_name=2, data_name=0,
                            data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data', 'x'],
                  out_=['z', 'y'])
        assert_array_equal(dset.data, z)
        assert_array_equal(dset.target, y)

        # by name
        dset = fetch_mldata(dataname, target_name='y', data_name='z',
                            data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data', 'x'],
                  out_=['z', 'y'])
        assert_array_equal(dset.data, z)
        assert_array_equal(dset.target, y)
    finally:
        mock.restaure()
