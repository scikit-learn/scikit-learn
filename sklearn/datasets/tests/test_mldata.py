"""Test functionality of mldata fetching utilities."""

from sklearn import datasets
from sklearn.datasets import mldata_filename, fetch_mldata
from sklearn.utils.testing import (assert_in, mock_urllib2)
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

    _urllib2_ref = datasets.mldata.urllib2
    datasets.mldata.urllib2 = mock_urllib2({'mock':
                                            {'label': sp.ones((150,)),
                                             'data': sp.ones((150, 4))}})
    try:
        mock = fetch_mldata('mock', data_home=tmpdir)
        assert_in(mock, in_=['COL_NAMES', 'DESCR', 'target', 'data'])

        assert_equal(mock.target.shape, (150,))
        assert_equal(mock.data.shape, (150, 4))

        assert_raises(datasets.mldata.urllib2.HTTPError,
                      fetch_mldata, 'not_existing_name')
    finally:
        datasets.mldata.urllib2 = _urllib2_ref


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_fetch_one_column():
    _urllib2_ref = datasets.mldata.urllib2
    try:
        dataname = 'onecol'
        # create fake data set in cache
        x = sp.arange(6).reshape(2, 3)
        datasets.mldata.urllib2 = mock_urllib2({dataname: {'x': x}})

        dset = fetch_mldata(dataname, data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'data'], out_=['target'])

        assert_equal(dset.data.shape, (2, 3))
        assert_array_equal(dset.data, x)

        # transposing the data array
        dset = fetch_mldata(dataname, transpose_data=False, data_home=tmpdir)
        assert_equal(dset.data.shape, (3, 2))
    finally:
        datasets.mldata.urllib2 = _urllib2_ref


@with_setup(setup_tmpdata, teardown_tmpdata)
def test_fetch_multiple_column():
    _urllib2_ref = datasets.mldata.urllib2
    try:
        # create fake data set in cache
        x = sp.arange(6).reshape(2, 3)
        y = sp.array([1, -1])
        z = sp.arange(12).reshape(4, 3)

        # by default
        dataname = 'threecol-default'
        datasets.mldata.urllib2 = mock_urllib2({dataname:
                                                ({'label': y,
                                                  'data': x,
                                                  'z': z},
                                                 ['z', 'data', 'label'])})

        dset = fetch_mldata(dataname, data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data', 'z'],
                  out_=['x', 'y'])

        assert_array_equal(dset.data, x)
        assert_array_equal(dset.target, y)
        assert_array_equal(dset.z, z.T)

        # by order
        dataname = 'threecol-order'
        datasets.mldata.urllib2 = mock_urllib2({dataname:
                                                ({'y': y,
                                                  'x': x,
                                                  'z': z},
                                                 ['y', 'x', 'z'])})

        dset = fetch_mldata(dataname, data_home=tmpdir)
        assert_in(dset, in_=['COL_NAMES', 'DESCR', 'target', 'data', 'z'],
                  out_=['x', 'y'])
        assert_array_equal(dset.data, x)
        assert_array_equal(dset.target, y)
        assert_array_equal(dset.z, z.T)

        # by number
        dataname = 'threecol-number'
        datasets.mldata.urllib2 = mock_urllib2({dataname:
                                                ({'y': y,
                                                  'x': x,
                                                  'z': z},
                                                 ['z', 'x', 'y'])})

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
        datasets.mldata.urllib2 = _urllib2_ref
