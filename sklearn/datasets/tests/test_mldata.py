"""Test functionality of mldata fetching utilities."""

import os
import scipy as sp
import shutil

from sklearn import datasets
from sklearn.datasets import mldata_filename, fetch_mldata

from sklearn.utils.testing import assert_in
from sklearn.utils.testing import assert_not_in
from sklearn.utils.testing import mock_mldata_urlopen
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_warns

import pytest


@pytest.fixture(scope='module')
def tmpdata(tmpdir_factory):
    tmpdir = tmpdir_factory.mktemp('tmp')
    tmpdir_path = str(tmpdir.join('mldata'))
    os.makedirs(tmpdir_path)
    yield str(tmpdir)
    shutil.rmtree(str(tmpdir))


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_mldata_filename():
    cases = [('datasets-UCI iris', 'datasets-uci-iris'),
             ('news20.binary', 'news20binary'),
             ('book-crossing-ratings-1.0', 'book-crossing-ratings-10'),
             ('Nile Water Level', 'nile-water-level'),
             ('MNIST (original)', 'mnist-original')]
    for name, desired in cases:
        assert_equal(mldata_filename(name), desired)


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_download(tmpdata):
    """Test that fetch_mldata is able to download and cache a data set."""
    _urlopen_ref = datasets.mldata.urlopen
    datasets.mldata.urlopen = mock_mldata_urlopen({
        'mock': {
            'label': sp.ones((150,)),
            'data': sp.ones((150, 4)),
        },
    })
    try:
        mock = assert_warns(DeprecationWarning, fetch_mldata,
                            'mock', data_home=tmpdata)
        for n in ["COL_NAMES", "DESCR", "target", "data"]:
            assert_in(n, mock)

        assert_equal(mock.target.shape, (150,))
        assert_equal(mock.data.shape, (150, 4))

        assert_raises(datasets.mldata.HTTPError,
                      assert_warns, DeprecationWarning,
                      fetch_mldata, 'not_existing_name')
    finally:
        datasets.mldata.urlopen = _urlopen_ref


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_fetch_one_column(tmpdata):
    _urlopen_ref = datasets.mldata.urlopen
    try:
        dataname = 'onecol'
        # create fake data set in cache
        x = sp.arange(6).reshape(2, 3)
        datasets.mldata.urlopen = mock_mldata_urlopen({dataname: {'x': x}})

        dset = fetch_mldata(dataname, data_home=tmpdata)
        for n in ["COL_NAMES", "DESCR", "data"]:
            assert_in(n, dset)
        assert_not_in("target", dset)

        assert_equal(dset.data.shape, (2, 3))
        assert_array_equal(dset.data, x)

        # transposing the data array
        dset = fetch_mldata(dataname, transpose_data=False, data_home=tmpdata)
        assert_equal(dset.data.shape, (3, 2))
    finally:
        datasets.mldata.urlopen = _urlopen_ref


@pytest.mark.filterwarnings('ignore::DeprecationWarning')
def test_fetch_multiple_column(tmpdata):
    _urlopen_ref = datasets.mldata.urlopen
    try:
        # create fake data set in cache
        x = sp.arange(6).reshape(2, 3)
        y = sp.array([1, -1])
        z = sp.arange(12).reshape(4, 3)

        # by default
        dataname = 'threecol-default'
        datasets.mldata.urlopen = mock_mldata_urlopen({
            dataname: (
                {
                    'label': y,
                    'data': x,
                    'z': z,
                },
                ['z', 'data', 'label'],
            ),
        })

        dset = fetch_mldata(dataname, data_home=tmpdata)
        for n in ["COL_NAMES", "DESCR", "target", "data", "z"]:
            assert_in(n, dset)
        assert_not_in("x", dset)
        assert_not_in("y", dset)

        assert_array_equal(dset.data, x)
        assert_array_equal(dset.target, y)
        assert_array_equal(dset.z, z.T)

        # by order
        dataname = 'threecol-order'
        datasets.mldata.urlopen = mock_mldata_urlopen({
            dataname: ({'y': y, 'x': x, 'z': z},
                       ['y', 'x', 'z']), })

        dset = fetch_mldata(dataname, data_home=tmpdata)
        for n in ["COL_NAMES", "DESCR", "target", "data", "z"]:
            assert_in(n, dset)
        assert_not_in("x", dset)
        assert_not_in("y", dset)

        assert_array_equal(dset.data, x)
        assert_array_equal(dset.target, y)
        assert_array_equal(dset.z, z.T)

        # by number
        dataname = 'threecol-number'
        datasets.mldata.urlopen = mock_mldata_urlopen({
            dataname: ({'y': y, 'x': x, 'z': z},
                       ['z', 'x', 'y']),
        })

        dset = fetch_mldata(dataname, target_name=2, data_name=0,
                            data_home=tmpdata)
        for n in ["COL_NAMES", "DESCR", "target", "data", "x"]:
            assert_in(n, dset)
        assert_not_in("y", dset)
        assert_not_in("z", dset)

        assert_array_equal(dset.data, z)
        assert_array_equal(dset.target, y)

        # by name
        dset = fetch_mldata(dataname, target_name='y', data_name='z',
                            data_home=tmpdata)
        for n in ["COL_NAMES", "DESCR", "target", "data", "x"]:
            assert_in(n, dset)
        assert_not_in("y", dset)
        assert_not_in("z", dset)

    finally:
        datasets.mldata.urlopen = _urlopen_ref
