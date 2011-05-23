"""Test functionality of mldata fetching utilities."""

from scikits.learn.datasets import mldata_filename, fetch_mldata
from nose.tools import assert_equal
import os
import shutil
import tempfile

def setup_module():
    pass

def teardown_module():
    pass

def test_mldata_filename():
    cases = [('datasets-UCI iris', 'datasets-uci-iris'),
             ('news20.binary', 'news20binary'),
             ('book-crossing-ratings-1.0', 'book-crossing-ratings-10'),
             ('Nile Water Level', 'nile-water-level'),
             ('MNIST (original)', 'mnist-original')]
    for name, desired in cases:
        assert_equal(mldata_filename(name), desired)

def test_fetch():
    tmpdir = tempfile.mkdtemp()
    print '----', tmpdir

    dataname = 'iris'
    iris = fetch_mldata(dataname, data_home=tmpdir)
    # cache created
    assert os.path.exists(os.path.join(tmpdir, 'mldata'))
    assert os.path.exists(os.path.join(tmpdir, 'mldata', 'iris.mat'))

    assert 'COL_NAMES' in iris
    assert 'DESCR' in iris
    assert 'target' in iris
    assert 'data' in iris

    assert_equal(iris.target.shape, (150,))
    assert_equal(iris.data.shape, (150, 4))

    # transposing the data array
    iris = fetch_mldata(dataname, transpose_data=False, data_home=tmpdir)
    assert_equal(iris.data.shape, (4, 150))

    # non-standard columns
    iris = fetch_mldata('datasets-UCI iris', target_name=1,
                        data_name=0, data_home=tmpdir)
    assert 'target' in iris
    assert 'data' in iris
    assert_equal(iris.target.shape, (150,))
    assert_equal(iris.data.shape, (150, 4))

    iris = fetch_mldata('datasets-UCI iris', target_name='class',
                        data_name='double0', data_home=tmpdir)
    assert 'target' in iris
    assert 'data' in iris
    assert_equal(iris.target.shape, (150,))
    assert_equal(iris.data.shape, (150, 4))

    # remove temporary dir
    shutil.rmtree(tmpdir)
