"""Test the 20news downloader, if the data is available."""
import numpy as np
import scipy.sparse as sp

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import SkipTest
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial

from sklearn import datasets


def test_20news():
    try:
        data = datasets.fetch_20newsgroups(
            subset='all', download_if_missing=False, shuffle=False)
    except IOError:
        raise SkipTest("Download 20 newsgroups to run this test")

    # Extract a reduced dataset
    data2cats = datasets.fetch_20newsgroups(
        subset='all', categories=data.target_names[-1:-3:-1], shuffle=False)
    # Check that the ordering of the target_names is the same
    # as the ordering in the full dataset
    assert_equal(data2cats.target_names,
                 data.target_names[-2:])
    # Assert that we have only 0 and 1 as labels
    assert_equal(np.unique(data2cats.target).tolist(), [0, 1])

    # Check that the number of filenames is consistent with data/target
    assert_equal(len(data2cats.filenames), len(data2cats.target))
    assert_equal(len(data2cats.filenames), len(data2cats.data))

    # Check that the first entry of the reduced dataset corresponds to
    # the first entry of the corresponding category in the full dataset
    entry1 = data2cats.data[0]
    category = data2cats.target_names[data2cats.target[0]]
    label = data.target_names.index(category)
    entry2 = data.data[np.where(data.target == label)[0][0]]
    assert_equal(entry1, entry2)


def test_20news_length_consistency():
    """Checks the length consistencies within the bunch

    This is a non-regression test for a bug present in 0.16.1.
    """
    try:
        data = datasets.fetch_20newsgroups(
            subset='all', download_if_missing=False, shuffle=False)
    except IOError:
        raise SkipTest("Download 20 newsgroups to run this test")
    # Extract the full dataset
    data = datasets.fetch_20newsgroups(subset='all')
    assert_equal(len(data['data']), len(data.data))
    assert_equal(len(data['target']), len(data.target))
    assert_equal(len(data['filenames']), len(data.filenames))


def test_20news_vectorized():
    try:
        datasets.fetch_20newsgroups(subset='all',
                                    download_if_missing=False)
    except IOError:
        raise SkipTest("Download 20 newsgroups to run this test")

    # test subset = train
    bunch = datasets.fetch_20newsgroups_vectorized(subset="train")
    assert sp.isspmatrix_csr(bunch.data)
    assert_equal(bunch.data.shape, (11314, 130107))
    assert_equal(bunch.target.shape[0], 11314)
    assert_equal(bunch.data.dtype, np.float64)

    # test subset = test
    bunch = datasets.fetch_20newsgroups_vectorized(subset="test")
    assert sp.isspmatrix_csr(bunch.data)
    assert_equal(bunch.data.shape, (7532, 130107))
    assert_equal(bunch.target.shape[0], 7532)
    assert_equal(bunch.data.dtype, np.float64)

    # test return_X_y option
    fetch_func = partial(datasets.fetch_20newsgroups_vectorized, subset='test')
    check_return_X_y(bunch, fetch_func)

    # test subset = all
    bunch = datasets.fetch_20newsgroups_vectorized(subset='all')
    assert sp.isspmatrix_csr(bunch.data)
    assert_equal(bunch.data.shape, (11314 + 7532, 130107))
    assert_equal(bunch.target.shape[0], 11314 + 7532)
    assert_equal(bunch.data.dtype, np.float64)
