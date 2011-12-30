"""Test the 20news downloader, if the data is available."""
import numpy as np
import scipy.sparse as sp
from nose.tools import assert_equal
from nose.tools import assert_true
from nose.plugins.skip import SkipTest

from sklearn import datasets


def test_20news():
    try:
        data = datasets.fetch_20newsgroups(subset='all',
                        download_if_missing=False,
                        shuffle=False)
    except IOError:
        raise SkipTest("Download 20 newsgroups to run this test")

    # Extract a reduced dataset
    data2cats = datasets.fetch_20newsgroups(subset='all',
                            categories=data.target_names[-1:-3:-1],
                            shuffle=False)
    # Check that the ordering of the target_names is the same
    # as the ordering in the full dataset
    assert_equal(data2cats.target_names,
                 data.target_names[-2:])
    # Assert that we have only 0 and 1 as labels
    assert_equal(np.unique(data2cats.target).tolist(), [0, 1])

    # Check that the first entry of the reduced dataset corresponds to
    # the first entry of the corresponding category in the full dataset
    entry1 = data2cats.data[0]
    category = data2cats.target_names[data2cats.target[0]]
    label = data.target_names.index(category)
    entry2 = data.data[np.where(data.target == label)[0][0]]
    assert_equal(entry1, entry2)


def test_20news_vectorized():
    # This test is slow.
    raise SkipTest

    bunch = datasets.fetch_20newsgroups_vectorized(subset="train")
    assert_true(sp.isspmatrix_csr(bunch.data))
    assert_equal(bunch.data.shape, (11314, 107428))
    assert_equal(bunch.target.shape[0], 11314)
    assert_equal(bunch.data.dtype, np.float64)

    bunch = datasets.fetch_20newsgroups_vectorized(subset="test")
    assert_true(sp.isspmatrix_csr(bunch.data))
    assert_equal(bunch.data.shape, (7532, 107428))
    assert_equal(bunch.target.shape[0], 7532)
    assert_equal(bunch.data.dtype, np.float64)

    bunch = datasets.fetch_20newsgroups_vectorized(subset="all")
    assert_true(sp.isspmatrix_csr(bunch.data))
    assert_equal(bunch.data.shape, (11314 + 7532, 107428))
    assert_equal(bunch.target.shape[0], 11314 + 7532)
    assert_equal(bunch.data.dtype, np.float64)
