"""Test the 20news downloader, if the data is available,
or if specifically requested via environment variable
(e.g. for travis cron job)."""
from functools import partial

import numpy as np
import scipy.sparse as sp

from sklearn.utils._testing import assert_allclose_dense_sparse
from sklearn.datasets.tests.test_common import check_return_X_y
from sklearn.preprocessing import normalize


def test_20news(fetch_20newsgroups_fxt):
    data = fetch_20newsgroups_fxt(subset='all', shuffle=False)

    # Extract a reduced dataset
    data2cats = fetch_20newsgroups_fxt(
        subset='all', categories=data.target_names[-1:-3:-1], shuffle=False)
    # Check that the ordering of the target_names is the same
    # as the ordering in the full dataset
    assert data2cats.target_names == data.target_names[-2:]
    # Assert that we have only 0 and 1 as labels
    assert np.unique(data2cats.target).tolist() == [0, 1]

    # Check that the number of filenames is consistent with data/target
    assert len(data2cats.filenames) == len(data2cats.target)
    assert len(data2cats.filenames) == len(data2cats.data)

    # Check that the first entry of the reduced dataset corresponds to
    # the first entry of the corresponding category in the full dataset
    entry1 = data2cats.data[0]
    category = data2cats.target_names[data2cats.target[0]]
    label = data.target_names.index(category)
    entry2 = data.data[np.where(data.target == label)[0][0]]
    assert entry1 == entry2

    # check that return_X_y option
    X, y = fetch_20newsgroups_fxt(subset='all', shuffle=False, return_X_y=True)
    assert len(X) == len(data.data)
    assert y.shape == data.target.shape


def test_20news_length_consistency(fetch_20newsgroups_fxt):
    """Checks the length consistencies within the bunch

    This is a non-regression test for a bug present in 0.16.1.
    """
    # Extract the full dataset
    data = fetch_20newsgroups_fxt(subset='all')
    assert len(data['data']) == len(data.data)
    assert len(data['target']) == len(data.target)
    assert len(data['filenames']) == len(data.filenames)


def test_20news_vectorized(fetch_20newsgroups_vectorized_fxt):
    # test subset = train
    bunch = fetch_20newsgroups_vectorized_fxt(subset="train")
    assert sp.isspmatrix_csr(bunch.data)
    assert bunch.data.shape == (11314, 130107)
    assert bunch.target.shape[0] == 11314
    assert bunch.data.dtype == np.float64

    # test subset = test
    bunch = fetch_20newsgroups_vectorized_fxt(subset="test")
    assert sp.isspmatrix_csr(bunch.data)
    assert bunch.data.shape == (7532, 130107)
    assert bunch.target.shape[0] == 7532
    assert bunch.data.dtype == np.float64

    # test return_X_y option
    fetch_func = partial(fetch_20newsgroups_vectorized_fxt, subset='test')
    check_return_X_y(bunch, fetch_func)

    # test subset = all
    bunch = fetch_20newsgroups_vectorized_fxt(subset='all')
    assert sp.isspmatrix_csr(bunch.data)
    assert bunch.data.shape == (11314 + 7532, 130107)
    assert bunch.target.shape[0] == 11314 + 7532
    assert bunch.data.dtype == np.float64


def test_20news_normalization(fetch_20newsgroups_vectorized_fxt):
    X = fetch_20newsgroups_vectorized_fxt(normalize=False)
    X_ = fetch_20newsgroups_vectorized_fxt(normalize=True)
    X_norm = X_['data'][:100]
    X = X['data'][:100]

    assert_allclose_dense_sparse(X_norm, normalize(X))
    assert np.allclose(np.linalg.norm(X_norm.todense(), axis=1), 1)
