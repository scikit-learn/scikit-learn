"""Test the 20news downloader, if the data is available,
or if specifically requested via environment variable
(e.g. for CI jobs)."""

import codecs
import os
import pickle
from functools import partial
from unittest.mock import patch

import numpy as np
import pytest
import scipy.sparse as sp

from sklearn.datasets import fetch_20newsgroups
from sklearn.datasets._base import _pkl_filepath
from sklearn.datasets._twenty_newsgroups import CACHE_NAME, _download_20newsgroups
from sklearn.datasets.tests.test_common import (
    check_as_frame,
    check_pandas_dependency_message,
    check_return_X_y,
)
from sklearn.preprocessing import normalize
from sklearn.utils import Bunch
from sklearn.utils._testing import assert_allclose_dense_sparse


def test_20news(fetch_20newsgroups_fxt):
    data = fetch_20newsgroups_fxt(subset="all", shuffle=False)
    assert data.DESCR.startswith(".. _20newsgroups_dataset:")

    # Extract a reduced dataset
    data2cats = fetch_20newsgroups_fxt(
        subset="all", categories=data.target_names[-1:-3:-1], shuffle=False
    )
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
    X, y = fetch_20newsgroups_fxt(subset="all", shuffle=False, return_X_y=True)
    assert len(X) == len(data.data)
    assert y.shape == data.target.shape


def test_20news_length_consistency(fetch_20newsgroups_fxt):
    """Checks the length consistencies within the bunch

    This is a non-regression test for a bug present in 0.16.1.
    """
    # Extract the full dataset
    data = fetch_20newsgroups_fxt(subset="all")
    assert len(data["data"]) == len(data.data)
    assert len(data["target"]) == len(data.target)
    assert len(data["filenames"]) == len(data.filenames)


def test_20news_vectorized(fetch_20newsgroups_vectorized_fxt):
    # test subset = train
    bunch = fetch_20newsgroups_vectorized_fxt(subset="train")
    assert sp.issparse(bunch.data) and bunch.data.format == "csr"
    assert bunch.data.shape == (11314, 130107)
    assert bunch.target.shape[0] == 11314
    assert bunch.data.dtype == np.float64
    assert bunch.DESCR.startswith(".. _20newsgroups_dataset:")

    # test subset = test
    bunch = fetch_20newsgroups_vectorized_fxt(subset="test")
    assert sp.issparse(bunch.data) and bunch.data.format == "csr"
    assert bunch.data.shape == (7532, 130107)
    assert bunch.target.shape[0] == 7532
    assert bunch.data.dtype == np.float64
    assert bunch.DESCR.startswith(".. _20newsgroups_dataset:")

    # test return_X_y option
    fetch_func = partial(fetch_20newsgroups_vectorized_fxt, subset="test")
    check_return_X_y(bunch, fetch_func)

    # test subset = all
    bunch = fetch_20newsgroups_vectorized_fxt(subset="all")
    assert sp.issparse(bunch.data) and bunch.data.format == "csr"
    assert bunch.data.shape == (11314 + 7532, 130107)
    assert bunch.target.shape[0] == 11314 + 7532
    assert bunch.data.dtype == np.float64
    assert bunch.DESCR.startswith(".. _20newsgroups_dataset:")


def test_20news_normalization(fetch_20newsgroups_vectorized_fxt):
    X = fetch_20newsgroups_vectorized_fxt(normalize=False)
    X_ = fetch_20newsgroups_vectorized_fxt(normalize=True)
    X_norm = X_["data"][:100]
    X = X["data"][:100]

    assert_allclose_dense_sparse(X_norm, normalize(X))
    assert np.allclose(np.linalg.norm(X_norm.todense(), axis=1), 1)


def test_20news_as_frame(fetch_20newsgroups_vectorized_fxt):
    pd = pytest.importorskip("pandas")

    bunch = fetch_20newsgroups_vectorized_fxt(as_frame=True)
    check_as_frame(bunch, fetch_20newsgroups_vectorized_fxt)

    frame = bunch.frame
    assert frame.shape == (11314, 130108)
    assert all([isinstance(col, pd.SparseDtype) for col in bunch.data.dtypes])

    # Check a small subset of features
    for expected_feature in [
        "beginner",
        "beginners",
        "beginning",
        "beginnings",
        "begins",
        "begley",
        "begone",
    ]:
        assert expected_feature in frame.keys()
    assert "category_class" in frame.keys()
    assert bunch.target.name == "category_class"


def test_as_frame_no_pandas(fetch_20newsgroups_vectorized_fxt, hide_available_pandas):
    check_pandas_dependency_message(fetch_20newsgroups_vectorized_fxt)


def test_outdated_pickle(fetch_20newsgroups_vectorized_fxt):
    with patch("os.path.exists") as mock_is_exist:
        with patch("joblib.load") as mock_load:
            # mock that the dataset was cached
            mock_is_exist.return_value = True
            # mock that we have an outdated pickle with only X and y returned
            mock_load.return_value = ("X", "y")
            err_msg = "The cached dataset located in"
            with pytest.raises(ValueError, match=err_msg):
                fetch_20newsgroups_vectorized_fxt(as_frame=True)


def _write_fake_20newsgroups_cache(data_home):
    """Write a minimal but structurally valid 20newsgroups cache to disk.

    This lets tests exercise `fetch_20newsgroups`'s post-processing (e.g. the
    `subset="all"` merging logic) without requiring network access.
    """
    cache = {
        "train": Bunch(
            data=["train doc 0", "train doc 1"],
            target=np.array([0, 1]),
            filenames=np.array(["train0", "train1"]),
            target_names=["alt.atheism", "sci.space"],
        ),
        "test": Bunch(
            data=["test doc 0"],
            target=np.array([1]),
            filenames=np.array(["test0"]),
            target_names=["alt.atheism", "sci.space"],
        ),
    }
    cache_path = _pkl_filepath(data_home, CACHE_NAME)
    compressed_content = codecs.encode(pickle.dumps(cache), "zlib_codec")
    with open(cache_path, "wb") as f:
        f.write(compressed_content)
    return cache


def test_20news_subset_all_merges_train_and_test(tmp_path):
    """Non-regression test for `fetch_20newsgroups(subset="all")`.

    A previous change to the `subset == "all"` branch introduced a loop
    variable that shadowed, without updating, a reference to the outer
    `subset` parameter, which made every `subset="all"` call raise
    `KeyError: 'all'`. This avoids the network dependency of `test_20news`
    so it always runs in CI.
    """
    _write_fake_20newsgroups_cache(tmp_path)

    data = fetch_20newsgroups(
        data_home=tmp_path,
        subset="all",
        shuffle=False,
        download_if_missing=False,
    )

    assert len(data.data) == 3
    assert list(data.data) == ["train doc 0", "train doc 1", "test doc 0"]
    assert_allclose_dense_sparse(data.target, np.array([0, 1, 1]))
    assert list(data.filenames) == ["train0", "train1", "test0"]


def test_download_20newsgroups_atomic_write(tmp_path, monkeypatch):
    """`_download_20newsgroups` must never leave a partially written cache.

    Non-regression test for a race condition (gh-32095) where concurrent
    callers (e.g. separate pytest-xdist workers) could observe a partially
    written cache file. The fix downloads and builds the cache in an
    isolated temporary directory and finalizes it with an atomic
    `os.replace`, mirroring the pattern already used by `fetch_openml` and
    `fetch_covtype`.
    """
    cache_path = _pkl_filepath(str(tmp_path), CACHE_NAME)

    fake_cache = {
        "train": Bunch(data=["doc"], target=np.array([0]), filenames=np.array(["f"])),
        "test": Bunch(data=["doc"], target=np.array([0]), filenames=np.array(["f"])),
    }

    def fake_fetch_remote(archive, dirname, n_retries, delay):
        # Simulate a downloaded archive without touching the network.
        archive_path = os.path.join(dirname, archive.filename)
        with open(archive_path, "wb") as f:
            f.write(b"not a real tar.gz, extraction is mocked below")
        return archive_path

    class _FakeTarfile:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    def fake_tarfile_open(*args, **kwargs):
        return _FakeTarfile()

    def fake_tarfile_extractall(fp, path):
        # Create the train/test folders that load_files below expects.
        os.makedirs(os.path.join(path, "20news-bydate-train"), exist_ok=True)
        os.makedirs(os.path.join(path, "20news-bydate-test"), exist_ok=True)

    def fake_load_files(path, encoding=None):
        return fake_cache["train"] if "train" in path else fake_cache["test"]

    monkeypatch.setattr(
        "sklearn.datasets._twenty_newsgroups._fetch_remote", fake_fetch_remote
    )
    monkeypatch.setattr(
        "sklearn.datasets._twenty_newsgroups.tarfile.open", fake_tarfile_open
    )
    monkeypatch.setattr(
        "sklearn.datasets._twenty_newsgroups.tarfile_extractall",
        fake_tarfile_extractall,
    )
    monkeypatch.setattr(
        "sklearn.datasets._twenty_newsgroups.load_files", fake_load_files
    )

    result = _download_20newsgroups(cache_path=cache_path, n_retries=3, delay=1)

    assert result["train"] is fake_cache["train"]
    assert os.path.exists(cache_path)
    # No leftover temporary directories or files should remain in data_home.
    assert os.listdir(tmp_path) == [os.path.basename(cache_path)]

    # The cache file that was written is a valid, fully-formed pickle.
    with open(cache_path, "rb") as f:
        uncompressed = codecs.decode(f.read(), "zlib_codec")
    reloaded = pickle.loads(uncompressed)
    assert set(reloaded) == {"train", "test"}
