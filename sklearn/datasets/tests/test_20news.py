"""Test the 20news downloader, if the data is available."""
import numpy as np
from nose.tools import assert_equal
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
    categories = ['alt.atheism', 'talk.religion.misc']
    try:
        X_train, y_train, X_test, y_test = \
            datasets.load_vectorized_20newsgroups(download_if_missing=False,
                                                  categories=categories,
                                                  shuffle=False)
    except IOError:
        raise SkipTest("Download 20 newsgroups to run this test")

    assert_equal(X_train.shape, (857, 16739))
    assert_equal(y_train.shape[0], 857)
    assert_equal(X_test.shape, (570, 16739))
    assert_equal(y_test.shape[0], 570)
