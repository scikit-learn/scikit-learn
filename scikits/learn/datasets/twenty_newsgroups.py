"""Caching loader for the 20 newsgroups text classification dataset


The description of the dataset is available on the official website at:

    http://people.csail.mit.edu/jrennie/20Newsgroups/

Quoting the introduction:

    The 20 Newsgroups data set is a collection of approximately 20,000
    newsgroup documents, partitioned (nearly) evenly across 20 different
    newsgroups. To the best of my knowledge, it was originally collected
    by Ken Lang, probably for his Newsweeder: Learning to filter netnews
    paper, though he does not explicitly mention this collection. The 20
    newsgroups collection has become a popular data set for experiments
    in text applications of machine learning techniques, such as text
    classification and text clustering.

This dataset loader will download the recommended "by date" variant of the
dataset and which features a point in time split between the train and
test sets. The compressed dataset size is around 14 Mb compressed. Once
uncompressed the train set is 52 MB and the test set is 34 MB.

The data is downloaded, extracted and cached in the '~/scikit_learn_data'
folder. However contrary to other datasets in the scikit, the data is
not vectorized into numpy arrays but the dataset list the filenames of
the posts and there categories as target signal.

The lack of vector feature extraction is intentional: there is no single
best way to turn text into vectors. Depending on the task various
preprocessing and text transformation are useful or not (n-grams,
lowercasing, stemming, stop-words filtering, TF-IDF weighting...).

"""
# Copyright (c) 2011 Olivier Grisel <olivier.grisel@ensta.org>
# License: Simplified BSD

import os
import urllib
import logging
import tarfile

from .base import get_data_home
from .base import load_filenames


URL = ("http://people.csail.mit.edu/jrennie/"
            "20Newsgroups/20news-bydate.tar.gz")
ARCHIVE_NAME = "20news-bydate.tar.gz"
TRAIN_FOLDER = "20news-bydate-train"
TEST_FOLDER = "20news-bydate-test"


def fetch_20newsgroups(data_home=None, subset='train', categories=None,
                      shuffle=True, rng=42, download_if_missing=True):
    """Load the filenames of the 20 newsgroups dataset

    Parameters
    ----------
    subset: optional, default: 'train'
        Select the dataset to load: 'train' for the training set, 'test'
        for the test set.

    data_home: optional, default: None
        Specify another download and cache folder for the datasets. By
        default all scikit learn data is stored in '~/scikit_learn_data'
        subfolders.

    categories: None or collection of string or unicode
        if None (default), load all the categories.
        if not Non, list of category names to load (other categories
        ignored)

    shuffle: True by default
        whether or not to shuffle the data: might be important for models
        that make the assumption that the samples are independent and
        identically distributed (i.i.d.) such as stochastic gradient
        descent for instance.

    rng: a numpy random number generator or a seed integer, 42 by default
        used to shuffle the dataset

    download_if_missing: optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.
    """

    data_home = get_data_home(data_home=data_home)
    twenty_home = os.path.join(data_home, "20news_home")
    archive_path = os.path.join(twenty_home, ARCHIVE_NAME)
    train_path = os.path.join(twenty_home, TRAIN_FOLDER)
    test_path = os.path.join(twenty_home, TEST_FOLDER)

    if not os.path.exists(twenty_home):
        os.makedirs(twenty_home)

    if not os.path.exists(train_path) or not os.path.exists(test_path):

        if not os.path.exists(archive_path):
            if download_if_missing:
                logging.warn("Downloading dataset from %s (14 MB)", URL)
                opener = urllib.urlopen(URL)
                open(archive_path, 'wb').write(opener.read())
            else:
                raise IOError("%s is missing" % archive_path)

        logging.info("Decompressing %s", archive_path)
        tarfile.open(archive_path, "r:gz").extractall(path=twenty_home)
        os.remove(archive_path)

    if subset == 'train':
        folder_path = train_path
    elif subset == 'test':
        folder_path = test_path
    else:
        raise ValueError(
            "subset can only be 'train' or 'test', got '%s'" % subset)

    description = subset + ' subset of the 20 newsgroups by date dataset'
    return load_filenames(folder_path, description=description,
                          categories=categories, shuffle=shuffle, rng=rng)


def load_20newsgroups(download_if_missing=False, **kwargs):
    """Alias for fetch_20newsgroups(download_if_missing=False)

    Check out fetch_20newsgroups.__doc__ for the documentation and parameters
    list.
    """
    return fetch_20newsgroups(download_if_missing=download_if_missing, **kwargs)
