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
import pickle
import shutil

import numpy as np

from .base import get_data_home
from .base import load_files
from ..utils import check_random_state, deprecated
from ..utils.fixes import in1d


logger = logging.getLogger(__name__)


URL = ("http://people.csail.mit.edu/jrennie/"
            "20Newsgroups/20news-bydate.tar.gz")
ARCHIVE_NAME = "20news-bydate.tar.gz"
CACHE_NAME = "20news-bydate.pkz"
TRAIN_FOLDER = "20news-bydate-train"
TEST_FOLDER = "20news-bydate-test"

def download_20newsgroups(target_dir, cache_path):
    """ Download the 20Newsgroups data and convert is in a zipped pickle
        storage.
    """
    archive_path = os.path.join(target_dir, ARCHIVE_NAME)
    train_path = os.path.join(target_dir, TRAIN_FOLDER)
    test_path = os.path.join(target_dir, TEST_FOLDER)

    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    if not os.path.exists(archive_path):
        logger.warn("Downloading dataset from %s (14 MB)", URL)
        opener = urllib.urlopen(URL)
        open(archive_path, 'wb').write(opener.read())

    logger.info("Decompressing %s", archive_path)
    tarfile.open(archive_path, "r:gz").extractall(path=target_dir)
    os.remove(archive_path)

    # Store a zipped pickle
    cache = dict(
            train=load_files(train_path),
            test=load_files(test_path)
        )
    open(cache_path, 'wb').write(pickle.dumps(cache).encode('zip'))
    shutil.rmtree(target_dir)
    return cache


def fetch_20newsgroups(data_home=None, subset='train', categories=None,
                      shuffle=True, random_state=42, download_if_missing=True):
    """Load the filenames of the 20 newsgroups dataset

    Parameters
    ----------
    subset: 'train' or 'test', 'all', optional
        Select the dataset to load: 'train' for the training set, 'test'
        for the test set, 'all' for both, with shuffled ordering.

    data_home: optional, default: None
        Specify an download and cache folder for the datasets. If None,
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    categories: None or collection of string or unicode
        If None (default), load all the categories.
        If not None, list of category names to load (other categories
        ignored).

    shuffle: bool, optional
        Whether or not to shuffle the data: might be important for models that
        make the assumption that the samples are independent and identically
        distributed (i.i.d.), such as stochastic gradient descent.

    random_state: numpy random number generator or seed integer
        Used to shuffle the dataset.

    download_if_missing: optional, True by default
        If False, raise an IOError if the data is not locally available
        instead of trying to download the data from the source site.
    """

    data_home = get_data_home(data_home=data_home)
    cache_path = os.path.join(data_home, CACHE_NAME)
    twenty_home = os.path.join(data_home, "20news_home")
    cache = None
    if os.path.exists(cache_path):
        try:
            cache = pickle.loads(open(cache_path, 'rb').read().decode('zip'))
        except Exception, e:
            print 80*'_'
            print 'Cache loading failed'
            print 80*'_'
            print e

    if cache is None:
        if download_if_missing:
            cache = download_20newsgroups(target_dir=twenty_home,
                                          cache_path=cache_path)
        else:
            raise IOError('20Newsgroups dataset not found')

    if subset in ('train', 'test'):
        data = cache[subset]
    elif subset == 'all':
        data_lst = list()
        target = list()
        for subset in ('train', 'test'):
            data = cache[subset]
            data_lst.extend(data.data)
            target.extend(data.target)

        data.data = data_lst
        data.target = np.array(target)
        data.description = 'the 20 newsgroups by date dataset'
    else:
        raise ValueError(
            "subset can only be 'train', 'test' or 'all', got '%s'" % subset)

    if categories is not None:
        labels = [(data.target_names.index(cat), cat) for cat in categories]
        # Sort the categories to have the ordering of the labels
        labels.sort()
        labels, categories = zip(*labels)
        mask = in1d(data.target, labels)
        data.target = data.target[mask]
        # searchsorted to have continuous labels
        data.target = np.searchsorted(labels, data.target)
        data.target_names = list(categories)
        # Use an object array to shuffle: avoids memory copy
        data_lst = np.array(data.data, dtype=object)
        data_lst = data_lst[mask]
        data.data = data_lst.tolist()

    if shuffle:
        random_state = check_random_state(random_state)
        indices = np.arange(data.target.shape[0])
        random_state.shuffle(indices)
        data.target = data.target[indices]
        # Use an object array to shuffle: avoids memory copy
        data_lst = np.array(data.data, dtype=object)
        data_lst = data_lst[indices]
        data.data = data_lst.tolist()

    return data


@deprecated("Use fetch_20newsgroups instead with download_if_missing=False")
def load_20newsgroups(download_if_missing=False, **kwargs):
    """Alias for fetch_20newsgroups(download_if_missing=False).

    See fetch_20newsgroups.__doc__ for documentation and parameter list.
    """
    return fetch_20newsgroups(download_if_missing=download_if_missing, **kwargs)
