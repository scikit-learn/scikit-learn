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
folder.

The `fetch_20newsgroups` function will not vectorize the data into numpy
arrays but the dataset lists the filenames of the posts and their categories
as target labels.

The `fetch_20newsgroups_tfidf` function will in addition do a simple tf-idf
vectorization step.

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
import scipy.sparse as sp

from .base import get_data_home
from .base import Bunch
from .base import load_files
from ..utils import check_random_state
from ..utils.fixes import in1d
from ..feature_extraction.text import CountVectorizer
from ..preprocessing import normalize
from ..externals import joblib


logger = logging.getLogger(__name__)


URL = ("http://people.csail.mit.edu/jrennie/"
            "20Newsgroups/20news-bydate.tar.gz")
ARCHIVE_NAME = "20news-bydate.tar.gz"
CACHE_NAME = "20news-bydate.pkz"
TRAIN_FOLDER = "20news-bydate-train"
TEST_FOLDER = "20news-bydate-test"


def download_20newsgroups(target_dir, cache_path):
    """Download the 20 newsgroups data and stored it as a zipped pickle."""
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
            train=load_files(train_path, charset='latin1'),
            test=load_files(test_path, charset='latin1')
        )
    open(cache_path, 'wb').write(pickle.dumps(cache).encode('zip'))
    shutil.rmtree(target_dir)
    return cache


def fetch_20newsgroups(data_home=None, subset='train', categories=None,
                      shuffle=True, random_state=42, download_if_missing=True):
    """Load the filenames of the 20 newsgroups dataset.

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
        except Exception as e:
            print 80 * '_'
            print 'Cache loading failed'
            print 80 * '_'
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
        filenames = list()
        for subset in ('train', 'test'):
            data = cache[subset]
            data_lst.extend(data.data)
            target.extend(data.target)
            filenames.extend(data.filenames)

        data.data = data_lst
        data.target = np.array(target)
        data.filenames = np.array(filenames)
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
        data.filenames = data.filenames[mask]
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
        data.filenames = data.filenames[indices]
        data.target = data.target[indices]
        # Use an object array to shuffle: avoids memory copy
        data_lst = np.array(data.data, dtype=object)
        data_lst = data_lst[indices]
        data.data = data_lst.tolist()

    return data


def fetch_20newsgroups_vectorized(subset="train", data_home=None):
    """Load the 20 newsgroups dataset and transform it into tf-idf vectors.

    This is a convenience function; the tf-idf transformation is done using the
    default settings for `sklearn.feature_extraction.text.Vectorizer`. For more
    advanced usage (stopword filtering, n-gram extraction, etc.), combine
    fetch_20newsgroups with a custom `Vectorizer` or `CountVectorizer`.

    Parameters
    ----------

    subset: 'train' or 'test', 'all', optional
        Select the dataset to load: 'train' for the training set, 'test'
        for the test set, 'all' for both, with shuffled ordering.

    data_home: optional, default: None
        Specify an download and cache folder for the datasets. If None,
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    Returns
    -------

    bunch : Bunch object
        bunch.data: sparse matrix, shape [n_samples, n_features]
        bunch.target: array, shape [n_samples]
        bunch.target_names: list, length [n_classes]
    """
    data_home = get_data_home(data_home=data_home)
    target_file = os.path.join(data_home, "20newsgroup_vectorized.pk")

    # we shuffle but use a fixed seed for the memoization
    data_train = fetch_20newsgroups(data_home=data_home,
                                    subset='train',
                                    categories=None,
                                    shuffle=True,
                                    random_state=12)

    data_test = fetch_20newsgroups(data_home=data_home,
                                   subset='test',
                                   categories=None,
                                   shuffle=True,
                                   random_state=12)

    if os.path.exists(target_file):
        X_train, X_test = joblib.load(target_file)
    else:
        vectorizer = CountVectorizer(dtype=np.int16)
        X_train = vectorizer.fit_transform(data_train.data).tocsr()
        X_test = vectorizer.transform(data_test.data).tocsr()
        joblib.dump((X_train, X_test), target_file, compress=9)

    # the data is stored as int16 for compactness
    # but normalize needs floats
    X_train = X_train.astype(np.float64)
    X_test = X_test.astype(np.float64)
    normalize(X_train, copy=False)
    normalize(X_test, copy=False)

    target_names = data_train.target_names

    if subset == "train":
        data = X_train
        target = data_train.target
    elif subset == "test":
        data = X_test
        target = data_test.target
    elif subset == "all":
        data = sp.vstack((X_train, X_test)).tocsr()
        target = np.concatenate((data_train.target, data_test.target))
    else:
        raise ValueError("%r is not a valid subset: should be one of "
                         "['train', 'test', 'all']" % subset)

    return Bunch(data=data, target=target, target_names=target_names)
