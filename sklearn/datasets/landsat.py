"""Landsat dataset.

One of the datasets used for comparison of different classification algorithms
in StatLog project. It contains 6435 samples with 36 dimensions. The task is to
predict one of 6 class labels.

The dataset was created from a part of the image of agricultural land in
Australia taken by Landsat satellite. Each pixel has a class label
corresponding to one of 6 types of terrain. And it is described by
gray-scale values of pixels in 3x3 neighborhood measured in 4 different
spectral bands (thus 36 features in total).

The dataset is available from UCI Machine Learning Repository

    https://archive.ics.uci.edu/ml/datasets/Statlog+(Landsat+Satellite)
"""

# Author: Nikolay Mayorov <n59_ru@hotmail.com> (based on covtype.py)
# License: BSD 3 clause

import sys
import errno
from io import BytesIO
import logging
import os
from os.path import exists, join
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

import numpy as np

from .base import get_data_home
from .base import Bunch
from ..externals import joblib
from ..utils import check_random_state


FOLDER_URL = ("https://archive.ics.uci.edu/"
              "ml/machine-learning-databases/statlog/satimage/")
TRAIN_URL = FOLDER_URL + "sat.trn"
TEST_URL = FOLDER_URL + "sat.tst"

logger = logging.getLogger()


def fetch_landsat(data_home=None, download_if_missing=True,
                  random_state=None, shuffle=False):
    """Load Landsat dataset, downloading it if necessary.

    Parameters
    ----------
    data_home : string, optional
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing : boolean, default=True
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    random_state : int, RandomState instance or None, optional (default=None)
        Random state for shuffling the dataset.
        If int, random_state is the seed used by the random number generator;
        If RandomState instance, random_state is the random number generator;
        If None, the random number generator is the RandomState instance used
        by `np.random`.

    shuffle : bool, default=False
        Whether to shuffle dataset.

    Returns
    -------
    dataset : dict-like object with the following attributes:

    dataset.data : array, shape (6435, 36)
        Each row corresponds to the 36 features in the dataset.

    dataset.target : array, shape (6435,)
        Each value corresponds to one of the 6 types of terrain. These types
        are coded by labels from [1, 2, 3, 4, 5, 7] (6 is missing).

    dataset.DESCR : string
        Description of the landsat satellite dataset.

    """

    data_home = get_data_home(data_home=data_home)
    if sys.version_info[0] == 3:
        # The zlib compression format use by joblib is not compatible when
        # switching from Python 2 to Python 3, let us use a separate folder
        # under Python 3:
        dir_suffix = "-py3"
    else:
        # Backward compat for Python 2 users
        dir_suffix = ""
    landsat_dir = join(data_home, "landsat" + dir_suffix)
    data_path = join(landsat_dir, "data")
    targets_path = join(landsat_dir, "targets")
    available = exists(data_path) and exists(targets_path)

    if download_if_missing and not available:
        _mkdirp(landsat_dir)
        logger.warning("Downloading %s" % TRAIN_URL)
        f = BytesIO(urlopen(TRAIN_URL).read())
        Xy = np.genfromtxt(f)
        logger.warning("Downloading %s" % TEST_URL)
        f = BytesIO(urlopen(TEST_URL).read())
        Xy = np.vstack((Xy, np.genfromtxt(f)))
        X = Xy[:, :-1]
        y = Xy[:, -1].astype(np.int32)
        joblib.dump(X, data_path, compress=9)
        joblib.dump(y, targets_path, compress=9)
    try:
        X, y
    except NameError:
        X = joblib.load(data_path)
        y = joblib.load(targets_path)

    if shuffle:
        ind = np.arange(X.shape[0])
        rng = check_random_state(random_state)
        rng.shuffle(ind)
        X = X[ind]
        y = y[ind]

    return Bunch(data=X, target=y, DESCR=__doc__)


def _mkdirp(d):
    """Ensure directory d exists (like mkdir -p on Unix)
    No guarantee that the directory is writable.
    """
    try:
        os.makedirs(d)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
