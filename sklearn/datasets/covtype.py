"""Forest covertype dataset.

A classic dataset for classification benchmarks, featuring categorical and
real-valued features.
"""

# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: 3-clause BSD.

import errno
from gzip import GzipFile
from io import BytesIO
import logging
import os
from os.path import exists, join
from urllib2 import urlopen

import numpy as np

from .base import get_data_home
from .base import Bunch
from ..externals import joblib
from ..utils import check_random_state


URL = ('http://archive.ics.uci.edu/ml/'
       'machine-learning-databases/covtype/covtype.data.gz')


logger = logging.getLogger()


def fetch_covtype(data_home=None, download_if_missing=True,
                  random_state=None, shuffle=False):
    """Load the covertype dataset, downloading it if necessary.

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
    """

    data_home = get_data_home(data_home=data_home)
    covtype_dir = join(data_home, "covertype")
    samples_path = join(covtype_dir, "samples")
    targets_path = join(covtype_dir, "targets")
    available = exists(samples_path)

    if download_if_missing and not available:
        _mkdirp(covtype_dir)
        logger.warn("Downloading %s" % URL)
        f = BytesIO(urlopen(URL).read())
        Xy = np.genfromtxt(GzipFile(fileobj=f), delimiter=',')

        X = Xy[:, :-1]
        y = Xy[:, -1].astype(np.int32)

        joblib.dump(X, samples_path, compress=9)
        joblib.dump(y, targets_path, compress=9)

    try:
        X, y
    except NameError:
        X = joblib.load(samples_path)
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
