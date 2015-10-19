"""Forest covertype dataset.

A classic dataset for classification benchmarks, featuring categorical and
real-valued features.

The dataset page is available from UCI Machine Learning Repository

    http://archive.ics.uci.edu/ml/datasets/Covertype

Courtesy of Jock A. Blackard and Colorado State University.
"""

# Author: Lars Buitinck <L.J.Buitinck@uva.nl>
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD 3 clause

from gzip import GzipFile
from io import BytesIO
import logging
from os.path import exists, join
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

import numpy as np

from .base import get_data_home
from .base import Bunch
from .base import _pkl_filepath
from ..utils.fixes import makedirs
from ..externals import joblib
from ..utils import check_random_state


URL = ('http://archive.ics.uci.edu/ml/'
       'machine-learning-databases/covtype/covtype.data.gz')


logger = logging.getLogger()


def fetch_covtype(data_home=None, download_if_missing=True,
                  random_state=None, shuffle=False):
    """Load the covertype dataset, downloading it if necessary.

    Read more in the :ref:`User Guide <datasets>`.

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

    dataset.data : numpy array of shape (581012, 54)
        Each row corresponds to the 54 features in the dataset.

    dataset.target : numpy array of shape (581012,)
        Each value corresponds to one of the 7 forest covertypes with values
        ranging between 1 to 7.

    dataset.DESCR : string
        Description of the forest covertype dataset.

    """

    data_home = get_data_home(data_home=data_home)
    covtype_dir = join(data_home, "covertype")
    samples_path = _pkl_filepath(covtype_dir, "samples")
    targets_path = _pkl_filepath(covtype_dir, "targets")
    available = exists(samples_path)

    if download_if_missing and not available:
        makedirs(covtype_dir, exist_ok=True)
        logger.warning("Downloading %s" % URL)
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
