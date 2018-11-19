"""Forest covertype dataset.

A classic dataset for classification benchmarks, featuring categorical and
real-valued features.

The dataset page is available from UCI Machine Learning Repository

    http://archive.ics.uci.edu/ml/datasets/Covertype

Courtesy of Jock A. Blackard and Colorado State University.
"""

# Author: Lars Buitinck
#         Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD 3 clause

from gzip import GzipFile
import logging
from os.path import dirname, exists, join
from os import remove

import numpy as np

from .base import get_data_home
from .base import _fetch_remote
from .base import RemoteFileMetadata
from ..utils import Bunch
from .base import _pkl_filepath
from ..utils.fixes import makedirs
from ..utils import _joblib
from ..utils import check_random_state

# The original data can be found in:
# http://archive.ics.uci.edu/ml/machine-learning-databases/covtype/covtype.data.gz
ARCHIVE = RemoteFileMetadata(
    filename='covtype.data.gz',
    url='https://ndownloader.figshare.com/files/5976039',
    checksum=('614360d0257557dd1792834a85a1cdeb'
              'fadc3c4f30b011d56afee7ffb5b15771'))

logger = logging.getLogger(__name__)


def fetch_covtype(data_home=None, download_if_missing=True,
                  random_state=None, shuffle=False, return_X_y=False):
    """Load the covertype dataset (classification).

    Download it if necessary.

    =================   ============
    Classes                        7
    Samples total             581012
    Dimensionality                54
    Features                     int
    =================   ============

    Read more in the :ref:`User Guide <covtype_dataset>`.

    Parameters
    ----------
    data_home : string, optional
        Specify another download and cache folder for the datasets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing : boolean, default=True
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    random_state : int, RandomState instance or None (default)
        Determines random number generation for dataset shuffling. Pass an int
        for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    shuffle : bool, default=False
        Whether to shuffle dataset.

    return_X_y : boolean, default=False.
        If True, returns ``(data.data, data.target)`` instead of a Bunch
        object.

        .. versionadded:: 0.20

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

    (data, target) : tuple if ``return_X_y`` is True

        .. versionadded:: 0.20
    """

    data_home = get_data_home(data_home=data_home)
    covtype_dir = join(data_home, "covertype")
    samples_path = _pkl_filepath(covtype_dir, "samples")
    targets_path = _pkl_filepath(covtype_dir, "targets")
    available = exists(samples_path)

    if download_if_missing and not available:
        if not exists(covtype_dir):
            makedirs(covtype_dir)
        logger.info("Downloading %s" % ARCHIVE.url)

        archive_path = _fetch_remote(ARCHIVE, dirname=covtype_dir)
        Xy = np.genfromtxt(GzipFile(filename=archive_path), delimiter=',')
        # delete archive
        remove(archive_path)

        X = Xy[:, :-1]
        y = Xy[:, -1].astype(np.int32)

        _joblib.dump(X, samples_path, compress=9)
        _joblib.dump(y, targets_path, compress=9)

    elif not available and not download_if_missing:
        raise IOError("Data not found and `download_if_missing` is False")
    try:
        X, y
    except NameError:
        X = _joblib.load(samples_path)
        y = _joblib.load(targets_path)

    if shuffle:
        ind = np.arange(X.shape[0])
        rng = check_random_state(random_state)
        rng.shuffle(ind)
        X = X[ind]
        y = y[ind]

    module_path = dirname(__file__)
    with open(join(module_path, 'descr', 'covtype.rst')) as rst_file:
        fdescr = rst_file.read()

    if return_X_y:
        return X, y

    return Bunch(data=X, target=y, DESCR=fdescr)
