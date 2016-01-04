"""California housing dataset.

The original database is available from StatLib

    http://lib.stat.cmu.edu/

The data contains 20,640 observations on 9 variables.

This dataset contains the average house value as target variable
and the following input variables (features): average income,
housing average age, average rooms, average bedrooms, population,
average occupation, latitude, and longitude in that order.

References
----------

Pace, R. Kelley and Ronald Barry, Sparse Spatial Autoregressions,
Statistics and Probability Letters, 33 (1997) 291-297.

"""
# Authors: Peter Prettenhofer
# License: BSD 3 clause

from io import BytesIO
import os
from os.path import exists
from os import makedirs
import tarfile

try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen

import numpy as np

from .base import get_data_home, Bunch
from .base import _pkl_filepath
from ..externals import joblib


DATA_URL = "http://www.dcc.fc.up.pt/~ltorgo/Regression/cal_housing.tgz"
TARGET_FILENAME = "cal_housing.pkz"

# Grab the module-level docstring to use as a description of the
# dataset
MODULE_DOCS = __doc__


def fetch_california_housing(data_home=None, download_if_missing=True):
    """Loader for the California housing dataset from StatLib.

    Read more in the :ref:`User Guide <datasets>`.

    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing: optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    Returns
    -------
    dataset : dict-like object with the following attributes:

    dataset.data : ndarray, shape [20640, 8]
        Each row corresponding to the 8 feature values in order.

    dataset.target : numpy array of shape (20640,)
        Each value corresponds to the average house value in units of 100,000.

    dataset.feature_names : array of length 8
        Array of ordered feature names used in the dataset.

    dataset.DESCR : string
        Description of the California housing dataset.

    Notes
    ------

    This dataset consists of 20,640 samples and 9 features.
    """
    data_home = get_data_home(data_home=data_home)
    if not exists(data_home):
        makedirs(data_home)
    filepath = _pkl_filepath(data_home, TARGET_FILENAME)
    if not exists(filepath):
        print('downloading Cal. housing from %s to %s' % (DATA_URL, data_home))
        archive_fileobj = BytesIO(urlopen(DATA_URL).read())
        fileobj = tarfile.open(
            mode="r:gz",
            fileobj=archive_fileobj).extractfile(
                'CaliforniaHousing/cal_housing.data')

        cal_housing = np.loadtxt(fileobj, delimiter=',')
        # Columns are not in the same order compared to the previous
        # URL resource on lib.stat.cmu.edu
        columns_index = [8, 7, 2, 3, 4, 5, 6, 1, 0]
        cal_housing = cal_housing[:, columns_index]
        joblib.dump(cal_housing, filepath, compress=6)
    else:
        cal_housing = joblib.load(filepath)

    feature_names = ["MedInc", "HouseAge", "AveRooms", "AveBedrms",
                     "Population", "AveOccup", "Latitude", "Longitude"]

    target, data = cal_housing[:, 0], cal_housing[:, 1:]

    # avg rooms = total rooms / households
    data[:, 2] /= data[:, 5]

    # avg bed rooms = total bed rooms / households
    data[:, 3] /= data[:, 5]

    # avg occupancy = population / housholds
    data[:, 5] = data[:, 4] / data[:, 5]

    # target in units of 100,000
    target = target / 100000.0

    return Bunch(data=data,
                 target=target,
                 feature_names=feature_names,
                 DESCR=MODULE_DOCS)
