"""Adult Census dataset.

The original database was provided by

    Ronny Kohavi and Barry Becker
    Data Mining and Visualization
    Silicon Graphics.

The version retrieved here comes from

    https://archive.ics.uci.edu/ml/datasets/Adult

Extraction was done by Barry Becker from the 1994 Census database.
A set of reasonably clean records was extracted using the following
conditions: ((AAGE>16) && (AGI>100) && (AFNLWGT>1)&& (HRSWK>0)).
Prediction task is to determine whether a person makes over 50K a year.

The retrieved dataset consisted of 48842 x 15
"""
# Authors: Herilalaina Rakotoarison
# License: BSD 3 clause

from io import BytesIO
from os.path import exists
from os import makedirs

try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen

import numpy as np

from .base import get_data_home
from ..utils import Bunch
from .base import _pkl_filepath
from ..externals import joblib


DATA_URL = "https://archive.ics.uci.edu/ml/machine-learning-databases/" \
           "adult/adult.data"
TARGET_FILENAME = "adult.data"

# Grab the module-level docstring to use as a description of the
# dataset
MODULE_DOCS = __doc__


def fetch_adult_census(data_home=None, download_if_missing=True):
    """Loader for the Adult census from Silicon Graphics.

    Read more in the :ref:`User Guide <adult_census>`.

    Parameters
    ----------
    data_home : optional, default: None
        Specify another download and cache folder for the datasets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    download_if_missing : optional, True by default
        If False, raise a IOError if the data is not locally available
        instead of trying to download the data from the source site.

    Returns
    -------
    dataset : dict-like object with the following attributes:

    dataset.data : ndarray, shape [48842, 14]
        Each row corresponding to the 14 feature values in order.

    dataset.target : numpy array of shape (48842,)
        Each value corresponds whether a person makes over 50K a year.

    dataset.feature_names : array of length 14
        Array of ordered feature names used in the dataset.

    dataset.DESCR : string
        Description of the Adult Census dataset.

    Notes
    ------

    This dataset consists of 48842 samples and 15 features.
    """
    data_home = get_data_home(data_home=data_home)
    if not exists(data_home):
        makedirs(data_home)

    feature_names = ["Age", "WorkClass", "Fnlwgt", "Education", "EducationNum",
                     "MaritalStatus", "Occupation", "Relationship", "Race",
                     "Sex", "CapitalGain", "CapitalLoss", "hoursPerWeek",
                     "NativeCountry"]
    target_name = "Class"

    filepath = _pkl_filepath(data_home, TARGET_FILENAME)
    if not exists(filepath):
        if not download_if_missing:
            raise IOError("Data not found and `download_if_missing` is False")

        print('downloading Adult census from %s to %s' % (DATA_URL, data_home))
        fileobj = BytesIO(urlopen(DATA_URL).read())
        adult_census = np.genfromtxt(fileobj, delimiter=',', dtype=None,
                                     names=feature_names + [target_name])
        joblib.dump(adult_census, filepath, compress=6)
    else:
        adult_census = joblib.load(filepath)

    target, data = adult_census[target_name], adult_census[feature_names]

    return Bunch(data=data,
                 target=target,
                 feature_names=feature_names,
                 DESCR=MODULE_DOCS)
