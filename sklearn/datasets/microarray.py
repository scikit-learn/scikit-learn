"""Caching loader of a microarray dataset.

The data contains microarray expression levels from human B-cells. It
comes from:

* Cheng, Y., & Church, G. M. (2000, August). Biclustering of expression
  data. In Ismb (Vol. 8, pp. 93-103).

and is based on Alizadeh et al. Nature 403: 503-511 (2000).

"""

# Copyright (c) 2011 Olivier Grisel <olivier.grisel@ensta.org>
# License: BSD 3 clause

import os
import pickle

import numpy as np

from .base import get_data_home
from ..utils import check_random_state
from ..externals import six

if six.PY3:
    from urllib.request import urlopen
else:
    from urllib2 import urlopen


URL = "http://arep.med.harvard.edu/biclustering/lymphoma.matrix"
CACHE_NAME = "lymphoma-micrarray.pkz"


def download_microarray(cache_path):
    """Download the microarray data and stored it as a zipped pickle."""
    lines = urlopen(URL).read().strip().split('\n')
    # insert a space before all negative signs
    lines = list(' -'.join(line.split('-')).split(' ') for line in lines)
    lines = list(list(int(i) for i in line if i) for line in lines)
    data = np.array(lines)

    # mask missing values
    data = np.ma.masked_equal(data, 999)

    # Store a zipped pickle
    open(cache_path, 'wb').write(pickle.dumps(data).encode('zip'))
    return data


def fetch_microarray(replace_missing=True, data_home=None,
                     random_state=0, download_if_missing=True):
    """Load the microarray data.

    Returns a ``numpy.ma.MaskedArray``, with missing values maksed.

    Parameters
    ----------
    replace_missing: optional, default: True
        Whether to replace missing values with random ones.

    data_home: optional, default: None
        Specify an download and cache folder for the datasets. If None,
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    random_state: numpy random number generator or seed integer
        Used to shuffle the dataset.

    download_if_missing: optional, True by default
        If False, raise an IOError if the data is not locally available
        instead of trying to download the data from the source site.

    """

    data_home = get_data_home(data_home=data_home)
    cache_path = os.path.join(data_home, CACHE_NAME)
    data = None
    if os.path.exists(cache_path):
        try:
            data = pickle.loads(open(cache_path, 'rb').read().decode('zip'))
        except Exception as e:
            print(80 * '_')
            print('Cache loading failed')
            print(80 * '_')
            print(e)

    if data is None:
        if download_if_missing:
            data = download_microarray(cache_path)
        else:
            raise IOError('microarray dataset not found')

    if replace_missing:
        generator = check_random_state(random_state)
        data[data.mask] = generator.randint(-800, 801, data.mask.sum())

    return data
