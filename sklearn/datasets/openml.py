import json
import numbers
import os
from os.path import join, exists
from http.client import IncompleteRead

try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen

from scipy.io.arff import loadarff
import numpy as np

from .base import get_data_home
# from ..externals.joblib import Memory
from ..externals.six import StringIO
from ..utils import Bunch

_SEARCH_NAME = "https://openml.org/api/v1/json/data/list/data_name/{}/limit/1"
_DATA_INFO = "https://openml.org/api/v1/json/data/{}"
_DATA_DOWNLOAD = "https://www.openml.org/data/download/{}"


def _get_data_info_by_name(name):
    json_string = urlopen(_SEARCH_NAME.format(name))
    json_data = json.load(json_string)
    return json_data['data']['dataset'][0]


def _get_data_description_by_id(data_id):
    json_string = urlopen(_DATA_INFO.format(data_id))
    json_data = json.load(json_string)
    return json_data['data_set_description']


def fetch_openml(name_or_id=None, version=1, data_home=None, memory=True):
    """Fetch dataset from openml by name or dataset id.

    Parameters
    ----------

    data_home : optional, default: None
        Specify another download and cache folder for the data sets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    Returns
    -------

    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the classification labels,
        'DESCR', the full description of the dataset, and
        'COL_NAMES', the original names of the dataset columns.
    """
    data_home = get_data_home(data_home=data_home)
    data_home = join(data_home, 'openml')
    # if memory:
    #     mem = Memory(join(data_home, 'cache'))
    #     _get_data_info_by_name = mem(_get_data_info_by_name)
    #     _get_data_description_by_id = mem(_get_data_description_by_id)
    #     _download_data = mem(_download_data)
    if not exists(data_home):
        os.makedirs(data_home)

    # check if dataset id is known

    if isinstance(name_or_id, numbers.Integral):
        if version != 1:
            raise ValueError(
                "Dataset id={} and version={} passed, but you can only "
                "specify a numeric id or a version, not both.".format(
                    name_or_id, version))
        data_id = name_or_id
    elif isinstance(name_or_id, str):
        data_info = _get_data_info_by_name(name_or_id)
        data_id = data_info['did']

    else:
        raise TypeError(
            "Invalid name_or_id {}, should be string or integer.".format(
                name_or_id))

    data_description = _get_data_description_by_id(data_id)
    target_name = data_description['default_target_attribute']

    # download actual data
    response = urlopen(_DATA_DOWNLOAD.format(data_id))
    # we need to catch IncompleteRead which is likely a server-side issue
    try:
        data_arff = response.read()
    except IncompleteRead as e:
        data_arff = e.partial
    # getting structured array and metadata
    data, meta = loadarff(StringIO(data_arff.decode("utf-8")))
    columns = np.array(meta.names())
    data_columns = columns[columns != target_name]
    # TODO: stacking the content of the structured array
    # this results in a copy. If the data was homogeneous
    # we could use a view instead.
    X = np.column_stack(data[c] for c in data_columns)
    y = data[target_name]

    description = "{}\n\nDownloaded from openml.org.".format(
        data_description['description'])

    bunch = Bunch(
        data=X, target=y, feature_names=data_columns,
        DESCR=description)

    return bunch
