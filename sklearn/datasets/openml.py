import json
import numbers
import sys
import os
from os.path import join, exists
from warnings import warn

try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen
    # from http.client import IncompleteRead

from scipy.io.arff import loadarff
import numpy as np

from .base import get_data_home
from ..externals.joblib import Memory
from ..externals.six import StringIO
from ..utils import Bunch

_SEARCH_NAME = "https://openml.org/api/v1/json/data/list/data_name/{}/limit/1"
_DATA_INFO = "https://openml.org/api/v1/json/data/{}"


def _get_data_info_by_name(name, version):
    if version == "active":
        json_string = urlopen(_SEARCH_NAME.format(name + "/status/active/"))
    else:
        # FIXME waiting for new filter mechanism
        json_string = urlopen(_SEARCH_NAME.format(name))
    json_data = json.load(json_string)
    return json_data['data']['dataset'][0]


def _get_data_description_by_id(data_id):
    json_string = urlopen(_DATA_INFO.format(data_id))
    json_data = json.load(json_string)
    return json_data['data_set_description']


def _download_data(url):
    response = urlopen(url)
    if sys.version_info[0] == 2:
        # Python2.7 numpy can't handle unicode?
        arff = loadarff(StringIO(response.read()))
    else:
        arff = loadarff(StringIO(response.read().decode('utf-8')))

    response.close()
    return arff


def fetch_openml(name_or_id=None, version='active', data_home=None,
                 target_column='default-target', memory=True):
    """Fetch dataset from openml by name or dataset id.

    Datasets are uniquely identified by either an integer ID or by a
    combination of name and version (i.e. there might be multiple
    versions of the 'iris' dataset).

    Parameters
    ----------
    name_or_id : string or integer
        Identifier of the dataset. If integer, assumed to be the id of the
        dataset on OpenML, if string, assumed to be the name of the dataset.

    version : integer or 'active', default='active'
        Version of the dataset. Only used if ``name_or_id`` is a string.
        If 'active' the oldest version that's still active is used.

    data_home : string or None, default None
        Specify another download and cache folder for the data sets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    target_column : string or None, default 'default-target'
        Specify the column name in the data to use as target. If
        'default-target', the standard target column a stored on the server
        is used. If ``None``, all columns are returned as data and the
        tharget is ``None``.

    memory : boolean, default=True
        Whether to store downloaded datasets using joblib.

    Returns
    -------

    data : Bunch
        Dictionary-like object, the interesting attributes are:
        'data', the data to learn, 'target', the regression target or
        classification labels, 'DESCR', the full description of the dataset,
        'feature_names', the original names of the dataset columns, and
        'details' which provide more information on the openml meta-data.
    """
    data_home = get_data_home(data_home=data_home)
    data_home = join(data_home, 'openml')
    if memory:
        mem = Memory(join(data_home, 'cache'), verbose=0).cache
    else:
        def mem(func):
            return func
    _get_data_info_by_name_ = mem(_get_data_info_by_name)
    _get_data_description_by_id_ = mem(_get_data_description_by_id)
    _download_data_ = mem(_download_data)

    if not exists(data_home):
        os.makedirs(data_home)

    # check if dataset id is known
    if isinstance(name_or_id, numbers.Integral):
        if version != "active":
            raise ValueError(
                "Dataset id={} and version={} passed, but you can only "
                "specify a numeric id or a version, not both.".format(
                    name_or_id, version))
        data_id = name_or_id
    elif isinstance(name_or_id, str):
        data_info = _get_data_info_by_name_(name_or_id, version)
        data_id = data_info['did']

    else:
        raise TypeError(
            "Invalid name_or_id {}, should be string or integer.".format(
                name_or_id))

    data_description = _get_data_description_by_id_(data_id)
    if data_description['status'] != "active":
        warn("Version {} of dataset {} is inactive, meaning that issues have"
             " been found in the dataset. Try using a newer version.".format(
                 data_description['name'], data_description['version']))
    if target_column == "default-target":
        target_column = data_description.get('default_target_attribute', None)

    # download actual data
    data, meta = _download_data_(data_description['url'])
    columns = np.array(meta.names())
    data_columns = columns[columns != target_column]
    # TODO: stacking the content of the structured array
    # this results in a copy. If the data was homogeneous
    # we could use a view instead.
    dtype = object if "nominal" in meta.types() else None
    X = np.array([data[c] for c in data_columns], dtype=dtype).T
    if target_column is not None:
        y = data[target_column]
    else:
        y = None

    description = u"{}\n\nDownloaded from openml.org.".format(
        data_description.pop('description'))

    bunch = Bunch(
        data=X, target=y, feature_names=data_columns,
        DESCR=description, details=data_description, meta=meta)

    return bunch
