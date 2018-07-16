import json
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


import numpy as np

from .base import get_data_home
from sklearn.externals.arff import loads
from ..externals.joblib import Memory
from ..externals.six.moves.urllib.error import HTTPError
from ..utils import Bunch

_SEARCH_NAME = "https://openml.org/api/v1/json/data/list/data_name/{}/limit/1"
_DATA_INFO = "https://openml.org/api/v1/json/data/{}"
_DATA_FEATURES = "https://openml.org/api/v1/json/data/features/{}"


def _openml_fileid_url(file_id):
    return "https://www.openml.org/data/v1/download/{}/".format(file_id)


def _convert_arff_data(arff_data):
    """
    converts the arff object into the appropriate matrix type (now: np.ndarray,
    later also: scipy.sparse.csr_matrix) based on the 'data part' (i.e., in the
    liac-arff dict, the object from the 'data' key)

    Parameters
    ----------
    arff_data : list or dict
        as obtained from liac-arff object

    returns : np.ndarray (or later also: scipy.sparse.csr_matrix)
    """
    if isinstance(arff_data, list):
        X = np.array(arff_data, dtype=object)
    # elif: extendable for sparse arff in future ()
    else:
        raise ValueError('Unexpected Data Type obtained from arff ' +
                         '(This should never happen).')
    return X


def _get_data_info_by_name(name, version):
    """
    Utilizes the openml dataset listing api to find a dataset by
    name/version
    OpenML api fn:
    https://www.openml.org/api_docs#!/data/get_data_list_data_name_data_name

    Parameters
    ----------
    name : str
        name of the dataset

    version : int or str
        If version is an integer, the exact name/version will be obtained from
        OpenML. If version is a string (value: "active") it will take the first
        version from OpenML that is annotated as active. Any other string
        values except "active" are treated as integer.

    Returns
    -------
    json_data['data']['dataset'][0]: json
        json representation of the first dataset object that adhired to the
        search criteria

    """
    data_found = True
    try:
        if version == "active":
            response = urlopen(_SEARCH_NAME.format(name) + "/status/active/")
        else:
            response = urlopen((_SEARCH_NAME +
                                "/data_version/{}").format(name, version))
    except HTTPError as error:
        if error.code == 412:
            data_found = False
        else:
            raise error

    if not data_found and version != "active":
        # might have been deactivated. will warn later
        data_found = True
        try:
            url = (_SEARCH_NAME +
                   "/data_version/{}/status/deactivated").format(name, version)
            response = urlopen(url)
        except HTTPError as error:
            # 412 is an OpenML specific error code, indicating a generic error
            # (e.g., data not found)
            if error.code == 412:
                data_found = False
            else:
                raise error

    if not data_found:
        # not in except for nicer traceback
        if version == "active":
            raise ValueError("No active dataset {} found.".format(name))
        raise ValueError("Dataset {} with version {}"
                         " not found.".format(name, version))

    json_data = json.loads(response.read().decode("utf-8"))
    return json_data['data']['dataset'][0]


def _get_data_description_by_id(data_id):
    data_found = True
    try:
        response = urlopen(_DATA_INFO.format(data_id))
    except HTTPError as error:
        if error.code == 412:
            data_found = False
    if not data_found:
        # not in except for nicer traceback
        raise ValueError("Dataset with id {} "
                         "not found.".format(data_id))
    json_data = json.loads(response.read().decode("utf-8"))
    response.close()
    return json_data['data_set_description']


def _get_data_features(data_id):
    data_found = True
    try:
        response = urlopen(_DATA_FEATURES.format(data_id))
    except HTTPError as error:
        # 412 is an OpenML specific error code, indicating a generic error
        # (e.g., data not found)
        if error.code == 412:
            data_found = False
        else:
            raise error
    if not data_found:
        # not in except for nicer traceback
        raise ValueError("Dataset with id {} "
                         "not found.".format(data_id))
    json_data = json.loads(response.read().decode("utf-8"))
    response.close()
    return json_data['data_features']['feature']


def _download_data(file_id):
    url = _openml_fileid_url(file_id)
    response = urlopen(url)
    if sys.version_info[0] == 2:
        # Python2.7 numpy can't handle unicode?
        arff = loads(response.read())
    else:
        arff = loads(response.read().decode('utf-8'))

    response.close()
    return arff


def _convert_numericals(data, name_feature):
    for feature in name_feature.values():
        if feature['data_type'] == "numeric":
            idx = int(feature['index'])
            data[:, idx] = data[:, idx].astype(np.float)
    return data


def fetch_openml(id=None, name=None, version='active', data_home=None,
                 target_column_name='default-target', cache=True):
    """Fetch dataset from openml by name or dataset id.

    Datasets are uniquely identified by either an integer ID or by a
    combination of name and version (i.e. there might be multiple
    versions of the 'iris' dataset). Please give either name or id
    (not both). In case a name is given, a version can also be
    provided.

    Parameters
    ----------
    id : int
        OpenML ID of the dataset. The most specific way of retrieving a
        dataset. If ID is not given, name (and potential version) are
        used to obtain a dataset.

    name : string
        String identifier of the dataset. Note that OpenML can have multiple
        datasets with the same name.

    version : integer or 'active', default='active'
        Version of the dataset. Can only be provided if also ``name`` is given.
        If 'active' the oldest version that's still active is used.

    data_home : string or None, default None
        Specify another download and cache folder for the data sets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    target_column_name : string or None, default 'default-target'
        Specify the column name in the data to use as target. If
        'default-target', the standard target column a stored on the server
        is used. If ``None``, all columns are returned as data and the
        target is ``None``.

    cache : boolean, default=True
        Whether to cache downloaded datasets using joblib.

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
    if cache:
        mem = Memory(join(data_home, 'cache'), verbose=0).cache
    else:
        def mem(func):
            return func
    _get_data_info_by_name_ = mem(_get_data_info_by_name)
    _get_data_description_by_id_ = mem(_get_data_description_by_id)
    _get_data_features_ = mem(_get_data_features)
    _download_data_ = mem(_download_data)

    if not exists(data_home):
        os.makedirs(data_home)

    # check legal function arguments. id XOR (name, version) should be provided
    if name is not None:
        if id is not None:
            raise ValueError(
                "Dataset id={} and name={} passed, but you can only "
                "specify a numeric id or a name, not both.".format(id, name))
        data_info = _get_data_info_by_name_(name, version)
        data_id = data_info['did']
    elif id is not None:
        # from the previous if statement, it is given that name is None
        if version is not "active":
            raise ValueError(
                "Dataset id={} and version={} passed, but you can only "
                "specify a numeric id or a version, not both.".format(id,
                                                                      name))
        data_id = id
    else:
        raise ValueError(
            "Neither name nor id are provided. Please provide name xor id.")

    data_description = _get_data_description_by_id_(data_id)
    if data_description['status'] != "active":
        warn("Version {} of dataset {} is inactive, meaning that issues have"
             " been found in the dataset. Try using a newer version.".format(
                 data_description['version'], data_description['name']))
    if target_column_name == "default-target":
        target_column_name = data_description.get('default_target_attribute',
                                                  None)

    # download actual data
    features = _get_data_features_(data_id)
    name_feature = {feature['name']: feature for feature in features}

    # TODO: stacking the content of the structured array
    # this results in a copy. If the data was homogeneous
    # and target at start or end, we could use a view instead.
    data_columns = []
    for feature in features:
        if (feature['name'] != target_column_name and feature['is_ignore'] ==
                'false' and feature['is_row_identifier'] == 'false'):
            data_columns.append(feature['name'])

    arff_data = _download_data_(data_description['file_id'])['data']
    data = _convert_arff_data(arff_data)
    data = _convert_numericals(data, name_feature)

    if target_column_name is not None:
        y = data[:, int(name_feature[target_column_name]['index'])]
    else:
        y = None

    if all([feature['data_type'] == "numeric" for feature in features
            if feature['name'] in data_columns]):
        dtype = None
    else:
        dtype = object
    col_slice = [int(name_feature[col_name]['index'])
                 for col_name in data_columns]
    X = data[:, col_slice].astype(dtype)

    description = u"{}\n\nDownloaded from openml.org.".format(
        data_description.pop('description'))

    bunch = Bunch(
        data=X, target=y, feature_names=data_columns,
        DESCR=description, details=data_description, features=features,
        url="https://www.openml.org/d/{}".format(data_id))

    return bunch
