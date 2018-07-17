import json
import sys
import os
from os.path import join, exists
from warnings import warn

try:
    # Python 3+
    from urllib.request import urlopen
except ImportError:
    # Python 2
    from urllib2 import urlopen


import numpy as np
import sklearn.externals.arff as arff

from .base import get_data_home
from ..externals.joblib import Memory
from ..externals.six.moves.urllib.error import HTTPError
from ..utils import Bunch

__all__ = ['fetch_openml']


_SEARCH_NAME = "https://openml.org/api/v1/json/data/list/data_name/{}/limit/1"
_DATA_INFO = "https://openml.org/api/v1/json/data/{}"
_DATA_FEATURES = "https://openml.org/api/v1/json/data/features/{}"
_DATA_FILE = "https://openml.org/data/v1/download/{}"


def _get_json_content_from_openml_api(url, error_message, raise_if_error):
    """
    Loads json data from the openml api

    Parameters
    ----------
    url : str
        The URL to load from. Should be an official OpenML endpoint

    error_message : str
        The error message to raise if an acceptable OpenML error is thrown
        (acceptable error is, e.g., data id not found. Other errors, like 404's
        will throw the native error message)

    raise_if_error : bool
        Whether to raise an error if OpenML returns an acceptable error (e.g.,
        date not found). If this argument is set to False, a None is returned
        in case of acceptable errors. Note that all other errors (e.g., 404)
        will still be raised as normal.

    Returns
    -------
    json_data : json or None
        the json result from the OpenML server if the call was successful;
        None otherwise iff raise_if_error was set to False and the error was
        ``acceptable``
    """
    data_found = True
    try:
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
        if raise_if_error:
            raise ValueError(error_message)
        else:
            return None
    json_data = json.loads(response.read().decode("utf-8"))
    response.close()
    return json_data


def _convert_arff_data(arff_data):
    """
    converts the arff object into the appropriate matrix type (now: np.ndarray,
    later also: scipy.sparse.csr_matrix) based on the 'data part' (i.e., in the
    liac-arff dict, the object from the 'data' key)

    Parameters
    ----------
    arff_data : list or dict
        as obtained from liac-arff object

    Returns
    -------
    X : np.ndarray
        (or later also: scipy.sparse.csr_matrix)
    """
    if isinstance(arff_data, list):
        X = np.array(arff_data, dtype=object)
    # elif: extendable for sparse arff in future ()
    else:
        # This should never happen
        raise ValueError('Unexpected Data Type obtained from arff.')
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
    if version == "active":
        # situation in which we return the oldest active version
        url = _SEARCH_NAME.format(name) + "/status/active/"
        error_msg = "No active dataset {} found.".format(name)
        json_data = _get_json_content_from_openml_api(url, error_msg, True)
        return json_data['data']['dataset'][0]

    # an integer version has been provided
    url = (_SEARCH_NAME + "/data_version/{}").format(name, version)
    json_data = _get_json_content_from_openml_api(url, None, False)
    if json_data is None:
        # we can do this in 1 fn call if OpenML does not require the
        # specification of the dataset status (i.e., return datasets with a
        # given name / version regardless of active, deactivated, etc. )
        # TODO: feature request OpenML.
        url += "/status/deactivated"
        error_msg = "Dataset {} with version {} not found.".format(name,
                                                                   version)
        json_data = _get_json_content_from_openml_api(url, error_msg, True)

    return json_data['data']['dataset'][0]


def _get_data_description_by_id(data_id):
    # OpenML API fn: https://www.openml.org/api_docs#!/data/get_data_id
    url = _DATA_INFO.format(data_id)
    error_message = "Dataset with id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)
    return json_data['data_set_description']


def _get_data_features(data_id):
    # OpenML fn: https://www.openml.org/api_docs#!/data/get_data_features_id
    url = _DATA_FEATURES.format(data_id)
    error_message = "Dataset with id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)
    return json_data['data_features']['feature']


def _download_data_arff(file_id):
    # Accesses an ARFF file on the OpenML server. Documentation:
    # https://www.openml.org/api_data_docs#!/data/get_download_id
    url = _DATA_FILE.format(file_id)
    response = urlopen(url)
    if sys.version_info[0] == 2:
        # Python2.7 numpy can't handle unicode?
        arff_file = arff.loads(response.read())
    else:
        arff_file = arff.loads(response.read().decode('utf-8'))

    response.close()
    return arff_file


def _convert_numericals(data, name_feature):
    # converts all numerical columns into floats
    for feature in name_feature.values():
        if feature['data_type'] == "numeric":
            idx = int(feature['index'])
            data[:, idx] = data[:, idx].astype(np.float)
    return data


def fetch_openml(name=None, version='active', data_id=None, data_home=None,
                 target_column_name='default-target', cache=True):
    """Fetch dataset from openml by name or dataset id.

    Datasets are uniquely identified by either an integer ID or by a
    combination of name and version (i.e. there might be multiple
    versions of the 'iris' dataset). Please give either name or id
    (not both). In case a name is given, a version can also be
    provided.

    Parameters
    ----------
    name : str or None
        String identifier of the dataset. Note that OpenML can have multiple
        datasets with the same name.

    version : integer or 'active', default='active'
        Version of the dataset. Can only be provided if also ``name`` is given.
        If 'active' the oldest version that's still active is used.

    data_id : int or None
        OpenML ID of the dataset. The most specific way of retrieving a
        dataset. If data_id is not given, name (and potential version) are
        used to obtain a dataset.

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
    cached_get_data_info_by_name = mem(_get_data_info_by_name)
    cached_get_data_description_by_id = mem(_get_data_description_by_id)
    cached_get_data_features = mem(_get_data_features)
    cached_download_data_arff = mem(_download_data_arff)

    if not exists(data_home):
        os.makedirs(data_home)

    # check legal function arguments. data_id XOR (name, version) should be
    # provided
    if name is not None:
        # OpenML is case-insensitive, but the caching mechanism is not
        # convert all data names (str) to lower case
        name = name.lower()
        if data_id is not None:
            raise ValueError(
                "Dataset data_id={} and name={} passed, but you can only "
                "specify a numeric data_id or a name, not "
                "both.".format(data_id, name))
        data_info = cached_get_data_info_by_name(name, version)
        data_id = data_info['did']
    elif data_id is not None:
        # from the previous if statement, it is given that name is None
        if version is not "active":
            raise ValueError(
                "Dataset data_id={} and version={} passed, but you can only "
                "specify a numeric data_id or a version, not "
                "both.".format(data_id, name))
    else:
        raise ValueError(
            "Neither name nor data_id are provided. Please provide name or "
            "data_id.")

    data_description = cached_get_data_description_by_id(data_id)
    if data_description['status'] != "active":
        warn("Version {} of dataset {} is inactive, meaning that issues have"
             " been found in the dataset. Try using a newer version.".format(
                 data_description['version'], data_description['name']))
    if target_column_name == "default-target":
        target_column_name = data_description.get('default_target_attribute',
                                                  None)

    # download actual data
    features = cached_get_data_features(data_id)
    name_feature = {feature['name']: feature for feature in features}

    # TODO: stacking the content of the structured array
    # this results in a copy. If the data was homogeneous
    # and target at start or end, we could use a view instead.
    data_columns = []
    for feature in features:
        if (feature['name'] != target_column_name and feature['is_ignore'] ==
                'false' and feature['is_row_identifier'] == 'false'):
            data_columns.append(feature['name'])

    arff_data = cached_download_data_arff(data_description['file_id'])['data']
    data = _convert_arff_data(arff_data)
    data = _convert_numericals(data, name_feature)

    if target_column_name is not None:
        # determine vector type
        if name_feature[target_column_name]['data_type'] == "numeric":
            dtype = np.float64
        else:
            dtype = object
        y = np.array(data[:, int(name_feature[target_column_name]['index'])],
                     dtype=dtype)
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
