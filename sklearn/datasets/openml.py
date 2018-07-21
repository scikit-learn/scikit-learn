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

from sklearn.externals import arff
from .base import get_data_home
from ..externals.six import string_types, PY2
from ..externals.six.moves.urllib.error import HTTPError
from ..utils import Bunch, Memory

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

    error_message : str or None
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
    X : np.array
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
    OpenML api function:
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
        # we can do this in 1 function call if OpenML does not require the
        # specification of the dataset status (i.e., return datasets with a
        # given name / version regardless of active, deactivated, etc. )
        # TODO: feature request OpenML.
        url += "/status/deactivated"
        error_msg = "Dataset {} with version {} not found.".format(name,
                                                                   version)
        json_data = _get_json_content_from_openml_api(url, error_msg, True)

    return json_data['data']['dataset'][0]


def _get_data_description_by_id(data_id):
    # OpenML API function: https://www.openml.org/api_docs#!/data/get_data_id
    url = _DATA_INFO.format(data_id)
    error_message = "Dataset with id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)
    return json_data['data_set_description']


def _get_data_features(data_id):
    # OpenML function:
    # https://www.openml.org/api_docs#!/data/get_data_features_id
    url = _DATA_FEATURES.format(data_id)
    error_message = "Dataset with id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)
    return json_data['data_features']['feature']


def _download_data_arff(file_id):
    # Accesses an ARFF file on the OpenML server. Documentation:
    # https://www.openml.org/api_data_docs#!/data/get_download_id
    url = _DATA_FILE.format(file_id)
    response = urlopen(url)
    if PY2:
        # Because Python2.7 numpy can't handle unicode
        arff_file = arff.loads(response.read())
    else:
        arff_file = arff.loads(response.read().decode('utf-8'))

    response.close()
    return arff_file


def _convert_numericals(data, name_feature):
    # converts all numerical columns in the numpy array into floats
    for feature in name_feature.values():
        if feature['data_type'] == "numeric":
            idx = int(feature['index'])
            data[:, idx] = data[:, idx].astype(np.float)
    return data


def _determine_default_target(name_feature):
    # determines the default target based on the data feature results
    # (which is currently more reliable than the data description;
    # see issue: https://github.com/openml/OpenML/issues/768)
    results = []
    for name, feature in name_feature.items():
        # note: string comparison (not boolean)
        if feature['is_target'] == "true":
            results.append(name)

    if len(results) == 0:
        return None
    elif len(results) == 1:
        return results[0]
    else:
        return results


def _determine_single_target_data_type(name_feature, target_column_name):
    # determine the data type of the y array in case there is a single target
    if not isinstance(target_column_name, string_types):
        raise ValueError('target_column_name should be of string type, '
                         'got: %s' % type(target_column_name))
    if target_column_name not in name_feature:
        raise KeyError('Could not find target_column_name={}')
    feature = name_feature[target_column_name]
    # note: we compare to a string, not boolean
    if feature['is_ignore'] == 'true':
        warn('target_column_name={} has flag is_ignore.'.format(
            target_column_name))
    if feature['is_row_identifier'] == 'true':
        warn('target_column_name={} has flag is_row_identifier.'.format(
            target_column_name))

    if feature['data_type'] == "numeric":
        return np.float64
    else:
        return object


def _determine_multi_target_data_type(name_feature, target_column_names):
    # determine the data type of the y array in case there are multiple targets
    # (throws an error if these targets do not comply with sklearn support)
    if not isinstance(target_column_names, list):
        raise ValueError('target_column_name should be list, '
                         'got: %s' % type(target_column_names))
    found_types = set()
    for target_column_name in target_column_names:
        if target_column_name not in name_feature:
            raise KeyError('Could not find target_column_name={}')
        if name_feature[target_column_name]['data_type'] == "numeric":
            found_types.add(np.float64)
        else:
            found_types.add(object)

        # note: we compare to a string, not boolean
        if name_feature[target_column_name]['is_ignore'] == 'true':
            warn('target_column_name={} has flag is_ignore.'.format(
                target_column_name))
        if name_feature[target_column_name]['is_row_identifier'] == 'true':
            warn('target_column_name={} has flag is_row_identifier.'.format(
                target_column_name))
    if len(found_types) != 1:
        raise ValueError('Can only handle homogeneous multi-target datasets, '
                         'i.e., all targets are either numeric or '
                         'categorical.')
    return list(found_types)[0]


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

    target_column_name : string, list or None, default 'default-target'
        Specify the column name in the data to use as target. If
        'default-target', the standard target column a stored on the server
        is used. If ``None``, all columns are returned as data and the
        target is ``None``. If list (of strings), all columns with these names
        are returned as multi-target (Note: not all scikit-learn classifiers
        can handle all types of multi-output combinations)

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

    # check valid function arguments. data_id XOR (name, version) should be
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

    # download data features, meta-info about column types
    features = cached_get_data_features(data_id)
    name_feature = {feature['name']: feature for feature in features}

    if target_column_name == "default-target":
        target_column_name = _determine_default_target(name_feature)

    # TODO: stacking the content of the structured array
    # this results in a copy. If the data was homogeneous
    # and target at start or end, we could use a view instead.
    data_columns = []
    for feature in features:
        # determine whether `feature` is a target
        if (isinstance(target_column_name, string_types) and
                feature['name'] == target_column_name):
            is_target = True
        elif (isinstance(target_column_name, list) and
                feature['name'] in target_column_name):
            is_target = True
        else:
            is_target = False

        if ((not is_target) and feature['is_ignore'] == 'false' and
                feature['is_row_identifier'] == 'false'):
            data_columns.append(feature['name'])

    arff_data = cached_download_data_arff(data_description['file_id'])['data']
    data = _convert_arff_data(arff_data)
    data = _convert_numericals(data, name_feature)

    if isinstance(target_column_name, string_types):
        # determine vector type
        dtype = _determine_single_target_data_type(name_feature,
                                                   target_column_name)
        y = np.array(data[:, int(name_feature[target_column_name]['index'])],
                     dtype=dtype)
    elif isinstance(target_column_name, list):
        dtype = _determine_multi_target_data_type(name_feature,
                                                  target_column_name)
        indices = [int(name_feature[col_name]['index'])
                   for col_name in target_column_name]
        y = np.array(data[:, indices], dtype=dtype)
    elif target_column_name is None:
        y = None
    else:
        # unexpected behaviour, this should never happen
        raise ValueError('Could not determine how to handle '
                         'target_column_name of type {}'.
                         format(type(target_column_name)))

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
