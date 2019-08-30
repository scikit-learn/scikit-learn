import gzip
import json
from os.path import join
from warnings import warn
from contextlib import closing
import itertools
from collections.abc import Generator
from collections import OrderedDict

from urllib.request import urlopen, Request

import numpy as np
import scipy.sparse
from joblib import Memory

from ..externals import _arff
from .base import get_data_home
from urllib.error import HTTPError
from ..utils import Bunch
from ..utils import get_chunk_n_rows
from ..utils import _chunk_generator
from ..utils import check_pandas_support  # noqa

__all__ = ['fetch_openml']

_OPENML_PREFIX = "https://openml.org/"
_SEARCH_NAME = "api/v1/json/data/list/data_name/{}/limit/2"
_DATA_INFO = "api/v1/json/data/{}"
_DATA_FEATURES = "api/v1/json/data/features/{}"
_DATA_QUALITIES = "api/v1/json/data/qualities/{}"
_DATA_FILE = "data/v1/download/{}"


class _CacheWithRetry:
    """Apply a function with a joblib.Memory cache

    If first call fails, likely due to pickling issue, clean the
    cache and try again.
    """
    def __init__(self, memory, data_home):
        self.memory = memory
        self.data_home = data_home

    def __call__(self, func, **kwargs):
        # TODO: the retry could be handled in joblib (maybe it already is?)
        memory = self.memory
        if memory is not None and self.data_home is not None:
            try:
                return memory.cache(func)(**kwargs)
            except HTTPError:
                raise
            except Exception:
                warn("Invalid cache, redownloading file", RuntimeWarning)

            # attemting to download data with cleared cache
            memory.clear()
            return memory.cache(func)(**kwargs)
        else:
            return func(**kwargs)


def _open_openml_url(openml_path):
    """
    Returns a resource from OpenML.org. Caches it to data_home if required.

    Parameters
    ----------
    openml_path : str
        OpenML URL that will be accessed. This will be prefixes with
        _OPENML_PREFIX

    Returns
    -------
    result : stream
        A stream to the OpenML resource
    """
    def is_gzip(_fsrc):
        return _fsrc.info().get('Content-Encoding', '') == 'gzip'

    req = Request(_OPENML_PREFIX + openml_path)
    req.add_header('Accept-encoding', 'gzip')

    fsrc = urlopen(req)
    if is_gzip(fsrc):
        return gzip.GzipFile(fileobj=fsrc, mode='rb')
    return fsrc


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

    def _load_json():
        with closing(_open_openml_url(url)) as response:
            return json.loads(response.read().decode("utf-8"))

    try:
        return _load_json()
    except HTTPError as error:
        # 412 is an OpenML specific error code, indicating a generic error
        # (e.g., data not found)
        if error.code != 412:
            raise error

    # 412 error, not in except for nicer traceback
    if raise_if_error:
        raise ValueError(error_message)
    return None


def _split_sparse_columns(arff_data, include_columns):
    """
    obtains several columns from sparse arff representation. Additionally, the
    column indices are re-labelled, given the columns that are not included.
    (e.g., when including [1, 2, 3], the columns will be relabelled to
    [0, 1, 2])

    Parameters
    ----------
    arff_data : tuple
        A tuple of three lists of equal size; first list indicating the value,
        second the x coordinate and the third the y coordinate.

    include_columns : list
        A list of columns to include.

    Returns
    -------
    arff_data_new : tuple
        Subset of arff data with only the include columns indicated by the
        include_columns argument.
    """
    arff_data_new = (list(), list(), list())
    reindexed_columns = {column_idx: array_idx for array_idx, column_idx
                         in enumerate(include_columns)}
    for val, row_idx, col_idx in zip(arff_data[0], arff_data[1], arff_data[2]):
        if col_idx in include_columns:
            arff_data_new[0].append(val)
            arff_data_new[1].append(row_idx)
            arff_data_new[2].append(reindexed_columns[col_idx])
    return arff_data_new


def _sparse_data_to_array(arff_data, include_columns):
    # turns the sparse data back into an array (can't use toarray() function,
    # as this does only work on numeric data)
    num_obs = max(arff_data[1]) + 1
    y_shape = (num_obs, len(include_columns))
    reindexed_columns = {column_idx: array_idx for array_idx, column_idx
                         in enumerate(include_columns)}
    # TODO: improve for efficiency
    y = np.empty(y_shape, dtype=np.float64)
    for val, row_idx, col_idx in zip(arff_data[0], arff_data[1], arff_data[2]):
        if col_idx in include_columns:
            y[row_idx, reindexed_columns[col_idx]] = val
    return y


def _convert_arff_data(arff_data, col_slice_x, col_slice_y, shape=None):
    """
    converts the arff object into the appropriate matrix type (np.array or
    scipy.sparse.csr_matrix) based on the 'data part' (i.e., in the
    liac-arff dict, the object from the 'data' key)

    Parameters
    ----------
    arff_data : list or dict
        as obtained from liac-arff object

    col_slice_x : list
        The column indices that are sliced from the original array to return
        as X data

    col_slice_y : list
        The column indices that are sliced from the original array to return
        as y data

    Returns
    -------
    X : np.array or scipy.sparse.csr_matrix
    y : np.array
    """
    if isinstance(arff_data, Generator):
        if shape[0] == -1:
            count = -1
        else:
            count = shape[0] * shape[1]
        data = np.fromiter(itertools.chain.from_iterable(arff_data),
                           dtype='float64', count=count)
        data = data.reshape(*shape)
        X = data[:, col_slice_x]
        y = data[:, col_slice_y]
        return X, y
    elif isinstance(arff_data, tuple):
        arff_data_X = _split_sparse_columns(arff_data, col_slice_x)
        num_obs = max(arff_data[1]) + 1
        X_shape = (num_obs, len(col_slice_x))
        X = scipy.sparse.coo_matrix(
            (arff_data_X[0], (arff_data_X[1], arff_data_X[2])),
            shape=X_shape, dtype=np.float64)
        X = X.tocsr()
        y = _sparse_data_to_array(arff_data, col_slice_y)
        return X, y
    else:
        # This should never happen
        raise ValueError('Unexpected Data Type obtained from arff.')


def _feature_to_dtype(feature):
    """Map feature to dtype for pandas DataFrame
    """
    if feature['data_type'] == 'string':
        return object
    elif feature['data_type'] == 'nominal':
        return 'category'
    # only numeric, integer, real are left
    elif (feature['number_of_missing_values'] != '0' or
          feature['data_type'] in ['numeric', 'real']):
        # cast to floats when there are any missing values
        return np.float64
    elif feature['data_type'] == 'integer':
        return np.int64
    raise ValueError('Unsupported feature: {}'.format(feature))


def _convert_arff_data_dataframe(arff, columns, features_dict):
    """Convert the ARFF object into a pandas DataFrame.

    Parameters
    ----------
    arff : dict
        As obtained from liac-arff object.

    columns : list
        Columns from dataframe to return.

    features_dict : dict
        Maps feature name to feature info from openml.

    Returns
    -------
    dataframe : pandas DataFrame
    """
    pd = check_pandas_support('fetch_openml with as_frame=True')

    attributes = OrderedDict(arff['attributes'])
    arff_columns = list(attributes)

    # calculate chunksize
    first_row = next(arff['data'])
    first_df = pd.DataFrame([first_row], columns=arff_columns)

    row_bytes = first_df.memory_usage(deep=True).sum()
    chunksize = get_chunk_n_rows(row_bytes)

    # read arff data with chunks
    columns_to_keep = [col for col in arff_columns if col in columns]
    dfs = []
    dfs.append(first_df[columns_to_keep])
    for data in _chunk_generator(arff['data'], chunksize):
        dfs.append(pd.DataFrame(data, columns=arff_columns)[columns_to_keep])
    df = pd.concat(dfs, ignore_index=True)

    for column in columns_to_keep:
        dtype = _feature_to_dtype(features_dict[column])
        if dtype == 'category':
            dtype = pd.api.types.CategoricalDtype(attributes[column])
        df[column] = df[column].astype(dtype, copy=False)
    return df


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
    first_dataset : json
        json representation of the first dataset object that adhired to the
        search criteria

    """
    if version == "active":
        # situation in which we return the oldest active version
        url = _SEARCH_NAME.format(name) + "/status/active/"
        error_msg = "No active dataset {} found.".format(name)
        json_data = _get_json_content_from_openml_api(url, error_msg, True)
        res = json_data['data']['dataset']
        if len(res) > 1:
            warn("Multiple active versions of the dataset matching the name"
                 " {name} exist. Versions may be fundamentally different, "
                 "returning version"
                 " {version}.".format(name=name, version=res[0]['version']))
        return res[0]

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
    error_message = "Dataset with data_id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)
    return json_data['data_set_description']


def _get_data_features(data_id):
    # OpenML function:
    # https://www.openml.org/api_docs#!/data/get_data_features_id
    url = _DATA_FEATURES.format(data_id)
    error_message = "Dataset with data_id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)
    return json_data['data_features']['feature']


def _get_data_qualities(data_id):
    # OpenML API function:
    # https://www.openml.org/api_docs#!/data/get_data_qualities_id
    url = _DATA_QUALITIES.format(data_id)
    error_message = "Dataset with data_id {} not found.".format(data_id)
    json_data = _get_json_content_from_openml_api(url, error_message, True)

    # the qualities might not be available, but we still try to process
    # the data
    return json_data.get('data_qualities', {}).get('quality', None)


def _get_num_samples(data_qualities):
    """Get the number of samples from data qualities.

    Parameters
    ----------
    data_qualities : list of dict
        Used to retrieve the number of instances (samples) in the dataset.

    Returns
    -------
    n_samples : int
        The number of samples in the dataset or -1 if data qualities are
        unavailable.
    """
    # If the data qualities are unavailable, we return -1
    default_n_samples = -1

    if data_qualities is None:
        return default_n_samples

    qualities = {d['name']: d['value'] for d in data_qualities}
    return int(float(qualities.get('NumberOfInstances', default_n_samples)))


def _download_data_arff(file_id, sparse, encode_nominal=True):
    # Accesses an ARFF file on the OpenML server. Documentation:
    # https://www.openml.org/api_data_docs#!/data/get_download_id
    # encode_nominal argument is to ensure unit testing, do not alter in
    # production!

    url = _DATA_FILE.format(file_id)

    with closing(_open_openml_url(url)) as response:
        if sparse is True:
            return_type = _arff.COO
        else:
            return_type = _arff.DENSE_GEN

        arff_file = _arff.loads(response.read().decode('utf-8'),
                                encode_nominal=encode_nominal,
                                return_type=return_type)
    return arff_file


def _verify_target_data_type(features_dict, target_columns):
    # verifies the data type of the y array in case there are multiple targets
    # (throws an error if these targets do not comply with sklearn support)
    if not isinstance(target_columns, list):
        raise ValueError('target_column should be list, '
                         'got: %s' % type(target_columns))
    found_types = set()
    for target_column in target_columns:
        if target_column not in features_dict:
            raise KeyError('Could not find target_column={}')
        if features_dict[target_column]['data_type'] == "numeric":
            found_types.add(np.float64)
        else:
            found_types.add(object)

        # note: we compare to a string, not boolean
        if features_dict[target_column]['is_ignore'] == 'true':
            warn('target_column={} has flag is_ignore.'.format(
                target_column))
        if features_dict[target_column]['is_row_identifier'] == 'true':
            warn('target_column={} has flag is_row_identifier.'.format(
                target_column))
    if len(found_types) > 1:
        raise ValueError('Can only handle homogeneous multi-target datasets, '
                         'i.e., all targets are either numeric or '
                         'categorical.')


def _valid_data_column_names(features_list, target_columns):
    # logic for determining on which columns can be learned. Note that from the
    # OpenML guide follows that columns that have the `is_row_identifier` or
    # `is_ignore` flag, these can not be learned on. Also target columns are
    # excluded.
    valid_data_column_names = []
    for feature in features_list:
        if (feature['name'] not in target_columns
                and feature['is_ignore'] != 'true'
                and feature['is_row_identifier'] != 'true'):
            valid_data_column_names.append(feature['name'])
    return valid_data_column_names


def fetch_openml(name=None, version='active', data_id=None, data_home=None,
                 target_column='default-target', cache=True, return_X_y=False,
                 as_frame=False):
    """Fetch dataset from openml by name or dataset id.

    Datasets are uniquely identified by either an integer ID or by a
    combination of name and version (i.e. there might be multiple
    versions of the 'iris' dataset). Please give either name or data_id
    (not both). In case a name is given, a version can also be
    provided.

    Read more in the :ref:`User Guide <openml>`.

    .. note:: EXPERIMENTAL

        The API is experimental (particularly the return value structure),
        and might have small backward-incompatible changes in future releases.

    Parameters
    ----------
    name : str or None
        String identifier of the dataset. Note that OpenML can have multiple
        datasets with the same name.

    version : integer or 'active', default='active'
        Version of the dataset. Can only be provided if also ``name`` is given.
        If 'active' the oldest version that's still active is used. Since
        there may be more than one active version of a dataset, and those
        versions may fundamentally be different from one another, setting an
        exact version is highly recommended.

    data_id : int or None
        OpenML ID of the dataset. The most specific way of retrieving a
        dataset. If data_id is not given, name (and potential version) are
        used to obtain a dataset.

    data_home : string or None, default None
        Specify another download and cache folder for the data sets. By default
        all scikit-learn data is stored in '~/scikit_learn_data' subfolders.

    target_column : string, list or None, default 'default-target'
        Specify the column name in the data to use as target. If
        'default-target', the standard target column a stored on the server
        is used. If ``None``, all columns are returned as data and the
        target is ``None``. If list (of strings), all columns with these names
        are returned as multi-target (Note: not all scikit-learn classifiers
        can handle all types of multi-output combinations)

    cache : boolean, default=True
        Whether to cache downloaded datasets using joblib.

    return_X_y : boolean, default=False.
        If True, returns ``(data, target)`` instead of a Bunch object. See
        below for more information about the `data` and `target` objects.

    as_frame : boolean, default=False
        If True, the data is a pandas DataFrame including columns with
        appropriate dtypes (numeric, string or categorical). The target is
        a pandas DataFrame or Series depending on the number of target_columns.
        The Bunch will contain a ``frame`` attribute with the target and the
        data. If ``return_X_y`` is True, then ``(data, target)`` will be pandas
        DataFrames or Series as describe above.

    Returns
    -------

    data : Bunch
        Dictionary-like object, with attributes:

        data : np.array, scipy.sparse.csr_matrix of floats, or pandas DataFrame
            The feature matrix. Categorical features are encoded as ordinals.
        target : np.array, pandas Series or DataFrame
            The regression target or classification labels, if applicable.
            Dtype is float if numeric, and object if categorical. If
            ``as_frame`` is True, ``target`` is a pandas object.
        DESCR : str
            The full description of the dataset
        feature_names : list
            The names of the dataset columns
        categories : dict or None
            Maps each categorical feature name to a list of values, such
            that the value encoded as i is ith in the list. If ``as_frame``
            is True, this is None.
        details : dict
            More metadata from OpenML
        frame : pandas DataFrame
            Only present when `as_frame=True`. DataFrame with ``data`` and
            ``target``.

    (data, target) : tuple if ``return_X_y`` is True

        .. note:: EXPERIMENTAL

            This interface is **experimental** and subsequent releases may
            change attributes without notice (although there should only be
            minor changes to ``data`` and ``target``).

        Missing values in the 'data' are represented as NaN's. Missing values
        in 'target' are represented as NaN's (numerical target) or None
        (categorical target)
    """
    data_home = get_data_home(data_home=data_home)
    data_home = join(data_home, 'openml')

    if cache and data_home is not None:
        # Use dataset specific cache location, so it can be separately cleaned
        memory = Memory(location=join(data_home, 'cache', str(data_id)),
                        verbose=0)
    else:
        memory = None
    _apply_cached_with_retry = _CacheWithRetry(memory, data_home)

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
        data_info = _apply_cached_with_retry(
            _get_data_info_by_name, name=name, version=version)
        data_id = data_info['did']
    elif data_id is not None:
        # from the previous if statement, it is given that name is None
        if version != "active":
            raise ValueError(
                "Dataset data_id={} and version={} passed, but you can only "
                "specify a numeric data_id or a version, not "
                "both.".format(data_id, name))
    else:
        raise ValueError(
            "Neither name nor data_id are provided. Please provide name or "
            "data_id.")

    bunch = _apply_cached_with_retry(
        _fetch_openml_inner, data_id=data_id, target_column=target_column,
        as_frame=as_frame)

    if return_X_y:
        return bunch.data, bunch.target
    else:
        return bunch


def _fetch_openml_inner(data_id=None, target_column='default-target',
                        as_frame=False):

    data_description = _get_data_description_by_id(data_id)
    if data_description['status'] != "active":
        warn("Version {} of dataset {} is inactive, meaning that issues have "
             "been found in the dataset. Try using a newer version from "
             "this URL: {}".format(
                data_description['version'],
                data_description['name'],
                data_description['url']))
    if 'error' in data_description:
        warn("OpenML registered a problem with the dataset. It might be "
             "unusable. Error: {}".format(data_description['error']))
    if 'warning' in data_description:
        warn("OpenML raised a warning on the dataset. It might be "
             "unusable. Warning: {}".format(data_description['warning']))

    return_sparse = False
    if data_description['format'].lower() == 'sparse_arff':
        return_sparse = True

    if as_frame and return_sparse:
        raise ValueError('Cannot return dataframe with sparse data')

    # download data features, meta-info about column types
    features_list = _get_data_features(data_id)

    if not as_frame:
        for feature in features_list:
            if 'true' in (feature['is_ignore'], feature['is_row_identifier']):
                continue
            if feature['data_type'] == 'string':
                raise ValueError('STRING attributes are not supported for '
                                 'array representation. Try as_frame=True')

    if target_column == "default-target":
        # determines the default target based on the data feature results
        # (which is currently more reliable than the data description;
        # see issue: https://github.com/openml/OpenML/issues/768)
        target_columns = [feature['name'] for feature in features_list
                          if feature['is_target'] == 'true']
    elif isinstance(target_column, str):
        # for code-simplicity, make target_column by default a list
        target_columns = [target_column]
    elif target_column is None:
        target_columns = []
    elif isinstance(target_column, list):
        target_columns = target_column
    else:
        raise TypeError("Did not recognize type of target_column"
                        "Should be str, list or None. Got: "
                        "{}".format(type(target_column)))
    data_columns = _valid_data_column_names(features_list,
                                            target_columns)

    # prepare which columns and data types should be returned for the X and y
    features_dict = {feature['name']: feature for feature in features_list}

    # XXX: col_slice_y should be all nominal or all numeric
    _verify_target_data_type(features_dict, target_columns)

    col_slice_y = [int(features_dict[col_name]['index'])
                   for col_name in target_columns]

    col_slice_x = [int(features_dict[col_name]['index'])
                   for col_name in data_columns]
    for col_idx in col_slice_y:
        feat = features_list[col_idx]
        nr_missing = int(feat['number_of_missing_values'])
        if nr_missing > 0:
            raise ValueError('Target column {} has {} missing values. '
                             'Missing values are not supported for target '
                             'columns. '.format(feat['name'], nr_missing))

    # determine arff encoding to return
    if not return_sparse:
        # The shape must include the ignored features to keep the right indexes
        # during the arff data conversion.
        data_qualities = _get_data_qualities(data_id)
        shape = _get_num_samples(data_qualities), len(features_list)
    else:
        shape = None

    # obtain the data
    arff = _download_data_arff(data_description['file_id'], return_sparse,
                               encode_nominal=not as_frame)

    description = "{}\n\nDownloaded from openml.org.".format(
        data_description.pop('description'))

    nominal_attributes = None
    frame = None
    if as_frame:
        columns = data_columns + target_columns
        frame = _convert_arff_data_dataframe(arff, columns, features_dict)
        X = frame[data_columns]
        if len(target_columns) >= 2:
            y = frame[target_columns]
        elif len(target_columns) == 1:
            y = frame[target_columns[0]]
        else:
            y = None
    else:
        # nominal attributes is a dict mapping from the attribute name to the
        # possible values. Includes also the target column (which will be
        # popped off below, before it will be packed in the Bunch object)
        nominal_attributes = {k: v for k, v in arff['attributes']
                              if isinstance(v, list) and
                              k in data_columns + target_columns}

        X, y = _convert_arff_data(arff['data'], col_slice_x,
                                  col_slice_y, shape)

        is_classification = {col_name in nominal_attributes
                             for col_name in target_columns}
        if not is_classification:
            # No target
            pass
        elif all(is_classification):
            y = np.hstack([
                np.take(
                    np.asarray(nominal_attributes.pop(col_name), dtype='O'),
                    y[:, i:i + 1].astype(int, copy=False))
                for i, col_name in enumerate(target_columns)
            ])
        elif any(is_classification):
            raise ValueError('Mix of nominal and non-nominal targets is not '
                             'currently supported')

        # reshape y back to 1-D array, if there is only 1 target column; back
        # to None if there are not target columns
        if y.shape[1] == 1:
            y = y.reshape((-1,))
        elif y.shape[1] == 0:
            y = None

    bunch = Bunch(
        data=X, target=y, frame=frame, feature_names=data_columns,
        DESCR=description, details=data_description,
        categories=nominal_attributes,
        url="https://www.openml.org/d/{}".format(data_id)
    )
    return bunch
