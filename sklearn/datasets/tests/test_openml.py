"""Test the openml loader.
"""
import gzip
import json
import numpy as np
import os
import re
import scipy.sparse
import sklearn
import pytest

from sklearn.datasets import fetch_openml
from sklearn.datasets.openml import (_open_openml_url,
                                     _get_data_description_by_id,
                                     _download_data_arff)
from sklearn.utils.testing import (assert_warns_message,
                                   assert_raise_message)
from sklearn.externals.six import string_types
from sklearn.externals.six.moves.urllib.error import HTTPError
from sklearn.datasets.tests.test_common import check_return_X_y
from functools import partial


currdir = os.path.dirname(os.path.abspath(__file__))
# if True, urlopen will be monkey patched to only use local files
test_offline = True


def _test_features_list(data_id):
    # XXX Test is intended to verify/ensure correct decoding behavior
    # Not usable with sparse data or datasets that have columns marked as
    # {row_identifier, ignore}
    def decode_column(data_bunch, col_idx):
        col_name = data_bunch.feature_names[col_idx]
        if col_name in data_bunch.categories:
            # XXX: This would be faster with np.take, although it does not
            # handle missing values fast (also not with mode='wrap')
            cat = data_bunch.categories[col_name]
            result = [cat[idx] if 0 <= idx < len(cat) else None for idx in
                      data_bunch.data[:, col_idx].astype(int)]
            return np.array(result, dtype='O')
        else:
            # non-nominal attribute
            return data_bunch.data[:, col_idx]

    data_bunch = fetch_openml(data_id=data_id, cache=False, target_column=None)

    # also obtain decoded arff
    data_description = _get_data_description_by_id(data_id, None)
    sparse = data_description['format'].lower() == 'sparse_arff'
    if sparse is True:
        raise ValueError('This test is not intended for sparse data, to keep '
                         'code relatively simple')
    data_arff = _download_data_arff(data_description['file_id'],
                                    sparse, None, False)
    data_downloaded = np.array(data_arff['data'], dtype='O')

    for i in range(len(data_bunch.feature_names)):
        # XXX: Test per column, as this makes it easier to avoid problems with
        # missing values

        np.testing.assert_array_equal(data_downloaded[:, i],
                                      decode_column(data_bunch, i))


def _fetch_dataset_from_openml(data_id, data_name, data_version,
                               target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               expected_data_dtype, expected_target_dtype,
                               expect_sparse, compare_default_target):
    # fetches a dataset in three various ways from OpenML, using the
    # fetch_openml function, and does various checks on the validity of the
    # result. Note that this function can be mocked (by invoking
    # _monkey_patch_webbased_functions before invoking this function)
    data_by_name_id = fetch_openml(name=data_name, version=data_version,
                                   cache=False)
    assert int(data_by_name_id.details['id']) == data_id

    fetch_openml(name=data_name, cache=False)
    # without specifying the version, there is no guarantee that the data id
    # will be the same

    # fetch with dataset id
    data_by_id = fetch_openml(data_id=data_id, cache=False,
                              target_column=target_column)
    assert data_by_id.details['name'] == data_name
    assert data_by_id.data.shape == (expected_observations, expected_features)
    if isinstance(target_column, str):
        # single target, so target is vector
        assert data_by_id.target.shape == (expected_observations, )
    elif isinstance(target_column, list):
        # multi target, so target is array
        assert data_by_id.target.shape == (expected_observations,
                                           len(target_column))
    assert data_by_id.data.dtype == np.float64
    assert data_by_id.target.dtype == expected_target_dtype
    assert len(data_by_id.feature_names) == expected_features
    for feature in data_by_id.feature_names:
        assert isinstance(feature, string_types)

    # TODO: pass in a list of expected nominal features
    for feature, categories in data_by_id.categories.items():
        feature_idx = data_by_id.feature_names.index(feature)
        values = np.unique(data_by_id.data[:, feature_idx])
        values = values[np.isfinite(values)]
        assert set(values) <= set(range(len(categories)))

    if compare_default_target:
        # check whether the data by id and data by id target are equal
        data_by_id_default = fetch_openml(data_id=data_id, cache=False)
        if data_by_id.data.dtype == np.float64:
            np.testing.assert_allclose(data_by_id.data,
                                       data_by_id_default.data)
        else:
            assert np.array_equal(data_by_id.data, data_by_id_default.data)
        if data_by_id.target.dtype == np.float64:
            np.testing.assert_allclose(data_by_id.target,
                                       data_by_id_default.target)
        else:
            assert np.array_equal(data_by_id.target, data_by_id_default.target)

    if expect_sparse:
        assert isinstance(data_by_id.data, scipy.sparse.csr_matrix)
    else:
        assert isinstance(data_by_id.data, np.ndarray)
        # np.isnan doesn't work on CSR matrix
        assert (np.count_nonzero(np.isnan(data_by_id.data)) ==
                expected_missing)

    # test return_X_y option
    fetch_func = partial(fetch_openml, data_id=data_id, cache=False,
                         target_column=target_column)
    check_return_X_y(data_by_id, fetch_func)
    return data_by_id


def _monkey_patch_webbased_functions(context,
                                     data_id,
                                     gzip_response):
    url_prefix_data_description = "https://openml.org/api/v1/json/data/"
    url_prefix_data_features = "https://openml.org/api/v1/json/data/features/"
    url_prefix_download_data = "https://openml.org/data/v1/"
    url_prefix_data_list = "https://openml.org/api/v1/json/data/list/"

    path_suffix = '.gz'
    read_fn = gzip.open

    class MockHTTPResponse(object):
        def __init__(self, data, is_gzip):
            self.data = data
            self.is_gzip = is_gzip

        def read(self, amt=-1):
            return self.data.read(amt)

        def tell(self):
            return self.data.tell()

        def seek(self, pos, whence=0):
            return self.data.seek(pos, whence)

        def close(self):
            self.data.close()

        def info(self):
            if self.is_gzip:
                return {'Content-Encoding': 'gzip'}
            return {}

    def _file_name(url, suffix):
        return (re.sub(r'\W', '-', url[len("https://openml.org/"):])
                + suffix + path_suffix)

    def _mock_urlopen_data_description(url, has_gzip_header):
        assert url.startswith(url_prefix_data_description)

        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            _file_name(url, '.json'))

        if has_gzip_header and gzip_response:
            fp = open(path, 'rb')
            return MockHTTPResponse(fp, True)
        else:
            fp = read_fn(path, 'rb')
            return MockHTTPResponse(fp, False)

    def _mock_urlopen_data_features(url, has_gzip_header):
        assert url.startswith(url_prefix_data_features)
        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            _file_name(url, '.json'))
        if has_gzip_header and gzip_response:
            fp = open(path, 'rb')
            return MockHTTPResponse(fp, True)
        else:
            fp = read_fn(path, 'rb')
            return MockHTTPResponse(fp, False)

    def _mock_urlopen_download_data(url, has_gzip_header):
        assert (url.startswith(url_prefix_download_data))

        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            _file_name(url, '.arff'))

        if has_gzip_header and gzip_response:
            fp = open(path, 'rb')
            return MockHTTPResponse(fp, True)
        else:
            fp = read_fn(path, 'rb')
            return MockHTTPResponse(fp, False)

    def _mock_urlopen_data_list(url, has_gzip_header):
        assert url.startswith(url_prefix_data_list)

        json_file_path = os.path.join(currdir, 'data', 'openml',
                                      str(data_id), _file_name(url, '.json'))
        # load the file itself, to simulate a http error
        json_data = json.loads(read_fn(json_file_path, 'rb').
                               read().decode('utf-8'))
        if 'error' in json_data:
            raise HTTPError(url=None, code=412,
                            msg='Simulated mock error',
                            hdrs=None, fp=None)

        if has_gzip_header:
            fp = open(json_file_path, 'rb')
            return MockHTTPResponse(fp, True)
        else:
            fp = read_fn(json_file_path, 'rb')
            return MockHTTPResponse(fp, False)

    def _mock_urlopen(request):
        url = request.get_full_url()
        has_gzip_header = request.get_header('Accept-encoding') == "gzip"
        if url.startswith(url_prefix_data_list):
            return _mock_urlopen_data_list(url, has_gzip_header)
        elif url.startswith(url_prefix_data_features):
            return _mock_urlopen_data_features(url, has_gzip_header)
        elif url.startswith(url_prefix_download_data):
            return _mock_urlopen_download_data(url, has_gzip_header)
        elif url.startswith(url_prefix_data_description):
            return _mock_urlopen_data_description(url, has_gzip_header)
        else:
            raise ValueError('Unknown mocking URL pattern: %s' % url)

    # XXX: Global variable
    if test_offline:
        context.setattr(sklearn.datasets.openml, 'urlopen', _mock_urlopen)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_iris(monkeypatch, gzip_response):
    # classification dataset with numeric only columns
    data_id = 61
    data_name = 'iris'
    data_version = 1
    target_column = 'class'
    expected_observations = 150
    expected_features = 4
    expected_missing = 0

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    assert_warns_message(
        UserWarning,
        "Multiple active versions of the dataset matching the name"
        " iris exist. Versions may be fundamentally different, "
        "returning version 1.",
        _fetch_dataset_from_openml,
        **{'data_id': data_id, 'data_name': data_name,
           'data_version': data_version,
           'target_column': target_column,
           'expected_observations': expected_observations,
           'expected_features': expected_features,
           'expected_missing': expected_missing,
           'expect_sparse': False,
           'expected_data_dtype': np.float64,
           'expected_target_dtype': object,
           'compare_default_target': True}
    )


def test_decode_iris(monkeypatch):
    data_id = 61
    _monkey_patch_webbased_functions(monkeypatch, data_id, False)
    _test_features_list(data_id)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_iris_multitarget(monkeypatch, gzip_response):
    # classification dataset with numeric only columns
    data_id = 61
    data_name = 'iris'
    data_version = 1
    target_column = ['sepallength', 'sepalwidth']
    expected_observations = 150
    expected_features = 3
    expected_missing = 0

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, np.float64, expect_sparse=False,
                               compare_default_target=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_anneal(monkeypatch, gzip_response):
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    target_column = 'class'
    # Not all original instances included for space reasons
    expected_observations = 11
    expected_features = 38
    expected_missing = 267
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, object, expect_sparse=False,
                               compare_default_target=True)


def test_decode_anneal(monkeypatch):
    data_id = 2
    _monkey_patch_webbased_functions(monkeypatch, data_id, False)
    _test_features_list(data_id)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_anneal_multitarget(monkeypatch, gzip_response):
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    target_column = ['class', 'product-type', 'shape']
    # Not all original instances included for space reasons
    expected_observations = 11
    expected_features = 36
    expected_missing = 267
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, object, expect_sparse=False,
                               compare_default_target=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_cpu(monkeypatch, gzip_response):
    # regression dataset with numeric and categorical columns
    data_id = 561
    data_name = 'cpu'
    data_version = 1
    target_column = 'class'
    expected_observations = 209
    expected_features = 7
    expected_missing = 0
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, np.float64, expect_sparse=False,
                               compare_default_target=True)


def test_decode_cpu(monkeypatch):
    data_id = 561
    _monkey_patch_webbased_functions(monkeypatch, data_id, False)
    _test_features_list(data_id)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_australian(monkeypatch, gzip_response):
    # sparse dataset
    # Australian is the only sparse dataset that is reasonably small
    # as it is inactive, we need to catch the warning. Due to mocking
    # framework, it is not deactivated in our tests
    data_id = 292
    data_name = 'Australian'
    data_version = 1
    target_column = 'Y'
    # Not all original instances included for space reasons
    expected_observations = 85
    expected_features = 14
    expected_missing = 0
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    assert_warns_message(
        UserWarning,
        "Version 1 of dataset Australian is inactive,",
        _fetch_dataset_from_openml,
        **{'data_id': data_id, 'data_name': data_name,
           'data_version': data_version,
           'target_column': target_column,
           'expected_observations': expected_observations,
           'expected_features': expected_features,
           'expected_missing': expected_missing,
           'expect_sparse': True,
           'expected_data_dtype': np.float64,
           'expected_target_dtype': object,
           'compare_default_target': False}  # numpy specific check
    )


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_miceprotein(monkeypatch, gzip_response):
    # JvR: very important check, as this dataset defined several row ids
    # and ignore attributes. Note that data_features json has 82 attributes,
    # and row id (1), ignore attributes (3) have been removed (and target is
    # stored in data.target)
    data_id = 40966
    data_name = 'MiceProtein'
    data_version = 4
    target_column = 'class'
    # Not all original instances included for space reasons
    expected_observations = 7
    expected_features = 77
    expected_missing = 7
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               np.float64, object, expect_sparse=False,
                               compare_default_target=True)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_emotions(monkeypatch, gzip_response):
    # classification dataset with multiple targets (natively)
    data_id = 40589
    data_name = 'emotions'
    data_version = 3
    target_column = ['amazed.suprised', 'happy.pleased', 'relaxing.calm',
                     'quiet.still', 'sad.lonely', 'angry.aggresive']
    expected_observations = 13
    expected_features = 72
    expected_missing = 0
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)

    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               np.float64, object, expect_sparse=False,
                               compare_default_target=True)


def test_decode_emotions(monkeypatch):
    data_id = 40589
    _monkey_patch_webbased_functions(monkeypatch, data_id, False)
    _test_features_list(data_id)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_open_openml_url_cache(monkeypatch, gzip_response):
    data_id = 61

    _monkey_patch_webbased_functions(
        monkeypatch, data_id, gzip_response)
    openml_path = sklearn.datasets.openml._DATA_FILE.format(data_id)
    test_directory = os.path.join(os.path.expanduser('~'), 'scikit_learn_data')
    # first fill the cache
    response1 = _open_openml_url(openml_path, test_directory)
    # assert file exists
    location = os.path.join(test_directory, 'openml.org', openml_path + '.gz')
    assert os.path.isfile(location)
    # redownload, to utilize cache
    response2 = _open_openml_url(openml_path, test_directory)
    assert response1.read() == response2.read()


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_notarget(monkeypatch, gzip_response):
    data_id = 61
    target_column = None
    expected_observations = 150
    expected_features = 5

    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    data = fetch_openml(data_id=data_id, target_column=target_column,
                        cache=False)
    assert data.data.shape == (expected_observations, expected_features)
    assert data.target is None


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_inactive(monkeypatch, gzip_response):
    # fetch inactive dataset by id
    data_id = 40675
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    glas2 = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        data_id=data_id, cache=False)
    # fetch inactive dataset by name and version
    assert glas2.data.shape == (163, 9)
    glas2_by_version = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        data_id=None, name="glass2", version=1, cache=False)
    assert int(glas2_by_version.details['id']) == data_id


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_nonexiting(monkeypatch, gzip_response):
    # there is no active version of glass2
    data_id = 40675
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    # Note that we only want to search by name (not data id)
    assert_raise_message(ValueError, "No active dataset glass2 found",
                         fetch_openml, name='glass2', cache=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_raises_illegal_multitarget(monkeypatch, gzip_response):
    data_id = 61
    targets = ['sepalwidth', 'class']
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    # Note that we only want to search by name (not data id)
    assert_raise_message(ValueError,
                         "Can only handle homogeneous multi-target datasets,",
                         fetch_openml, data_id=data_id,
                         target_column=targets, cache=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_warn_ignore_attribute(monkeypatch, gzip_response):
    data_id = 40966
    expected_row_id_msg = "target_column={} has flag is_row_identifier."
    expected_ignore_msg = "target_column={} has flag is_ignore."
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    # single column test
    assert_warns_message(UserWarning, expected_row_id_msg.format('MouseID'),
                         fetch_openml, data_id=data_id,
                         target_column='MouseID',
                         cache=False)
    assert_warns_message(UserWarning, expected_ignore_msg.format('Genotype'),
                         fetch_openml, data_id=data_id,
                         target_column='Genotype',
                         cache=False)
    # multi column test
    assert_warns_message(UserWarning, expected_row_id_msg.format('MouseID'),
                         fetch_openml, data_id=data_id,
                         target_column=['MouseID', 'class'],
                         cache=False)
    assert_warns_message(UserWarning, expected_ignore_msg.format('Genotype'),
                         fetch_openml, data_id=data_id,
                         target_column=['Genotype', 'class'],
                         cache=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_string_attribute(monkeypatch, gzip_response):
    data_id = 40945
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    # single column test
    assert_raise_message(ValueError,
                         'STRING attributes are not yet supported',
                         fetch_openml, data_id=data_id, cache=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_illegal_column(monkeypatch, gzip_response):
    data_id = 61
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    assert_raise_message(KeyError, "Could not find target_column=",
                         fetch_openml, data_id=data_id,
                         target_column='undefined', cache=False)

    assert_raise_message(KeyError, "Could not find target_column=",
                         fetch_openml, data_id=data_id,
                         target_column=['undefined', 'class'],
                         cache=False)


@pytest.mark.parametrize('gzip_response', [True, False])
def test_fetch_openml_raises_missing_values_target(monkeypatch, gzip_response):
    data_id = 2
    _monkey_patch_webbased_functions(monkeypatch, data_id, gzip_response)
    assert_raise_message(ValueError, "Target column ",
                         fetch_openml, data_id=data_id, target_column='family')


def test_fetch_openml_raises_illegal_argument():
    assert_raise_message(ValueError, "Dataset data_id=",
                         fetch_openml, data_id=-1, name="name")

    assert_raise_message(ValueError, "Dataset data_id=",
                         fetch_openml, data_id=-1, name=None,
                         version="version")

    assert_raise_message(ValueError, "Dataset data_id=",
                         fetch_openml, data_id=-1, name="name",
                         version="version")

    assert_raise_message(ValueError, "Neither name nor data_id are provided. "
                         "Please provide name or data_id.", fetch_openml)
