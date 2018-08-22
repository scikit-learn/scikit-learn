"""Test the openml loader.
"""
import gzip
import json
import numpy as np
import os
import re

import pytest
import scipy.sparse
import sklearn

from sklearn.datasets import fetch_openml
from sklearn.datasets.openml import (_open_openml_url,
                                     _get_data_description_by_id,
                                     _download_data_arff)
from sklearn.utils.testing import (assert_warns_message,
                                   assert_raise_message)
from sklearn.externals.six import string_types
from sklearn.externals.six.moves.urllib.error import HTTPError


currdir = os.path.dirname(os.path.abspath(__file__))
# if True, urlopen will be monkey patched to only use local files
test_offline = True
test_gzip = True


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
    return data_by_id


def _monkey_patch_webbased_functions(context, data_id, gziped_files,
                                     falsify_checksum=False):
    url_prefix_data_description = "https://openml.org/api/v1/json/data/"
    url_prefix_data_features = "https://openml.org/api/v1/json/data/features/"
    url_prefix_download_data = "https://openml.org/data/v1/"
    url_prefix_data_list = "https://openml.org/api/v1/json/data/list/"

    path_suffix = ''
    read_fn = open
    if gziped_files:
        path_suffix = '.gz'
        read_fn = gzip.open

    def _file_name(url, suffix):
        return (re.sub(r'\W', '-', url[len("https://openml.org/"):])
                + suffix + path_suffix)

    def _mock_urlopen_data_description(url):
        assert url.startswith(url_prefix_data_description)

        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            _file_name(url, '.json'))
        return read_fn(path, 'rb')

    def _mock_urlopen_data_features(url):
        assert url.startswith(url_prefix_data_features)

        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            _file_name(url, '.json'))
        return read_fn(path, 'rb')

    def _mock_urlopen_download_data(url):
        assert (url.startswith(url_prefix_download_data))

        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            _file_name(url, '.arff'))
        return read_fn(path, 'rb')

    def _mock_urlopen_data_list(url):
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
        return read_fn(json_file_path, 'rb')

    def _mock_urlopen(url):
        if url.startswith(url_prefix_data_list):
            return _mock_urlopen_data_list(url)
        elif url.startswith(url_prefix_data_features):
            return _mock_urlopen_data_features(url)
        elif url.startswith(url_prefix_download_data):
            return _mock_urlopen_download_data(url)
        elif url.startswith(url_prefix_data_description):
            return _mock_urlopen_data_description(url)
        else:
            raise ValueError('Unknown mocking URL pattern: %s' % url)

    # XXX: Global variable
    if test_offline:
        context.setattr(sklearn.datasets.openml, 'urlopen', _mock_urlopen)


def test_fetch_openml_iris(monkeypatch):
    # classification dataset with numeric only columns
    data_id = 61
    data_name = 'iris'
    data_version = 1
    target_column = 'class'
    expected_observations = 150
    expected_features = 4
    expected_missing = 0

    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               np.float64, object, expect_sparse=False,
                               compare_default_target=True)


def test_decode_iris():
    data_id = 61
    _test_features_list(data_id)


def test_fetch_openml_iris_multitarget(monkeypatch):
    # classification dataset with numeric only columns
    data_id = 61
    data_name = 'iris'
    data_version = 1
    target_column = ['sepallength', 'sepalwidth']
    expected_observations = 150
    expected_features = 3
    expected_missing = 0

    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, np.float64, expect_sparse=False,
                               compare_default_target=False)


def test_fetch_openml_anneal(monkeypatch):
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    target_column = 'class'
    # Not all original instances included for space reasons
    expected_observations = 11
    expected_features = 38
    expected_missing = 267
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, object, expect_sparse=False,
                               compare_default_target=True)


def test_decode_anneal():
    data_id = 2
    _test_features_list(data_id)


def test_fetch_openml_anneal_multitarget(monkeypatch):
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    target_column = ['class', 'product-type', 'shape']
    # Not all original instances included for space reasons
    expected_observations = 11
    expected_features = 36
    expected_missing = 267
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, object, expect_sparse=False,
                               compare_default_target=False)


def test_fetch_openml_cpu(monkeypatch):
    # regression dataset with numeric and categorical columns
    data_id = 561
    data_name = 'cpu'
    data_version = 1
    target_column = 'class'
    expected_observations = 209
    expected_features = 7
    expected_missing = 0
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               object, np.float64, expect_sparse=False,
                               compare_default_target=True)


def test_decode_cpu():
    data_id = 561
    _test_features_list(data_id)


def test_fetch_openml_australian(monkeypatch):
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
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
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


def test_fetch_openml_miceprotein(monkeypatch):
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
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               np.float64, object, expect_sparse=False,
                               compare_default_target=True)


def test_fetch_openml_emotions(monkeypatch):
    # classification dataset with multiple targets (natively)
    data_id = 40589
    data_name = 'emotions'
    data_version = 3
    target_column = ['amazed.suprised', 'happy.pleased', 'relaxing.calm',
                     'quiet.still', 'sad.lonely', 'angry.aggresive']
    expected_observations = 13
    expected_features = 72
    expected_missing = 0
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)

    _fetch_dataset_from_openml(data_id, data_name, data_version, target_column,
                               expected_observations, expected_features,
                               expected_missing,
                               np.float64, object, expect_sparse=False,
                               compare_default_target=True)


def test_decode_emotions():
    data_id = 40589
    _test_features_list(data_id)


def test_open_openml_url_cache(monkeypatch):
    data_id = 61

    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
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


def test_fetch_openml_notarget(monkeypatch):
    data_id = 61
    target_column = None
    expected_observations = 150
    expected_features = 5

    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    data = fetch_openml(data_id=data_id, target_column=target_column,
                        cache=False)
    assert data.data.shape == (expected_observations, expected_features)
    assert data.target is None


def test_fetch_openml_inactive(monkeypatch):
    # fetch inactive dataset by id
    data_id = 40675
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    glas2 = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        data_id=data_id, cache=False)
    # fetch inactive dataset by name and version
    assert glas2.data.shape == (163, 9)
    glas2_by_version = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        data_id=None, name="glass2", version=1, cache=False)
    assert int(glas2_by_version.details['id']) == data_id


def test_fetch_nonexiting(monkeypatch):
    # there is no active version of glass2
    data_id = 40675
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    # Note that we only want to search by name (not data id)
    assert_raise_message(ValueError, "No active dataset glass2 found",
                         fetch_openml, name='glass2', cache=False)


def test_raises_illegal_multitarget(monkeypatch):
    data_id = 61
    targets = ['sepalwidth', 'class']
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    # Note that we only want to search by name (not data id)
    assert_raise_message(ValueError,
                         "Can only handle homogeneous multi-target datasets,",
                         fetch_openml, data_id=data_id,
                         target_column=targets, cache=False)


def test_warn_ignore_attribute(monkeypatch):
    data_id = 40966
    expected_row_id_msg = "target_column={} has flag is_row_identifier."
    expected_ignore_msg = "target_column={} has flag is_ignore."
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
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


def test_string_attribute(monkeypatch):
    data_id = 40945
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    # single column test
    assert_raise_message(ValueError,
                         'STRING attributes are not yet supported',
                         fetch_openml, data_id=data_id, cache=False)


def test_illegal_column(monkeypatch):
    data_id = 61
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
    assert_raise_message(KeyError, "Could not find target_column=",
                         fetch_openml, data_id=data_id,
                         target_column='undefined', cache=False)

    assert_raise_message(KeyError, "Could not find target_column=",
                         fetch_openml, data_id=data_id,
                         target_column=['undefined', 'class'],
                         cache=False)


def test_fetch_openml_raises_missing_values_target(monkeypatch):
    data_id = 2
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)
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


###########################################################################
    # Validate MD5 Checksum
###########################################################################

def _monkey_patch_checksum_data_description(context, data_id,
                                            gziped_files, false_checksum):

    path_suffix = ''
    read_fn = open
    if gziped_files:
        path_suffix = '.gz'
        read_fn = gzip.open

    def mock_data_description(data_id, data_home):
        path = os.path.join(currdir, 'data', 'openml', str(data_id),
                            'api-v1-json-data-' + str(data_id) + '.json' +
                            path_suffix)

        fp = read_fn(path)

        json_data = json.load(fp)
        data_description = json_data['data_set_description']
        # Falsify the checksum to purposefully get a warning
        data_description['md5_checksum'] = false_checksum

        return data_description

    context.setattr(sklearn.datasets.openml, '_get_data_description_by_id',
                    mock_data_description)


@pytest.mark.parametrize("data_id, true_checksum, false_checksum, cache", [
    (61, 'ad484452702105cbf3d30f8deaba39a9',
     'ad484452702105cbf3d30f8deaba39a8', True),

    (61, 'ad484452702105cbf3d30f8deaba39a9',
     'ad484452702105cbf3d30f8deaba39a8', False),

    (561, 'e1c69097976ecd20de7d215919130ccc',
     'e1c69097976ecd20de7d215919130ccd', True),

    (561, 'e1c69097976ecd20de7d215919130ccc',
     'e1c69097976ecd20de7d215919130ccd', False),
])
def test_fetch_openml_checksum_invalid(monkeypatch, tmpdir, data_id,
                                       true_checksum, false_checksum, cache):

    warn_message = 'Data set file hash {} does not match the checksum {}.'
    warn_message = warn_message.format(true_checksum, false_checksum)

    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip,
                                     falsify_checksum=True)
    _monkey_patch_checksum_data_description(monkeypatch, data_id, test_gzip,
                                            false_checksum)
    assert_warns_message(UserWarning, warn_message, fetch_openml,
                         data_id=data_id, data_home=str(tmpdir), cache=cache)


@pytest.mark.parametrize('data_id, cache', [
    (61, True),
    (61, False),
    (561, True),
    (561, False)
])
def test_fetch_openml_checksum_valid(monkeypatch, tmpdir, data_id, cache):
    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip)

    # Capture all warnings
    with pytest.warns(None) as records:
        fetch_openml(data_id=data_id, data_home=str(tmpdir), cache=cache)
        # assert no warnings
        assert not records


@pytest.mark.parametrize('data_id, cache', [
    (61, True),
    (61, False),
    (561, True),
    (561, False)
])
def test_fetch_openml_checksum_invalid_no_verification(monkeypatch, tmpdir,
                                                       data_id, cache):

    _monkey_patch_webbased_functions(monkeypatch, data_id, test_gzip,
                                     falsify_checksum=True)
    # Capture all warnings
    with pytest.warns(None) as records:
        fetch_openml(data_id=data_id, cache=cache, data_home=str(tmpdir),
                     verify_checksum=False)
        # assert no warnings
        assert not records
