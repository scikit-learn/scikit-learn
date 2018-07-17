"""Test the openml loader.
"""
import json
import numpy as np
import os
import scipy.sparse
import sklearn

from sklearn.datasets import fetch_openml
from sklearn.datasets.openml import _get_data_features
from sklearn.utils.testing import (assert_warns_message,
                                   assert_raise_message)
from sklearn.externals.six.moves.urllib.error import HTTPError


try:
    # Python 2
    from urllib2 import urlopen
except ImportError:
    # Python 3+
    from urllib.request import urlopen


def fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse):
    data_by_name_id = fetch_openml(name=data_name, version=data_version,
                                   cache=False)
    assert int(data_by_name_id.details['id']) == data_id

    fetch_openml(name=data_name, cache=False)
    # without specifying the version, there is no guarantee that the data id
    # will be the same

    # fetch with dataset id
    data_by_id = fetch_openml(data_id=data_id, cache=False)
    assert data_by_id.details['name'] == data_name
    assert data_by_id.data.shape == (expected_observations, expected_features)
    assert data_by_id.target.shape == (expected_observations, )
    if expect_sparse:
        assert isinstance(data_by_id.data, scipy.sparse.csr_matrix)
    else:
        assert isinstance(data_by_id.data, np.ndarray)

    # check numeric features. Note that the response of _get_data_features is
    # mocked too.
    feature_name_type = {feature['name']: feature['data_type']
                         for feature in _get_data_features(data_id)}
    for idx, feature_name in enumerate(data_by_id.feature_names):
        if feature_name_type[feature_name] == 'numeric':
            # casting trick according to Jaime at
            # stackoverflow.com/questions/19486283/
            # how-do-i-quickly-check-if-all-elements-of-numpy-array-are-floats
            assert np.issubdtype(np.array(list(data_by_id.data[:, idx])).dtype,
                                 np.number)

    if 'default_target_attribute' in data_by_id.details:
        target = data_by_id.details['default_target_attribute']
        if feature_name_type[target] == 'numeric':
            assert np.issubdtype(np.array(list(data_by_id.data[:, idx])).dtype,
                                 np.number)
    return data_by_id


def _monkey_patch_webbased_functions(context, data_id):
    testdir_path = os.path.dirname(os.path.realpath(__file__))
    url_prefix_data_description = "https://openml.org/api/v1/json/data/"
    url_prefix_data_features = "https://openml.org/api/v1/json/data/features/"
    url_prefix_download_data = "https://openml.org/data/v1/"
    url_prefix_data_list = "https://openml.org/api/v1/json/data/list/"

    def _mock_urlopen_data_description(url):
        assert (url.startswith(url_prefix_data_description))

        path = os.path.join(testdir_path,
                            'mock_openml/%d/data_description.json' % data_id)
        return urlopen('file://' + path)

    def _mock_urlopen_data_features(url):
        assert (url.startswith(url_prefix_data_features))

        path = os.path.join(testdir_path,
                            'mock_openml/%d/data_features.json' % data_id)
        return urlopen('file://' + path)

    def _mock_urlopen_download_data(url):
        assert (url.startswith(url_prefix_download_data))

        path = os.path.join(testdir_path, 'mock_openml/%d/data.arff' % data_id)
        return urlopen('file://' + path)

    def _mock_urlopen_data_list(url):
        # url contains key value pairs of attributes, e.g.,
        # openml.org/api/v1/json/data_name/iris/data_version/1 should
        # ideally become {data_name: 'iris', data_version: '1'}
        assert(url.startswith(url_prefix_data_list))
        att_list = url[len(url_prefix_data_list):].split('/')
        key_val_dict = dict(zip(att_list[::2], att_list[1::2]))
        # add defaults, so we can make assumptions about the content
        if 'data_version' not in key_val_dict:
            key_val_dict['data_version'] = None
        if 'status' not in key_val_dict:
            key_val_dict['status'] = "active"
        mock_file = "%s_%s_%s.json" % (key_val_dict['data_name'],
                                       key_val_dict['data_version'],
                                       key_val_dict['status'])
        json_file_path = os.path.join(testdir_path, 'mock_openml',
                                      str(data_id), mock_file)
        # load the file itself, to simulate a http error
        json_data = json.loads(open(json_file_path, 'r').read())
        if 'error' in json_data:
            raise HTTPError(url=None, code=412,
                            msg='Simulated mock error',
                            hdrs=None, fp=None)
        return urlopen('file://' + json_file_path)

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
            # This should never happen
            raise ValueError('Unknown mocking URL pattern: %s' % url)

    context.setattr(sklearn.datasets.openml, 'urlopen', _mock_urlopen)


def test_fetch_openml_iris(monkeypatch):
    # classification dataset with numeric only columns
    data_id = 61
    data_name = 'iris'
    data_version = 1
    expected_observations = 150
    expected_features = 4

    _monkey_patch_webbased_functions(monkeypatch, data_id)
    fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse=False)


def test_fetch_openml_anneal(monkeypatch):
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    # Not all original instances included for space reasons
    expected_observations = 11
    expected_features = 38
    _monkey_patch_webbased_functions(monkeypatch, data_id)
    fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse=False)


def test_fetch_openml_cpu(monkeypatch):
    # regression dataset with numeric and categorical columns
    data_id = 561
    data_name = 'cpu'
    data_version = 1
    expected_observations = 209
    expected_features = 7
    _monkey_patch_webbased_functions(monkeypatch, data_id)
    fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse=False)


def test_fetch_openml_australian(monkeypatch):
    # sparse dataset
    # Australian is the only sparse dataset that is reasonably small
    # as it is inactive, we need to catch the warning. Due to mocking
    # framework, it is not deactivated in our tests
    data_id = 292
    data_name = 'Australian'
    data_version = 1
    # Not all original instances included for space reasons
    expected_observations = 85
    expected_features = 14
    _monkey_patch_webbased_functions(monkeypatch, data_id)
    assert_warns_message(
        UserWarning,
        "Version 1 of dataset Australian is inactive,",
        fetch_dataset_from_openml,
        **{'data_id': data_id, 'data_name': data_name,
           'data_version': data_version,
           'expected_observations': expected_observations,
           'expected_features': expected_features,
           'expect_sparse': False}
    )


def test_fetch_openml_miceprotein(monkeypatch):
    data_id = 40966
    data_name = 'MiceProtein'
    data_version = 4
    # Not all original instances included for space reasons
    expected_observations = 7
    expected_features = 77
    _monkey_patch_webbased_functions(monkeypatch, data_id)
    data = fetch_dataset_from_openml(data_id, data_name, data_version,
                                     expected_observations, expected_features,
                                     expect_sparse=False)
    # JvR: very important check, as this dataset defined several row ids
    # and ignore attributes. Note that data_features json has 82 attributes,
    # and row id (1), ignore attributes (3) have been removed (and target is
    # stored in data.target)
    assert (data.data.shape == (expected_observations, expected_features))
    assert (data.target.shape == (expected_observations, ))
    assert (data.data.dtype == np.float64)


def test_fetch_openml_inactive(monkeypatch):
    # makes contact with openml server. not mocked
    # fetch inactive dataset by id
    data_id = 40675
    _monkey_patch_webbased_functions(monkeypatch, data_id)
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
    _monkey_patch_webbased_functions(monkeypatch, data_id)
    # Note that we only want to search by name (not data id)
    assert_raise_message(ValueError, "No active dataset glass2 found",
                         fetch_openml, name='glass2', cache=False)


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
