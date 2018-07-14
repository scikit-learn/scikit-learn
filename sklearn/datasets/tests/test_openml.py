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
from sklearn.externals._liacarff.arff import load


def fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse):
    # fetch with version
    data_by_name_id = fetch_openml(name=data_name, version=data_version)
    assert int(data_by_name_id.details['id']) == data_id

    # fetch without version
    fetch_openml(name=data_name)
    # without specifying the version, there is no guarantee that the data id
    # will be the same

    # fetch with dataset id
    data_by_id = fetch_openml(data_id)
    assert data_by_id.details['name'] == data_name
    assert data_by_id.data.shape == (expected_observations, expected_features)
    assert data_by_id.target.shape == (expected_observations, )
    if expect_sparse:
        assert isinstance(data_by_id.data, scipy.sparse.csr_matrix)
    else:
        assert isinstance(data_by_id.data, np.ndarray)

    # check numeric features:
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

    def _mock_data_description(id):
        path = os.path.join(testdir_path,
                           'mock_openml/%d/data_description.json' % id)
        description = open(path, 'r').read()
        return json.loads(description)['data_set_description']

    def _mock_data_features(id):
        path = os.path.join(testdir_path,
                            'mock_openml/%d/data_features.json' % id)
        features = open(path, 'r').read()
        return json.loads(features)['data_features']['feature']

    def _mock_download_data(_):
        path = os.path.join(testdir_path, 'mock_openml/%d/data.arff' % data_id)
        arff_fp = open(path, 'r').read()
        return load(arff_fp)

    def _mock_get_data_info_by_name(name_or_id, version):
        path = os.path.join(testdir_path,
                            'mock_openml/%d/%s_%s.json' % (data_id,
                                                           name_or_id,
                                                           version))
        data_info = open(path, 'r').read()
        data_info_json = json.loads(data_info)
        return data_info_json['data']['dataset'][0]

    context.setattr(sklearn.datasets.openml,
                    '_get_data_description_by_id',
                    _mock_data_description)
    context.setattr(sklearn.datasets.openml,
                    '_get_data_features',
                    _mock_data_features)
    context.setattr(sklearn.datasets.openml, '_download_data',
                    _mock_download_data)
    context.setattr(sklearn.datasets.openml,
                    '_get_data_info_by_name',
                    _mock_get_data_info_by_name)


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


def test_fetch_openml_anneal():
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    expected_observations = 898
    expected_features = 38
    fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse=False)


def test_fetch_openml_cpu():
    # regression dataset with numeric and categorical columns
    data_id = 561
    data_name = 'cpu'
    data_version = 1
    expected_observations = 209
    expected_features = 7
    fetch_dataset_from_openml(data_id, data_name, data_version,
                              expected_observations, expected_features,
                              expect_sparse=False)


def test_fetch_openml_australian():
    # sparse dataset
    # Australian is the only sparse dataset that is reasonably small
    # as it is inactive, we need to catch the warning
    data_id = 292
    data_name = 'Australian'
    data_version = 1
    expected_observations = 690
    expected_features = 14
    assert_warns_message(
        UserWarning,
        "Version 1 of dataset Australian is inactive,",
        fetch_dataset_from_openml,
        **{'data_id': data_id, 'data_name': data_name,
           'data_version': data_version,
           'expected_observations': expected_observations,
           'expected_features': expected_features,
           'expect_sparse': False}
        # Sadly, due to a bug in liac-arff library, the data
        # is always returned as a dense array.
        # discussion in OpenML library:
        # https://github.com/openml/openml-python/issues/487
    )


def test_fetch_openml_inactive():
    # fetch inactive dataset by id
    glas2 = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        40675)
    # fetch inactive dataset by name and version
    assert glas2.data.shape == (163, 9)
    glas2_by_version = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        None, "glass2", 1)
    assert glas2_by_version.details['id'] == '40675'


def test_fetch_nonexiting():
    # there is no active version of glass2
    assert_raise_message(ValueError, "No active dataset glass2 found",
                         fetch_openml, None, 'glass2')


def test_fetch_openml_raises_illegal_argument():
    assert_raise_message(ValueError, "Dataset id=",
                         fetch_openml, -1, "name")

    assert_raise_message(ValueError, "Dataset id=",
                         fetch_openml, -1, None, "version")

    assert_raise_message(ValueError, "Dataset id=",
                         fetch_openml, -1, "name", "version")

    assert_raise_message(ValueError, "Neither name nor id are provided. " +
                         "Please provide name xor id.", fetch_openml)
