"""Test the openml loader.
"""
import numpy as np

from sklearn.datasets import fetch_openml
from sklearn.datasets.openml import _get_data_features
from sklearn.utils.testing import (assert_warns_message,
                                   assert_raise_message)


def fetch_dataset_from_openml(data_id, data_name, data_version, expected_observations, expected_features):
    # fetch with version
    data_by_name_id = fetch_openml(data_name, version=data_version)
    assert int(data_by_name_id.details['id']) == data_id

    # fetch without version
    data_by_name = fetch_openml(data_name)
    assert int(data_by_name.details['id']) == data_id

    # fetch with dataset id
    data_by_id = fetch_openml(data_id)
    assert data_by_id.details['name'] == data_name
    assert data_by_id.data.shape == (expected_observations, expected_features)
    assert data_by_id.target.shape == (expected_observations, )

    # check numeric features:
    feature_name_type = {feature['name']: feature['data_type'] for feature in _get_data_features(data_id)}
    for idx, feature_name in enumerate(data_by_id.feature_names):
        if feature_name_type[feature_name] == 'numeric':
            # casting trick according to Jaime at
            # stackoverflow.com/questions/19486283/how-do-i-quickly-check-if-all-elements-of-numpy-array-are-floats
            assert np.issubdtype(np.array(list(data_by_id.data[:, idx])).dtype, np.number)

    if 'default_target_attribute' in data_by_id.details:
        target = data_by_id.details['default_target_attribute']
        if feature_name_type[target] == 'numeric':
            assert np.issubdtype(np.array(list(data_by_id.data[:, idx])).dtype, np.number)


def test_fetch_openml_iris():
    # classification dataset with numeric only columns
    data_id = 61
    data_name = 'iris'
    data_version = 1
    expected_observations = 150
    expected_features = 4
    fetch_dataset_from_openml(data_id, data_name, data_version, expected_observations, expected_features)


def test_fetch_openml_anneal():
    # classification dataset with numeric and categorical columns
    data_id = 2
    data_name = 'anneal'
    data_version = 1
    expected_observations = 898
    expected_features = 38
    fetch_dataset_from_openml(data_id, data_name, data_version, expected_observations, expected_features)


def test_fetch_openml_cpu():
    # regression dataset with numeric and categorical columns
    data_id = 561
    data_name = 'cpu'
    data_version = 1
    expected_observations = 209
    expected_features = 7
    fetch_dataset_from_openml(data_id, data_name, data_version, expected_observations, expected_features)


def test_fetch_openml_inactive():
    # fetch inactive dataset by id
    glas2 = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        40675)
    # fetch inactive dataset by name and version
    assert glas2.data.shape == (163, 9)
    glas2_by_version = assert_warns_message(
        UserWarning, "Version 1 of dataset glass2 is inactive,", fetch_openml,
        "glass2", 1)
    # there is no active version of glass2
    assert glas2_by_version.details['id'] == '40675'
    assert_raise_message(ValueError, "No active dataset glass2 found",
                         fetch_openml, 'glass2')
