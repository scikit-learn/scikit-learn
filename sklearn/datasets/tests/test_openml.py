"""Test the openml loader.

Skipped on travis.
"""

from sklearn.datasets import fetch_openml
from sklearn.utils.testing import (assert_warns_message,
                                   assert_raise_message)


def test_fetch_openml():
    # check_skip_travis()
    # fetch with version
    iris_1 = fetch_openml("iris", version=1)
    assert iris_1.details['id'] == '61'
    # fetch without version
    iris_1 = fetch_openml("iris")
    assert iris_1.details['id'] == '61'
    # fetch with dataset id
    iris_by_id = fetch_openml(61)
    assert iris_by_id.details['name'] == "iris"
    assert iris_by_id.data.shape == (150, 4)
    assert iris_by_id.target.shape == (150,)
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
