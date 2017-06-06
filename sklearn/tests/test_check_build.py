"""
Smoke Test the check_build module
"""

# Author: G Varoquaux
# License: BSD 3 clause

from sklearn.__check_build import raise_build_error

from sklearn.utils.testing import assert_raises


def test_raise_build_error():
    assert_raises(ImportError, raise_build_error, ImportError())
