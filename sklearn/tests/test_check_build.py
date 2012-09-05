"""
Smoke Test the check_build module
"""
# This

# Author: G Varoquaux
# License: BSD style

from sklearn.__check_build import raise_build_error

from nose.tools import assert_raises


def test_raise_build_error():
    assert_raises(ImportError, raise_build_error, ImportError())
