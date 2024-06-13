"""
Smoke Test the check_build module
"""

# Author: G Varoquaux
# SPDX-License-Identifier: BSD-3-Clause

import pytest

from sklearn.__check_build import raise_build_error


def test_raise_build_error():
    with pytest.raises(ImportError):
        raise_build_error(ImportError())
