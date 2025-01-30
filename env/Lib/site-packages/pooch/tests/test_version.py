# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
Test the version.
"""
from packaging.version import Version

import pooch


def test_version():
    "Check there's a usable version number in the usual __version__"
    assert pooch.__version__.startswith("v")
    # Check that it's PEP440 compliant (will raise an exception otherwise)
    Version(pooch.__version__)
