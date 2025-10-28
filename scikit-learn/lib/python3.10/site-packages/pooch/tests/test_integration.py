# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
# pylint: disable=redefined-outer-name
"""
Test the entire process of creating a Pooch and using it.
"""
import os
import shutil
from pathlib import Path

import pytest

from .. import create, os_cache
from .. import __version__ as full_version
from .utils import check_tiny_data, capture_log


@pytest.mark.network
def test_create_and_fetch():
    "Fetch a data file from the local storage"
    path = os_cache("pooch-testing")
    if path.exists():
        shutil.rmtree(str(path))
    pup = create(
        path=path,
        base_url="https://github.com/fatiando/pooch/raw/{version}/data/",
        version=full_version,
        version_dev="main",
        env="POOCH_DATA_DIR",
    )
    # Make sure the storage isn't created until a download is required
    assert not pup.abspath.exists()
    pup.load_registry(Path(os.path.dirname(__file__), "data", "registry.txt"))
    for target in ["tiny-data.txt", "subdir/tiny-data.txt"]:
        with capture_log() as log_file:
            fname = pup.fetch(target)
            assert log_file.getvalue().split()[0] == "Downloading"
        check_tiny_data(fname)
        # Now modify the file to trigger an update on the next fetch
        with open(fname, "w", encoding="utf-8") as fin:
            fin.write("The data is now different")
        with capture_log() as log_file:
            fname = pup.fetch(target)
            assert log_file.getvalue().split()[0] == "Updating"
        check_tiny_data(fname)
