# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
Test the utility functions.
"""
import os
import shutil
import time
from pathlib import Path
import tempfile
from tempfile import TemporaryDirectory
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

import pytest

from ..utils import (
    parse_url,
    make_local_storage,
    temporary_file,
    unique_file_name,
)


def test_unique_name_long():
    "The file name should never be longer than 255 characters"
    url = f"https://www.something.com/data{'a' * 500}.txt"
    assert len(url) > 255
    fname = unique_file_name(url)
    assert len(fname) == 255
    assert fname[-10:] == "aaaaaa.txt"
    assert fname.split("-")[1][:10] == "aaaaaaaaaa"


@pytest.mark.parametrize(
    "pool",
    [ThreadPoolExecutor, ProcessPoolExecutor],
    ids=["threads", "processes"],
)
def test_make_local_storage_parallel(pool, monkeypatch):
    "Try to create the cache folder in parallel"
    # Can cause multiple attempts at creating the folder which leads to an
    # exception. Check that this doesn't happen.
    # See https://github.com/fatiando/pooch/issues/170

    # Monkey path makedirs to make it delay before creating the directory.
    # Otherwise, the dispatch is too fast and the directory will exist before
    # another process tries to create it.

    # Need to keep a reference to the original function to avoid infinite
    # recursions from the monkey patching.
    makedirs = os.makedirs

    def mockmakedirs(path, exist_ok=False):  # pylint: disable=unused-argument
        "Delay before calling makedirs"
        time.sleep(1.5)
        makedirs(path, exist_ok=exist_ok)

    monkeypatch.setattr(os, "makedirs", mockmakedirs)

    data_cache = os.path.join(os.curdir, "test_parallel_cache")
    assert not os.path.exists(data_cache)

    try:
        with pool() as executor:
            futures = [
                executor.submit(make_local_storage, data_cache) for i in range(4)
            ]
            for future in futures:
                future.result()
            assert os.path.exists(data_cache)
    finally:
        if os.path.exists(data_cache):
            shutil.rmtree(data_cache)


def test_local_storage_makedirs_permissionerror(monkeypatch):
    "Should warn the user when can't create the local data dir"

    def mockmakedirs(path, exist_ok=False):  # pylint: disable=unused-argument
        "Raise an exception to mimic permission issues"
        raise PermissionError("Fake error")

    data_cache = os.path.join(os.curdir, "test_permission")
    assert not os.path.exists(data_cache)

    monkeypatch.setattr(os, "makedirs", mockmakedirs)

    with pytest.raises(PermissionError) as error:
        make_local_storage(
            path=data_cache,
            env="SOME_VARIABLE",
        )
        assert "Pooch could not create data cache" in str(error)
        assert "'SOME_VARIABLE'" in str(error)


def test_local_storage_newfile_permissionerror(monkeypatch):
    "Should warn the user when can't write to the local data dir"
    # This is a separate function because there should be a warning if the data
    # dir already exists but we can't write to it.

    def mocktempfile(**kwargs):  # pylint: disable=unused-argument
        "Raise an exception to mimic permission issues"
        raise PermissionError("Fake error")

    with TemporaryDirectory() as data_cache:
        os.makedirs(os.path.join(data_cache, "1.0"))
        assert os.path.exists(data_cache)

        monkeypatch.setattr(tempfile, "NamedTemporaryFile", mocktempfile)

        with pytest.raises(PermissionError) as error:
            make_local_storage(
                path=data_cache,
                env="SOME_VARIABLE",
            )
            assert "Pooch could not write to data cache" in str(error)
            assert "'SOME_VARIABLE'" in str(error)


@pytest.mark.parametrize(
    "url,output",
    [
        (
            "http://127.0.0.1:8080/test.nc",
            {"protocol": "http", "netloc": "127.0.0.1:8080", "path": "/test.nc"},
        ),
        (
            "ftp://127.0.0.1:8080/test.nc",
            {"protocol": "ftp", "netloc": "127.0.0.1:8080", "path": "/test.nc"},
        ),
        (
            "doi:10.6084/m9.figshare.923450.v1/dike.json",
            {
                "protocol": "doi",
                "netloc": "10.6084/m9.figshare.923450.v1",
                "path": "/dike.json",
            },
        ),
        (
            r"doi:10.5281/zenodo.7632643/santisoler/pooch-test-data-v1.zip",
            {
                "protocol": "doi",
                "netloc": "10.5281/zenodo.7632643",
                "path": "/santisoler/pooch-test-data-v1.zip",
            },
        ),
    ],
    ids=["http", "ftp", "doi", "zenodo-doi-with-slash"],
)
def test_parse_url(url, output):
    "Parse URL into 3 components"
    assert parse_url(url) == output


def test_parse_url_invalid_doi():
    "Should fail if we forget to not include // in the DOI link"
    with pytest.raises(ValueError):
        parse_url("doi://XXX/XXX/fname.txt")


def test_temporary_file():
    "Make sure the file is writable and cleaned up in the end"
    with temporary_file() as tmp:
        assert Path(tmp).exists()
        with open(tmp, "w", encoding="utf-8") as outfile:
            outfile.write("Meh")
        with open(tmp, encoding="utf-8") as infile:
            assert infile.read().strip() == "Meh"
    assert not Path(tmp).exists()


def test_temporary_file_path():
    "Make sure the file is writable and cleaned up in the end when given a dir"
    with TemporaryDirectory() as path:
        with temporary_file(path) as tmp:
            assert Path(tmp).exists()
            assert path in tmp
            with open(tmp, "w", encoding="utf-8") as outfile:
                outfile.write("Meh")
            with open(tmp, encoding="utf-8") as infile:
                assert infile.read().strip() == "Meh"
        assert not Path(tmp).exists()


def test_temporary_file_exception():
    "Make sure the file is writable and cleaned up when there is an exception"
    try:
        with temporary_file() as tmp:
            assert Path(tmp).exists()
            raise ValueError("Nooooooooo!")
    except ValueError:
        assert not Path(tmp).exists()
