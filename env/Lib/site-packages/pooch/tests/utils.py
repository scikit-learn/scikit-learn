# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
Utilities for testing code.
"""
import os
import io
import logging
import shutil
import stat
from pathlib import Path
from contextlib import contextmanager

from .. import __version__ as full_version
from ..utils import check_version, get_logger


def check_tiny_data(fname):
    """
    Load the tiny-data.txt file and check that the contents are correct.
    """
    assert os.path.exists(fname)
    with open(fname, encoding="utf-8") as tinydata:
        content = tinydata.read()
    true_content = "\n".join(
        ["# A tiny data file for test purposes only", "1  2  3  4  5  6"]
    )
    assert content.strip() == true_content


def check_large_data(fname):
    """
    Load the large-data.txt file and check that the contents are correct.
    """
    assert os.path.exists(fname)
    with open(fname, encoding="utf-8") as data:
        content = data.read()
    true_content = ["# A larer data file for test purposes only"]
    true_content.extend(["1  2  3  4  5  6"] * 6002)
    assert content.strip() == "\n".join(true_content)


def pooch_test_url():
    """
    Get the base URL for the test data used in Pooch itself.

    The URL is a GitHub raw link to the ``pooch/tests/data`` directory from the
    `GitHub repository <https://github.com/fatiando/pooch>`__. It matches the
    pooch version specified in ``pooch.version.full_version``.

    Returns
    -------
    url
        The versioned URL for pooch's test data.

    """
    version = check_version(full_version, fallback="main")
    url = f"https://github.com/fatiando/pooch/raw/{version}/pooch/tests/data/"
    return url


def pooch_test_figshare_url():
    """
    Get the base URL for the test data stored in figshare.

    The URL contains the DOI for the figshare dataset using the appropriate
    version for this version of Pooch.

    Returns
    -------
    url
        The URL for pooch's test data.

    """
    url = "doi:10.6084/m9.figshare.14763051.v1/"
    return url


def pooch_test_zenodo_url():
    """
    Get the base URL for the test data stored in Zenodo.

    The URL contains the DOI for the Zenodo dataset using the appropriate
    version for this version of Pooch.

    Returns
    -------
    url
        The URL for pooch's test data.

    """
    url = "doi:10.5281/zenodo.4924875/"
    return url


def pooch_test_zenodo_with_slash_url():
    """
    Get base URL for test data in Zenodo, where the file name contains a slash

    The URL contains the DOI for the Zenodo dataset that has a slash in the
    filename (created with the GitHub-Zenodo integration service), using the
    appropriate version for this version of Pooch.

    Returns
    -------
    url
        The URL for pooch's test data.

    """
    url = "doi:10.5281/zenodo.7632643/"
    return url


def pooch_test_dataverse_url():
    """
    Get the base URL for the test data stored on a DataVerse instance.

    Returns
    -------
    url
        The URL for pooch's test data.
    """
    url = "doi:10.11588/data/TKCFEF/"
    return url


def pooch_test_registry():
    """
    Get a registry for the test data used in Pooch itself.

    Returns
    -------
    registry
        Dictionary with pooch's test data files and their hashes.

    """
    registry = {
        "tiny-data.txt": "baee0894dba14b12085eacb204284b97e362f4f3e5a5807693cc90ef415c1b2d",
        "large-data.txt": "98de171fb320da82982e6bf0f3994189fff4b42b23328769afce12bdd340444a",
        "subdir/tiny-data.txt": "baee0894dba14b12085eacb204284b97e362f4f3e5a5807693cc90ef415c1b2d",
        "tiny-data.zip": "0d49e94f07bc1866ec57e7fd1b93a351fba36842ec9b13dd50bf94e8dfa35cbb",
        "store.zip": "0498d2a001e71051bbd2acd2346f38da7cbd345a633cb7bf0f8a20938714b51a",
        "tiny-data.tar.gz": "41503f083814f43a01a8e9a30c28d7a9fe96839a99727a7fdd0acf7cd5bab63b",
        "store.tar.gz": "088c7f4e0f1859b1c769bb6065de24376f366374817ede8691a6ac2e49f29511",
        "tiny-data.txt.bz2": "753663687a4040c90c8578061867d1df623e6aa8011c870a5dbd88ee3c82e306",
        "tiny-data.txt.gz": "2e2da6161291657617c32192dba95635706af80c6e7335750812907b58fd4b52",
        "tiny-data.txt.xz": "99dcb5c32a6e916344bacb4badcbc2f2b6ee196977d1d8187610c21e7e607765",
    }
    return registry


@contextmanager
def capture_log(level=logging.DEBUG):
    """
    Create a context manager for reading from the logs.

    Yields
    ------
    log_file : StringIO
        a file-like object to which the logs were written
    """
    log_file = io.StringIO()
    handler = logging.StreamHandler(log_file)
    handler.setLevel(level)
    get_logger().addHandler(handler)
    yield log_file
    get_logger().removeHandler(handler)


@contextmanager
def data_over_ftp(server, fname):
    """
    Add a test data file to the test FTP server and clean it up afterwards.

    Parameters
    ----------
    server
        The ``ftpserver`` fixture provided by pytest-localftpserver.
    fname : str
        The name of a file *relative* to the test data folder of the package
        (usually just the file name, not the full path).

    Yields
    ------
    url : str
        The download URL of the data file from the test FTP server.

    """
    package_path = str(Path(__file__).parent / "data" / fname)
    server_path = os.path.join(server.anon_root, fname)
    try:
        shutil.copyfile(package_path, server_path)
        url = f"ftp://localhost/{fname}"
        yield url
    finally:
        if os.path.exists(server_path):
            os.remove(server_path)


def _recursive_chmod_directories(root, mode):
    """
    Recursively change the permissions on the child directories using a bitwise
    OR operation.
    """
    for item in root.iterdir():
        if item.is_dir():
            item.chmod(item.stat().st_mode | mode)
            _recursive_chmod_directories(item, mode)


def mirror_directory(source, destination):
    """
    Copy contents of the source directory into destination and fix permissions.

    Parameters
    ----------
    source : str, :class:`pathlib.Path`
        Source data directory.
    destination : str, :class:`pathlib.Path`
        Destination directory that will contain the copy of source. The actual
        source directory (not just it's contents) is copied.

    Returns
    -------
    mirror : :class:`pathlib.Path`
        The path of the mirrored output directory.

    """
    source = Path(source)
    mirror = Path(destination) / source.name
    shutil.copytree(source, mirror)
    _recursive_chmod_directories(mirror, mode=stat.S_IWUSR)
    return mirror
