# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
# pylint: disable=missing-docstring,import-outside-toplevel,import-self
#
# Import functions/classes to make the API
from .core import Pooch, create, retrieve
from .utils import os_cache, check_version, get_logger
from .hashes import file_hash, make_registry
from .downloaders import (
    HTTPDownloader,
    FTPDownloader,
    SFTPDownloader,
    DOIDownloader,
)
from .processors import Unzip, Untar, Decompress

# This file is generated automatically by setuptools_scm
from . import _version


# Add a "v" to the version number
__version__ = f"v{_version.version}"


def test(doctest=True, verbose=True, coverage=False):
    """
    Run the test suite.

    Uses `py.test <http://pytest.org/>`__ to discover and run the tests.

    Parameters
    ----------

    doctest : bool
        If ``True``, will run the doctests as well (code examples that start
        with a ``>>>`` in the docs).
    verbose : bool
        If ``True``, will print extra information during the test run.
    coverage : bool
        If ``True``, will run test coverage analysis on the code as well.
        Requires ``pytest-cov``.

    Raises
    ------

    AssertionError
        If pytest returns a non-zero error code indicating that some tests have
        failed.

    """
    import pytest

    package = __name__
    args = []
    if verbose:
        args.append("-vv")
    if coverage:
        args.append(f"--cov={package}")
        args.append("--cov-report=term-missing")
    if doctest:
        args.append("--doctest-modules")
    args.append("--pyargs")
    args.append(package)
    status = pytest.main(args)
    assert status == 0, "Some tests have failed."
