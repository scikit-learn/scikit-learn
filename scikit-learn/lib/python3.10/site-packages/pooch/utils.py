# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
Misc utilities
"""
import logging
import os
import tempfile
import hashlib
from pathlib import Path
from urllib.parse import urlsplit
from contextlib import contextmanager
import warnings

import platformdirs
from packaging.version import Version


LOGGER = logging.Logger("pooch")
LOGGER.addHandler(logging.StreamHandler())


def file_hash(*args, **kwargs):
    """
    WARNING: Importing this function from pooch.utils is DEPRECATED.
    Please import from the top-level namespace (`from pooch import file_hash`)
    instead, which is fully backwards compatible with pooch >= 0.1.

    Examples
    --------

    >>> fname = "test-file-for-hash.txt"
    >>> with open(fname, "w") as f:
    ...     __ = f.write("content of the file")
    >>> print(file_hash(fname))
    0fc74468e6a9a829f103d069aeb2bb4f8646bad58bf146bb0e3379b759ec4a00
    >>> import os
    >>> os.remove(fname)

    """
    # pylint: disable=import-outside-toplevel
    from .hashes import file_hash as new_file_hash

    message = """
    Importing file_hash from pooch.utils is DEPRECATED. Please import from the
    top-level namespace (`from pooch import file_hash`) instead, which is fully
    backwards compatible with pooch >= 0.1.
    """
    warnings.warn(message, DeprecationWarning, stacklevel=2)
    return new_file_hash(*args, **kwargs)


def get_logger():
    r"""
    Get the default event logger.

    The logger records events like downloading files, unzipping archives, etc.
    Use the method :meth:`logging.Logger.setLevel` of this object to adjust the
    verbosity level from Pooch.

    Returns
    -------
    logger : :class:`logging.Logger`
        The logger object for Pooch
    """
    return LOGGER


def os_cache(project):
    r"""
    Default cache location based on the operating system.

    The folder locations are defined by the ``platformdirs``  package
    using the ``user_cache_dir`` function.
    Usually, the locations will be following (see the
    `platformdirs documentation <https://platformdirs.readthedocs.io>`__):

    * Mac: ``~/Library/Caches/<AppName>``
    * Unix: ``~/.cache/<AppName>`` or the value of the ``XDG_CACHE_HOME``
      environment variable, if defined.
    * Windows: ``C:\Users\<user>\AppData\Local\<AppAuthor>\<AppName>\Cache``

    Parameters
    ----------
    project : str
        The project name.

    Returns
    -------
    cache_path : :class:`pathlib.Path`
        The default location for the data cache. User directories (``'~'``) are
        not expanded.

    """
    return Path(platformdirs.user_cache_dir(project))


def check_version(version, fallback="master"):
    """
    Check if a version is PEP440 compliant and there are no unreleased changes.

    For example, ``version = "0.1"`` will be returned as is but ``version =
    "0.1+10.8dl8dh9"`` will return the fallback. This is the convention used by
    `versioneer <https://github.com/warner/python-versioneer>`__ to mark that
    this version is 10 commits ahead of the last release.

    Parameters
    ----------
    version : str
        A version string.
    fallback : str
        What to return if the version string has unreleased changes.

    Returns
    -------
    version : str
        If *version* is PEP440 compliant and there are unreleased changes, then
        return *version*. Otherwise, return *fallback*.

    Raises
    ------
    InvalidVersion
        If *version* is not PEP440 compliant.

    Examples
    --------

    >>> check_version("0.1")
    '0.1'
    >>> check_version("0.1a10")
    '0.1a10'
    >>> check_version("0.1+111.9hdg36")
    'master'
    >>> check_version("0.1+111.9hdg36", fallback="dev")
    'dev'

    """
    parse = Version(version)
    if parse.local is not None:
        return fallback
    return version


def parse_url(url):
    """
    Parse a URL into 3 components:

    <protocol>://<netloc>/<path>

    Example URLs:

    * http://127.0.0.1:8080/test.nc
    * ftp://127.0.0.1:8080/test.nc
    * doi:10.6084/m9.figshare.923450.v1/test.nc

    The DOI is a special case. The protocol will be "doi", the netloc will be
    the DOI, and the path is what comes after the last "/".
    The only exception are Zenodo dois: the protocol will be "doi", the netloc
    will be composed by the "prefix/suffix" and the path is what comes after
    the second "/". This allows to support special cases of Zenodo dois where
    the path contains forward slashes "/", created by the GitHub-Zenodo
    integration service.

    Parameters
    ----------
    url : str
        The URL.

    Returns
    -------
    parsed_url : dict
        Three components of a URL (e.g.,
        ``{'protocol':'http', 'netloc':'127.0.0.1:8080','path': '/test.nc'}``).

    """
    if url.startswith("doi://"):
        raise ValueError(
            f"Invalid DOI link '{url}'. You must not use '//' after 'doi:'."
        )
    if url.startswith("doi:"):
        protocol = "doi"
        parts = url[4:].split("/")
        if "zenodo" in parts[1].lower():
            netloc = "/".join(parts[:2])
            path = "/" + "/".join(parts[2:])
        else:
            netloc = "/".join(parts[:-1])
            path = "/" + parts[-1]
    else:
        parsed_url = urlsplit(url)
        protocol = parsed_url.scheme or "file"
        netloc = parsed_url.netloc
        path = parsed_url.path
    return {"protocol": protocol, "netloc": netloc, "path": path}


def cache_location(path, env=None, version=None):
    """
    Location of the cache given a base path and optional configuration.

    Checks for the environment variable to overwrite the path of the local
    cache. Optionally add *version* to the path if given.

    Parameters
    ----------
    path : str, PathLike, list or tuple
        The path to the local data storage folder. If this is a list or tuple,
        we'll join the parts with the appropriate separator. Use
        :func:`pooch.os_cache` for a sensible default.
    version : str or None
        The version string for your project. Will be appended to given path if
        not None.
    env : str or None
        An environment variable that can be used to overwrite *path*. This
        allows users to control where they want the data to be stored. We'll
        append *version* to the end of this value as well.

    Returns
    -------
    local_path : PathLike
        The path to the local directory.

    """
    if env is not None and env in os.environ and os.environ[env]:
        path = os.environ[env]
    if isinstance(path, (list, tuple)):
        path = os.path.join(*path)
    if version is not None:
        path = os.path.join(str(path), version)
    path = os.path.expanduser(str(path))
    return Path(path)


def make_local_storage(path, env=None):
    """
    Create the local cache directory and make sure it's writable.

    Parameters
    ----------
    path : str or PathLike
        The path to the local data storage folder.
    env : str or None
        An environment variable that can be used to overwrite *path*. Only used
        in the error message in case the folder is not writable.
    """
    path = str(path)
    # Check that the data directory is writable
    if not os.path.exists(path):
        action = "create"
    else:
        action = "write to"

    try:
        if action == "create":
            # When running in parallel, it's possible that multiple jobs will
            # try to create the path at the same time. Use exist_ok to avoid
            # raising an error.
            os.makedirs(path, exist_ok=True)
        else:
            with tempfile.NamedTemporaryFile(dir=path):
                pass
    except PermissionError as error:
        message = [
            str(error),
            f"| Pooch could not {action} data cache folder '{path}'.",
            "Will not be able to download data files.",
        ]
        if env is not None:
            message.append(
                f"Use environment variable '{env}' to specify a different location."
            )
        raise PermissionError(" ".join(message)) from error


@contextmanager
def temporary_file(path=None):
    """
    Create a closed and named temporary file and make sure it's cleaned up.

    Using :class:`tempfile.NamedTemporaryFile` will fail on Windows if trying
    to open the file a second time (when passing its name to Pooch function,
    for example). This context manager creates the file, closes it, yields the
    file path, and makes sure it's deleted in the end.

    Parameters
    ----------
    path : str or PathLike
        The directory in which the temporary file will be created.

    Yields
    ------
    fname : str
        The path to the temporary file.

    """
    tmp = tempfile.NamedTemporaryFile(delete=False, dir=path)
    # Close the temp file so that it can be opened elsewhere
    tmp.close()
    try:
        yield tmp.name
    finally:
        if os.path.exists(tmp.name):
            os.remove(tmp.name)


def unique_file_name(url):
    """
    Create a unique file name based on the given URL.

    The file name will be unique to the URL by prepending the name with the MD5
    hash (hex digest) of the URL. The name will also include the last portion
    of the URL.

    The format will be: ``{md5}-{filename}.{ext}``

    The file name will be cropped so that the entire name (including the hash)
    is less than 255 characters long (the limit on most file systems).

    Parameters
    ----------
    url : str
        The URL with a file name at the end.

    Returns
    -------
    fname : str
        The file name, unique to this URL.

    Examples
    --------

    >>> print(unique_file_name("https://www.some-server.org/2020/data.txt"))
    02ddee027ce5ebb3d7059fb23d210604-data.txt
    >>> print(unique_file_name("https://www.some-server.org/2019/data.txt"))
    9780092867b497fca6fc87d8308f1025-data.txt
    >>> print(unique_file_name("https://www.some-server.org/2020/data.txt.gz"))
    181a9d52e908219c2076f55145d6a344-data.txt.gz

    """
    md5 = hashlib.md5(url.encode()).hexdigest()
    fname = parse_url(url)["path"].split("/")[-1]
    # Crop the start of the file name to fit 255 characters including the hash
    # and the :
    fname = fname[-(255 - len(md5) - 1) :]
    unique_name = f"{md5}-{fname}"
    return unique_name
