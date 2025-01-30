# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
The main Pooch class and a factory function for it.
"""
import os
import time
import contextlib
from pathlib import Path
import shlex
import shutil


from .hashes import hash_matches, file_hash
from .utils import (
    check_version,
    get_logger,
    make_local_storage,
    cache_location,
    temporary_file,
    os_cache,
    unique_file_name,
)
from .downloaders import DOIDownloader, choose_downloader, doi_to_repository


def retrieve(
    url,
    known_hash,
    fname=None,
    path=None,
    processor=None,
    downloader=None,
    progressbar=False,
):
    """
    Download and cache a single file locally.

    Uses HTTP or FTP by default, depending on the protocol in the given *url*.
    Other download methods can be controlled through the *downloader* argument
    (see below).

    The file will be downloaded to a temporary location first and its hash will
    be compared to the given *known_hash*. This is done to ensure that the
    download happened correctly and securely. If the hash doesn't match, the
    file will be deleted and an exception will be raised.

    If the file already exists locally, its hash will be compared to
    *known_hash*. If they are not the same, this is interpreted as the file
    needing to be updated and it will be downloaded again.

    You can bypass these checks by passing ``known_hash=None``. If this is
    done, the SHA256 hash of the downloaded file will be logged to the screen.
    It is highly recommended that you copy and paste this hash as *known_hash*
    so that future downloads are guaranteed to be the exact same file. This is
    crucial for reproducible computations.

    If the file exists in the given *path* with the given *fname* and the hash
    matches, it will not be downloaded and the absolute path to the file will
    be returned.

    .. note::

        This function is meant for downloading single files. If you need to
        manage the download and caching of several files, with versioning, use
        :func:`pooch.create` and :class:`pooch.Pooch` instead.

    Parameters
    ----------
    url : str
        The URL to the file that is to be downloaded. Ideally, the URL should
        end in a file name.
    known_hash : str or None
        A known hash (checksum) of the file. Will be used to verify the
        download or check if an existing file needs to be updated. By default,
        will assume it's a SHA256 hash. To specify a different hashing method,
        prepend the hash with ``algorithm:``, for example
        ``md5:pw9co2iun29juoh`` or ``sha1:092odwhi2ujdp2du2od2odh2wod2``. If
        None, will NOT check the hash of the downloaded file or check if an
        existing file needs to be updated.
    fname : str or None
        The name that will be used to save the file. Should NOT include the
        full path, just the file name (it will be appended to *path*). If
        None, will create a unique file name using a combination of the last
        part of the URL (assuming it's the file name) and the MD5 hash of the
        URL. For example, ``81whdo2d2e928yd1wi22-data-file.csv``. This ensures
        that files from different URLs never overwrite each other, even if they
        have the same name.
    path : str or PathLike or None
        The location of the cache folder on disk. This is where the file will
        be saved. If None, will save to a ``pooch`` folder in the default cache
        location for your operating system (see :func:`pooch.os_cache`).
    processor : None or callable
        If not None, then a function (or callable object) that will be called
        before returning the full path and after the file has been downloaded
        (if required). See :ref:`processors` for details.
    downloader : None or callable
        If not None, then a function (or callable object) that will be called
        to download a given URL to a provided local file name. See
        :ref:`downloaders` for details.
    progressbar : bool or an arbitrary progress bar object
        If True, will print a progress bar of the download to standard error
        (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to be
        installed. Alternatively, an arbitrary progress bar object can be
        passed. See :ref:`custom-progressbar` for details.

    Returns
    -------
    full_path : str
        The absolute path (including the file name) of the file in the local
        storage.

    Examples
    --------

    Download one of the data files from the Pooch repository on GitHub:

    >>> import os
    >>> from pooch import __version__, check_version, retrieve
    >>> # Make a URL for the version of pooch we have installed
    >>> url = "https://github.com/fatiando/pooch/raw/{}/data/tiny-data.txt"
    >>> url = url.format(check_version(__version__, fallback="main"))
    >>> # Download the file and save it locally. Will check the MD5 checksum of
    >>> # the downloaded file against the given value to make sure it's the
    >>> # right file. You can use other hashes by specifying different
    >>> # algorithm names (sha256, sha1, etc).
    >>> fname = retrieve(
    ...     url, known_hash="md5:70e2afd3fd7e336ae478b1e740a5f08e",
    ... )
    >>> with open(fname) as f:
    ...     print(f.read().strip())
    # A tiny data file for test purposes only
    1  2  3  4  5  6
    >>> # Running again won't trigger a download and only return the path to
    >>> # the existing file.
    >>> fname2 = retrieve(
    ...     url, known_hash="md5:70e2afd3fd7e336ae478b1e740a5f08e",
    ... )
    >>> print(fname2 == fname)
    True
    >>> os.remove(fname)

    Files that are compressed with gzip, xz/lzma, or bzip2 can be automatically
    decompressed by passing using the :class:`pooch.Decompress` processor:

    >>> from pooch import Decompress
    >>> # URLs to a gzip compressed version of the data file.
    >>> url = ("https://github.com/fatiando/pooch/raw/{}/"
    ...        + "pooch/tests/data/tiny-data.txt.gz")
    >>> url = url.format(check_version(__version__, fallback="main"))
    >>> # By default, you would have to decompress the file yourself
    >>> fname = retrieve(
    ...     url,
    ...     known_hash="md5:8812ba10b6c7778014fdae81b03f9def",
    ... )
    >>> print(os.path.splitext(fname)[1])
    .gz
    >>> # Use the processor to decompress after download automatically and
    >>> # return the path to the decompressed file instead.
    >>> fname2 = retrieve(
    ...     url,
    ...     known_hash="md5:8812ba10b6c7778014fdae81b03f9def",
    ...     processor=Decompress(),
    ... )
    >>> print(fname2 == fname)
    False
    >>> with open(fname2) as f:
    ...     print(f.read().strip())
    # A tiny data file for test purposes only
    1  2  3  4  5  6
    >>> os.remove(fname)
    >>> os.remove(fname2)

    When downloading archives (zip or tar), it can be useful to unpack them
    after download to avoid having to do that yourself. Use the processors
    :class:`pooch.Unzip` or :class:`pooch.Untar` to do this automatically:

    >>> from pooch import Unzip
    >>> # URLs to a zip archive with a single data file.
    >>> url = ("https://github.com/fatiando/pooch/raw/{}/"
    ...        + "pooch/tests/data/tiny-data.zip")
    >>> url = url.format(check_version(__version__, fallback="main"))
    >>> # By default, you would get the path to the archive
    >>> fname = retrieve(
    ...     url,
    ...     known_hash="md5:e9592cb46cf3514a1079051f8a148148",
    ... )
    >>> print(os.path.splitext(fname)[1])
    .zip
    >>> os.remove(fname)
    >>> # Using the processor, the archive will be unzipped and a list with the
    >>> # path to every file will be returned instead of a single path.
    >>> fnames = retrieve(
    ...     url,
    ...     known_hash="md5:e9592cb46cf3514a1079051f8a148148",
    ...     processor=Unzip(),
    ... )
    >>> # There was only a single file in our archive.
    >>> print(len(fnames))
    1
    >>> with open(fnames[0]) as f:
    ...     print(f.read().strip())
    # A tiny data file for test purposes only
    1  2  3  4  5  6
    >>> for f in fnames:
    ...     os.remove(f)


    """
    if path is None:
        path = os_cache("pooch")
    if fname is None:
        fname = unique_file_name(url)
    # Make the path absolute.
    path = cache_location(path, env=None, version=None)

    full_path = path.resolve() / fname
    action, verb = download_action(full_path, known_hash)

    if action in ("download", "update"):
        # We need to write data, so create the local data directory if it
        # doesn't already exist.
        make_local_storage(path)

        get_logger().info(
            "%s data from '%s' to file '%s'.",
            verb,
            url,
            str(full_path),
        )

        if downloader is None:
            downloader = choose_downloader(url, progressbar=progressbar)

        stream_download(url, full_path, known_hash, downloader, pooch=None)

        if known_hash is None:
            get_logger().info(
                "SHA256 hash of downloaded file: %s\n"
                "Use this value as the 'known_hash' argument of 'pooch.retrieve'"
                " to ensure that the file hasn't changed if it is downloaded again"
                " in the future.",
                file_hash(str(full_path)),
            )

    if processor is not None:
        return processor(str(full_path), action, None)

    return str(full_path)


def create(
    path,
    base_url,
    version=None,
    version_dev="master",
    env=None,
    registry=None,
    urls=None,
    retry_if_failed=0,
    allow_updates=True,
):
    """
    Create a :class:`~pooch.Pooch` with sensible defaults to fetch data files.

    If a version string is given, the Pooch will be versioned, meaning that the
    local storage folder and the base URL depend on the project version. This
    is necessary if your users have multiple versions of your library installed
    (using virtual environments) and you updated the data files between
    versions. Otherwise, every time a user switches environments would trigger
    a re-download of the data. The version string will be appended to the local
    storage path (for example, ``~/.mypooch/cache/v0.1``) and inserted into the
    base URL (for example,
    ``https://github.com/fatiando/pooch/raw/v0.1/data``). If the version string
    contains ``+XX.XXXXX``, it will be interpreted as a development version.

    Does **not** create the local data storage folder. The folder will only be
    created the first time a download is attempted with
    :meth:`pooch.Pooch.fetch`. This makes it safe to use this function at the
    module level (so it's executed on ``import`` and the resulting
    :class:`~pooch.Pooch` is a global variable).

    Parameters
    ----------
    path : str, PathLike, list or tuple
        The path to the local data storage folder. If this is a list or tuple,
        we'll join the parts with the appropriate separator. The *version* will
        be appended to the end of this path. Use :func:`pooch.os_cache` for a
        sensible default.
    base_url : str
        Base URL for the remote data source. All requests will be made relative
        to this URL. The string should have a ``{version}`` formatting mark in
        it. We will call ``.format(version=version)`` on this string. If the
        URL does not end in a ``'/'``, a trailing ``'/'`` will be added
        automatically.
    version : str or None
        The version string for your project. Should be PEP440 compatible. If
        None is given, will not attempt to format *base_url* and no subfolder
        will be appended to *path*.
    version_dev : str
        The name used for the development version of a project. If your data is
        hosted on Github (and *base_url* is a Github raw link), then
        ``"master"`` is a good choice (default). Ignored if *version* is None.
    env : str or None
        An environment variable that can be used to overwrite *path*. This
        allows users to control where they want the data to be stored. We'll
        append *version* to the end of this value as well.
    registry : dict or None
        A record of the files that are managed by this Pooch. Keys should be
        the file names and the values should be their hashes. Only files
        in the registry can be fetched from the local storage. Files in
        subdirectories of *path* **must use Unix-style separators** (``'/'``)
        even on Windows.
    urls : dict or None
        Custom URLs for downloading individual files in the registry. A
        dictionary with the file names as keys and the custom URLs as values.
        Not all files in *registry* need an entry in *urls*. If a file has an
        entry in *urls*, the *base_url* will be ignored when downloading it in
        favor of ``urls[fname]``.
    retry_if_failed : int
        Retry a file download the specified number of times if it fails because
        of a bad connection or a hash mismatch. By default, downloads are only
        attempted once (``retry_if_failed=0``). Initially, will wait for 1s
        between retries and then increase the wait time by 1s with each retry
        until a maximum of 10s.
    allow_updates : bool or str
        Whether existing files in local storage that have a hash mismatch with
        the registry are allowed to update from the remote URL. If a string is
        passed, we will assume it's the name of an environment variable that
        will be checked for the true/false value. If ``False``, any mismatch
        with hashes in the registry will result in an error. Defaults to
        ``True``.

    Returns
    -------
    pooch : :class:`~pooch.Pooch`
        The :class:`~pooch.Pooch` initialized with the given arguments.

    Examples
    --------

    Create a :class:`~pooch.Pooch` for a release (v0.1):

    >>> pup = create(path="myproject",
    ...              base_url="http://some.link.com/{version}/",
    ...              version="v0.1",
    ...              registry={"data.txt": "9081wo2eb2gc0u..."})
    >>> print(pup.path.parts)  # The path is a pathlib.Path
    ('myproject', 'v0.1')
    >>> # The local folder is only created when a dataset is first downloaded
    >>> print(pup.path.exists())
    False
    >>> print(pup.base_url)
    http://some.link.com/v0.1/
    >>> print(pup.registry)
    {'data.txt': '9081wo2eb2gc0u...'}
    >>> print(pup.registry_files)
    ['data.txt']

    If this is a development version (12 commits ahead of v0.1), then the
    ``version_dev`` will be used (defaults to ``"master"``):

    >>> pup = create(path="myproject",
    ...              base_url="http://some.link.com/{version}/",
    ...              version="v0.1+12.do9iwd")
    >>> print(pup.path.parts)
    ('myproject', 'master')
    >>> print(pup.base_url)
    http://some.link.com/master/

    Versioning is optional (but highly encouraged):

    >>> pup = create(path="myproject",
    ...              base_url="http://some.link.com/",
    ...              registry={"data.txt": "9081wo2eb2gc0u..."})
    >>> print(pup.path.parts)  # The path is a pathlib.Path
    ('myproject',)
    >>> print(pup.base_url)
    http://some.link.com/

    To place the storage folder at a subdirectory, pass in a list and we'll
    join the path for you using the appropriate separator for your operating
    system:

    >>> pup = create(path=["myproject", "cache", "data"],
    ...              base_url="http://some.link.com/{version}/",
    ...              version="v0.1")
    >>> print(pup.path.parts)
    ('myproject', 'cache', 'data', 'v0.1')

    The user can overwrite the storage path by setting an environment variable:

    >>> # The variable is not set so we'll use *path*
    >>> pup = create(path=["myproject", "not_from_env"],
    ...              base_url="http://some.link.com/{version}/",
    ...              version="v0.1",
    ...              env="MYPROJECT_DATA_DIR")
    >>> print(pup.path.parts)
    ('myproject', 'not_from_env', 'v0.1')
    >>> # Set the environment variable and try again
    >>> import os
    >>> os.environ["MYPROJECT_DATA_DIR"] = os.path.join("myproject", "env")
    >>> pup = create(path=["myproject", "not_env"],
    ...              base_url="http://some.link.com/{version}/",
    ...              version="v0.1",
    ...              env="MYPROJECT_DATA_DIR")
    >>> print(pup.path.parts)
    ('myproject', 'env', 'v0.1')

    """
    if version is not None:
        version = check_version(version, fallback=version_dev)
        base_url = base_url.format(version=version)
    # Don't create the cache folder here! This function is usually called in
    # the module context (at import time), so touching the file system is not
    # recommended. It could cause crashes when multiple processes/threads try
    # to import at the same time (which would try to create the folder several
    # times at once).
    path = cache_location(path, env, version)
    if isinstance(allow_updates, str):
        allow_updates = os.environ.get(allow_updates, "true").lower() != "false"
    # add trailing "/"
    base_url = base_url.rstrip("/") + "/"
    pup = Pooch(
        path=path,
        base_url=base_url,
        registry=registry,
        urls=urls,
        retry_if_failed=retry_if_failed,
        allow_updates=allow_updates,
    )
    return pup


class Pooch:
    """
    Manager for a local data storage that can fetch from a remote source.

    Avoid creating ``Pooch`` instances directly. Use :func:`pooch.create`
    instead.

    Parameters
    ----------
    path : str
        The path to the local data storage folder. The path must exist in the
        file system.
    base_url : str
        Base URL for the remote data source. All requests will be made relative
        to this URL.
    registry : dict or None
        A record of the files that are managed by this good boy. Keys should be
        the file names and the values should be their hashes. Only files
        in the registry can be fetched from the local storage. Files in
        subdirectories of *path* **must use Unix-style separators** (``'/'``)
        even on Windows.
    urls : dict or None
        Custom URLs for downloading individual files in the registry. A
        dictionary with the file names as keys and the custom URLs as values.
        Not all files in *registry* need an entry in *urls*. If a file has an
        entry in *urls*, the *base_url* will be ignored when downloading it in
        favor of ``urls[fname]``.
    retry_if_failed : int
        Retry a file download the specified number of times if it fails because
        of a bad connection or a hash mismatch. By default, downloads are only
        attempted once (``retry_if_failed=0``). Initially, will wait for 1s
        between retries and then increase the wait time by 1s with each retry
        until a maximum of 10s.
    allow_updates : bool
        Whether existing files in local storage that have a hash mismatch with
        the registry are allowed to update from the remote URL. If ``False``,
        any mismatch with hashes in the registry will result in an error.
        Defaults to ``True``.

    """

    def __init__(
        self,
        path,
        base_url,
        registry=None,
        urls=None,
        retry_if_failed=0,
        allow_updates=True,
    ):
        self.path = path
        self.base_url = base_url
        if registry is None:
            registry = {}
        self.registry = registry
        if urls is None:
            urls = {}
        self.urls = dict(urls)
        self.retry_if_failed = retry_if_failed
        self.allow_updates = allow_updates

    @property
    def abspath(self):
        "Absolute path to the local storage"
        return Path(os.path.abspath(os.path.expanduser(str(self.path))))

    @property
    def registry_files(self):
        "List of file names on the registry"
        return list(self.registry)

    def fetch(self, fname, processor=None, downloader=None, progressbar=False):
        """
        Get the absolute path to a file in the local storage.

        If it's not in the local storage, it will be downloaded. If the hash of
        the file in local storage doesn't match the one in the registry, will
        download a new copy of the file. This is considered a sign that the
        file was updated in the remote storage. If the hash of the downloaded
        file still doesn't match the one in the registry, will raise an
        exception to warn of possible file corruption.

        Post-processing actions sometimes need to be taken on downloaded files
        (unzipping, conversion to a more efficient format, etc). If these
        actions are time or memory consuming, it would be best to do this only
        once right after the file is downloaded. Use the *processor* argument
        to specify a function that is executed after the download to perform
        these actions. See :ref:`processors` for details.

        Custom file downloaders can be provided through the *downloader*
        argument. By default, Pooch will determine the download protocol from
        the URL in the registry. If the server for a given file requires
        authentication (username and password), use a downloader that support
        these features. Downloaders can also be used to print custom messages
        (like a progress bar), etc. See :ref:`downloaders` for details.

        Parameters
        ----------
        fname : str
            The file name (relative to the *base_url* of the remote data
            storage) to fetch from the local storage.
        processor : None or callable
            If not None, then a function (or callable object) that will be
            called before returning the full path and after the file has been
            downloaded. See :ref:`processors` for details.
        downloader : None or callable
            If not None, then a function (or callable object) that will be
            called to download a given URL to a provided local file name. See
            :ref:`downloaders` for details.
        progressbar : bool or an arbitrary progress bar object
            If True, will print a progress bar of the download to standard
            error (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to
            be installed. Alternatively, an arbitrary progress bar object can
            be passed. See :ref:`custom-progressbar` for details.

        Returns
        -------
        full_path : str
            The absolute path (including the file name) of the file in the
            local storage.

        """
        self._assert_file_in_registry(fname)

        url = self.get_url(fname)
        full_path = self.abspath / fname
        known_hash = self.registry[fname]
        action, verb = download_action(full_path, known_hash)

        if action == "update" and not self.allow_updates:
            raise ValueError(
                f"{fname} needs to update {full_path} but updates are disallowed."
            )

        if action in ("download", "update"):
            # We need to write data, so create the local data directory if it
            # doesn't already exist.
            make_local_storage(str(self.abspath))

            get_logger().info(
                "%s file '%s' from '%s' to '%s'.",
                verb,
                fname,
                url,
                str(self.abspath),
            )

            if downloader is None:
                downloader = choose_downloader(url, progressbar=progressbar)

            stream_download(
                url,
                full_path,
                known_hash,
                downloader,
                pooch=self,
                retry_if_failed=self.retry_if_failed,
            )

        if processor is not None:
            return processor(str(full_path), action, self)

        return str(full_path)

    def _assert_file_in_registry(self, fname):
        """
        Check if a file is in the registry and raise :class:`ValueError` if
        it's not.
        """
        if fname not in self.registry:
            raise ValueError(f"File '{fname}' is not in the registry.")

    def get_url(self, fname):
        """
        Get the full URL to download a file in the registry.

        Parameters
        ----------
        fname : str
            The file name (relative to the *base_url* of the remote data
            storage) to fetch from the local storage.

        """
        self._assert_file_in_registry(fname)
        return self.urls.get(fname, "".join([self.base_url, fname]))

    def load_registry(self, fname):
        """
        Load entries from a file and add them to the registry.

        Use this if you are managing many files.

        Each line of the file should have file name and its hash separated by
        a space. Hash can specify checksum algorithm using "alg:hash" format.
        In case no algorithm is provided, SHA256 is used by default.
        Only one file per line is allowed. Custom download URLs for individual
        files can be specified as a third element on the line. Line comments
        can be added and must be prepended with ``#``.

        Parameters
        ----------
        fname : str | fileobj
            Path (or open file object) to the registry file.

        """
        with contextlib.ExitStack() as stack:
            if hasattr(fname, "read"):
                # It's a file object
                fin = fname
            else:
                # It's a file path
                fin = stack.enter_context(open(fname, encoding="utf-8"))

            for linenum, line in enumerate(fin):
                if isinstance(line, bytes):
                    line = line.decode("utf-8")

                line = line.strip()
                # skip line comments
                if line.startswith("#"):
                    continue

                elements = shlex.split(line)
                if not len(elements) in [0, 2, 3]:
                    raise OSError(
                        f"Invalid entry in Pooch registry file '{fname}': "
                        f"expected 2 or 3 elements in line {linenum + 1} but got "
                        f"{len(elements)}. Offending entry: '{line}'"
                    )
                if elements:
                    file_name = elements[0]
                    file_checksum = elements[1]
                    if len(elements) == 3:
                        file_url = elements[2]
                        self.urls[file_name] = file_url
                    self.registry[file_name] = file_checksum.lower()

    def load_registry_from_doi(self):
        """
        Populate the registry using the data repository API

        Fill the registry with all the files available in the data repository,
        along with their hashes. It will make a request to the data repository
        API to retrieve this information. No file is downloaded during this
        process.

        .. important::

            This method is intended to be used only when the ``base_url`` is
            a DOI.
        """

        # Ensure that this is indeed a DOI-based pooch
        downloader = choose_downloader(self.base_url)
        if not isinstance(downloader, DOIDownloader):
            raise ValueError(
                f"Invalid base_url '{self.base_url}': "
                + "Pooch.load_registry_from_doi is only implemented for DOIs"
            )

        # Create a repository instance
        doi = self.base_url.replace("doi:", "")
        repository = doi_to_repository(doi)

        # Call registry population for this repository
        return repository.populate_registry(self)

    def is_available(self, fname, downloader=None):
        """
        Check availability of a remote file without downloading it.

        Use this method when working with large files to check if they are
        available for download.

        Parameters
        ----------
        fname : str
            The file name (relative to the *base_url* of the remote data
            storage).
        downloader : None or callable
            If not None, then a function (or callable object) that will be
            called to check the availability of the file on the server. See
            :ref:`downloaders` for details.

        Returns
        -------
        status : bool
            True if the file is available for download. False otherwise.

        """
        self._assert_file_in_registry(fname)
        url = self.get_url(fname)
        if downloader is None:
            downloader = choose_downloader(url)
        try:
            available = downloader(url, None, self, check_only=True)
        except TypeError as error:
            error_msg = (
                f"Downloader '{str(downloader)}' does not support availability checks."
            )
            raise NotImplementedError(error_msg) from error
        return available


def download_action(path, known_hash):
    """
    Determine the action that is needed to get the file on disk.

    Parameters
    ----------
    path : PathLike
        The path to the file on disk.
    known_hash : str
        A known hash (checksum) of the file. Will be used to verify the
        download or check if an existing file needs to be updated. By default,
        will assume it's a SHA256 hash. To specify a different hashing method,
        prepend the hash with ``algorithm:``, for example
        ``md5:pw9co2iun29juoh`` or ``sha1:092odwhi2ujdp2du2od2odh2wod2``.

    Returns
    -------
    action, verb : str
        The action that must be taken and the English verb (infinitive form of
        *action*) used in the log:
        * ``'download'``: File does not exist locally and must be downloaded.
        * ``'update'``: File exists locally but needs to be updated.
        * ``'fetch'``: File exists locally and only need to inform its path.


    """
    if not path.exists():
        action = "download"
        verb = "Downloading"
    elif not hash_matches(str(path), known_hash):
        action = "update"
        verb = "Updating"
    else:
        action = "fetch"
        verb = "Fetching"
    return action, verb


def stream_download(url, fname, known_hash, downloader, pooch=None, retry_if_failed=0):
    """
    Stream the file and check that its hash matches the known one.

    The file is first downloaded to a temporary file name in the cache folder.
    It will be moved to the desired file name only if the hash matches the
    known hash. Otherwise, the temporary file is deleted.

    If the download fails for either a bad connection or a hash mismatch, we
    will retry the download the specified number of times in case the failure
    was due to a network error.
    """
    # Lazy import requests to speed up import time
    import requests.exceptions  # pylint: disable=C0415

    # Ensure the parent directory exists in case the file is in a subdirectory.
    # Otherwise, move will cause an error.
    if not fname.parent.exists():
        os.makedirs(str(fname.parent))
    download_attempts = 1 + retry_if_failed
    max_wait = 10
    for i in range(download_attempts):
        try:
            # Stream the file to a temporary so that we can safely check its
            # hash before overwriting the original.
            with temporary_file(path=str(fname.parent)) as tmp:
                downloader(url, tmp, pooch)
                hash_matches(tmp, known_hash, strict=True, source=str(fname.name))
                shutil.move(tmp, str(fname))
            break
        except (ValueError, requests.exceptions.RequestException):
            if i == download_attempts - 1:
                raise
            retries_left = download_attempts - (i + 1)
            get_logger().info(
                "Failed to download '%s'. "
                "Will attempt the download again %d more time%s.",
                str(fname.name),
                retries_left,
                "s" if retries_left > 1 else "",
            )
            time.sleep(min(i + 1, max_wait))
