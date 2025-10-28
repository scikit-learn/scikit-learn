# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
The classes that actually handle the downloads.
"""
import os
import sys
import ftplib

import warnings

from .utils import parse_url

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None

try:
    import paramiko
except ImportError:
    paramiko = None


# Set the default timeout in seconds so it can be configured in a pinch for the
# methods that don't or can't expose a way set it at runtime.
# See https://github.com/fatiando/pooch/issues/409
DEFAULT_TIMEOUT = 30


def choose_downloader(url, progressbar=False):
    """
    Choose the appropriate downloader for the given URL based on the protocol.

    Parameters
    ----------
    url : str
        A URL (including protocol).
    progressbar : bool or an arbitrary progress bar object
        If True, will print a progress bar of the download to standard error
        (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to be
        installed. Alternatively, an arbitrary progress bar object can be
        passed. See :ref:`custom-progressbar` for details.

    Returns
    -------
    downloader
        A downloader class, like :class:`pooch.HTTPDownloader`,
        :class:`pooch.FTPDownloader`, or :class: `pooch.SFTPDownloader`.

    Examples
    --------

    >>> downloader = choose_downloader("http://something.com")
    >>> print(downloader.__class__.__name__)
    HTTPDownloader
    >>> downloader = choose_downloader("https://something.com")
    >>> print(downloader.__class__.__name__)
    HTTPDownloader
    >>> downloader = choose_downloader("ftp://something.com")
    >>> print(downloader.__class__.__name__)
    FTPDownloader
    >>> downloader = choose_downloader("doi:DOI/filename.csv")
    >>> print(downloader.__class__.__name__)
    DOIDownloader

    """
    known_downloaders = {
        "ftp": FTPDownloader,
        "https": HTTPDownloader,
        "http": HTTPDownloader,
        "sftp": SFTPDownloader,
        "doi": DOIDownloader,
    }

    parsed_url = parse_url(url)
    if parsed_url["protocol"] not in known_downloaders:
        raise ValueError(
            f"Unrecognized URL protocol '{parsed_url['protocol']}' in '{url}'. "
            f"Must be one of {known_downloaders.keys()}."
        )
    downloader = known_downloaders[parsed_url["protocol"]](progressbar=progressbar)
    return downloader


class HTTPDownloader:  # pylint: disable=too-few-public-methods
    """
    Download manager for fetching files over HTTP/HTTPS.

    When called, downloads the given file URL into the specified local file.
    Uses the :mod:`requests` library to manage downloads.

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to customize
    the download of files (for example, to use authentication or print a
    progress bar).

    Parameters
    ----------
    progressbar : bool or an arbitrary progress bar object
        If True, will print a progress bar of the download to standard error
        (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to be
        installed. Alternatively, an arbitrary progress bar object can be
        passed. See :ref:`custom-progressbar` for details.
    chunk_size : int
        Files are streamed *chunk_size* bytes at a time instead of loading
        everything into memory at one. Usually doesn't need to be changed.
    **kwargs
        All keyword arguments given when creating an instance of this class
        will be passed to :func:`requests.get`.

    Examples
    --------

    Download one of the data files from the Pooch repository:

    >>> import os
    >>> from pooch import __version__, check_version
    >>> url = "https://github.com/fatiando/pooch/raw/{}/data/tiny-data.txt"
    >>> url = url.format(check_version(__version__, fallback="main"))
    >>> downloader = HTTPDownloader()
    >>> # Not using with Pooch.fetch so no need to pass an instance of Pooch
    >>> downloader(url=url, output_file="tiny-data.txt", pooch=None)
    >>> os.path.exists("tiny-data.txt")
    True
    >>> with open("tiny-data.txt") as f:
    ...     print(f.read().strip())
    # A tiny data file for test purposes only
    1  2  3  4  5  6
    >>> os.remove("tiny-data.txt")

    Authentication can be handled by passing a user name and password to
    :func:`requests.get`. All arguments provided when creating an instance of
    the class are forwarded to :func:`requests.get`. We'll use
    ``auth=(username, password)`` to use basic HTTPS authentication. The
    https://httpbin.org website allows us to make a fake a login request using
    whatever username and password we provide to it:

    >>> user = "doggo"
    >>> password = "goodboy"
    >>> # httpbin will ask for the user and password we provide in the URL
    >>> url = f"https://httpbin.org/basic-auth/{user}/{password}"
    >>> # Trying without the login credentials causes an error
    >>> downloader = HTTPDownloader()
    >>> try:
    ...     downloader(url=url, output_file="tiny-data.txt", pooch=None)
    ... except Exception:
    ...     print("There was an error!")
    There was an error!
    >>> # Pass in the credentials to HTTPDownloader
    >>> downloader = HTTPDownloader(auth=(user, password))
    >>> downloader(url=url, output_file="tiny-data.txt", pooch=None)
    >>> with open("tiny-data.txt") as f:
    ...     for line in f:
    ...         print(line.rstrip())
    {
      "authenticated": true,
      "user": "doggo"
    }
    >>> os.remove("tiny-data.txt")

    """

    def __init__(self, progressbar=False, chunk_size=1024, **kwargs):
        self.kwargs = kwargs
        self.progressbar = progressbar
        self.chunk_size = chunk_size
        if self.progressbar is True and tqdm is None:
            raise ValueError("Missing package 'tqdm' required for progress bars.")

    def __call__(
        self, url, output_file, pooch, check_only=False
    ):  # pylint: disable=R0914
        """
        Download the given URL over HTTP to the given output file.

        Uses :func:`requests.get`.

        Parameters
        ----------
        url : str
            The URL to the file you want to download.
        output_file : str or file-like object
            Path (and file name) to which the file will be downloaded.
        pooch : :class:`~pooch.Pooch`
            The instance of :class:`~pooch.Pooch` that is calling this method.
        check_only : bool
            If True, will only check if a file exists on the server and
            **without downloading the file**. Will return ``True`` if the file
            exists and ``False`` otherwise.

        Returns
        -------
        availability : bool or None
            If ``check_only==True``, returns a boolean indicating if the file
            is available on the server. Otherwise, returns ``None``.

        """
        # Lazy import requests to speed up import time
        import requests  # pylint: disable=C0415

        if check_only:
            timeout = self.kwargs.get("timeout", DEFAULT_TIMEOUT)
            response = requests.head(url, timeout=timeout, allow_redirects=True)
            available = bool(response.status_code == 200)
            return available

        kwargs = self.kwargs.copy()
        timeout = kwargs.pop("timeout", DEFAULT_TIMEOUT)
        kwargs.setdefault("stream", True)
        ispath = not hasattr(output_file, "write")
        if ispath:
            # pylint: disable=consider-using-with
            output_file = open(output_file, "w+b")
            # pylint: enable=consider-using-with
        try:
            response = requests.get(url, timeout=timeout, **kwargs)
            response.raise_for_status()
            content = response.iter_content(chunk_size=self.chunk_size)
            total = int(response.headers.get("content-length", 0))
            if self.progressbar is True:
                # Need to use ascii characters on Windows because there isn't
                # always full unicode support
                # (see https://github.com/tqdm/tqdm/issues/454)
                use_ascii = bool(sys.platform == "win32")
                progress = tqdm(
                    total=total,
                    ncols=79,
                    ascii=use_ascii,
                    unit="B",
                    unit_scale=True,
                    leave=True,
                )
            elif self.progressbar:
                progress = self.progressbar
                progress.total = total
            for chunk in content:
                if chunk:
                    output_file.write(chunk)
                    output_file.flush()
                    if self.progressbar:
                        # Use the chunk size here because chunk may be much
                        # larger if the data are decompressed by requests after
                        # reading (happens with text files).
                        progress.update(self.chunk_size)
            # Make sure the progress bar gets filled even if the actual number
            # is chunks is smaller than expected. This happens when streaming
            # text files that are compressed by the server when sending (gzip).
            # Binary files don't experience this.
            if self.progressbar:
                progress.reset()
                progress.update(total)
                progress.close()
        finally:
            if ispath:
                output_file.close()
        return None


class FTPDownloader:  # pylint: disable=too-few-public-methods
    """
    Download manager for fetching files over FTP.

    When called, downloads the given file URL into the specified local file.
    Uses the :mod:`ftplib` module to manage downloads.

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to customize
    the download of files (for example, to use authentication or print a
    progress bar).

    Parameters
    ----------
    port : int
        Port used for the FTP connection.
    username : str
        User name used to login to the server. Only needed if the server
        requires authentication (i.e., no anonymous FTP).
    password : str
        Password used to login to the server. Only needed if the server
        requires authentication (i.e., no anonymous FTP). Use the empty string
        to indicate no password is required.
    account : str
        Some servers also require an "account" name for authentication.
    timeout : int
        Timeout in seconds for ftp socket operations, use None to mean no
        timeout.
    progressbar : bool
        If True, will print a progress bar of the download to standard error
        (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to be
        installed. **Custom progress bars are not yet supported.**
    chunk_size : int
        Files are streamed *chunk_size* bytes at a time instead of loading
        everything into memory at one. Usually doesn't need to be changed.

    """

    def __init__(
        self,
        port=21,
        username="anonymous",
        password="",
        account="",
        timeout=None,
        progressbar=False,
        chunk_size=1024,
    ):
        self.port = port
        self.username = username
        self.password = password
        self.account = account
        self.timeout = timeout
        self.progressbar = progressbar
        self.chunk_size = chunk_size
        if self.progressbar is True and tqdm is None:
            raise ValueError("Missing package 'tqdm' required for progress bars.")

    def __call__(self, url, output_file, pooch, check_only=False):
        """
        Download the given URL over FTP to the given output file.

        Parameters
        ----------
        url : str
            The URL to the file you want to download.
        output_file : str or file-like object
            Path (and file name) to which the file will be downloaded.
        pooch : :class:`~pooch.Pooch`
            The instance of :class:`~pooch.Pooch` that is calling this method.
        check_only : bool
            If True, will only check if a file exists on the server and
            **without downloading the file**. Will return ``True`` if the file
            exists and ``False`` otherwise.

        Returns
        -------
        availability : bool or None
            If ``check_only==True``, returns a boolean indicating if the file
            is available on the server. Otherwise, returns ``None``.

        """
        parsed_url = parse_url(url)
        ftp = ftplib.FTP(timeout=self.timeout)
        ftp.connect(host=parsed_url["netloc"], port=self.port)

        if check_only:
            directory, file_name = os.path.split(parsed_url["path"])
            try:
                ftp.login(user=self.username, passwd=self.password, acct=self.account)
                available = file_name in ftp.nlst(directory)
            finally:
                ftp.close()
            return available

        ispath = not hasattr(output_file, "write")
        if ispath:
            # pylint: disable=consider-using-with
            output_file = open(output_file, "w+b")
            # pylint: enable=consider-using-with
        try:
            ftp.login(user=self.username, passwd=self.password, acct=self.account)
            command = f"RETR {parsed_url['path']}"
            if self.progressbar:
                # Make sure the file is set to binary mode, otherwise we can't
                # get the file size. See: https://stackoverflow.com/a/22093848
                ftp.voidcmd("TYPE I")
                use_ascii = bool(sys.platform == "win32")
                progress = tqdm(
                    total=int(ftp.size(parsed_url["path"])),
                    ncols=79,
                    ascii=use_ascii,
                    unit="B",
                    unit_scale=True,
                    leave=True,
                )
                with progress:

                    def callback(data):
                        "Update the progress bar and write to output"
                        progress.update(len(data))
                        output_file.write(data)

                    ftp.retrbinary(command, callback, blocksize=self.chunk_size)
            else:
                ftp.retrbinary(command, output_file.write, blocksize=self.chunk_size)
        finally:
            ftp.quit()
            if ispath:
                output_file.close()
        return None


class SFTPDownloader:  # pylint: disable=too-few-public-methods
    """
    Download manager for fetching files over SFTP.

    When called, downloads the given file URL into the specified local file.
    Requires `paramiko <https://github.com/paramiko/paramiko>`__ to be
    installed.

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to customize
    the download of files (for example, to use authentication or print a
    progress bar).

    Parameters
    ----------
    port : int
        Port used for the SFTP connection.
    username : str
        User name used to login to the server. Only needed if the server
        requires authentication (i.e., no anonymous SFTP).
    password : str
        Password used to login to the server. Only needed if the server
        requires authentication (i.e., no anonymous SFTP). Use the empty
        string to indicate no password is required.
    timeout : int
        Timeout in seconds for sftp socket operations, use None to mean no
        timeout.
    progressbar : bool or an arbitrary progress bar object
        If True, will print a progress bar of the download to standard
        error (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to
        be installed.

    """

    def __init__(
        self,
        port=22,
        username="anonymous",
        password="",
        account="",
        timeout=None,
        progressbar=False,
    ):
        self.port = port
        self.username = username
        self.password = password
        self.account = account
        self.timeout = timeout
        self.progressbar = progressbar
        # Collect errors and raise only once so that both missing packages are
        # captured. Otherwise, the user is only warned of one of them at a
        # time (and we can't test properly when they are both missing).
        errors = []
        if self.progressbar and tqdm is None:
            errors.append("Missing package 'tqdm' required for progress bars.")
        if paramiko is None:
            errors.append("Missing package 'paramiko' required for SFTP downloads.")
        if errors:
            raise ValueError(" ".join(errors))

    def __call__(self, url, output_file, pooch):
        """
        Download the given URL over SFTP to the given output file.

        The output file must be given as a string (file name/path) and not an
        open file object! Otherwise, paramiko cannot save to that file.

        Parameters
        ----------
        url : str
            The URL to the file you want to download.
        output_file : str
            Path (and file name) to which the file will be downloaded. **Cannot
            be a file object**.
        pooch : :class:`~pooch.Pooch`
            The instance of :class:`~pooch.Pooch` that is calling this method.
        """
        parsed_url = parse_url(url)
        connection = paramiko.Transport(sock=(parsed_url["netloc"], self.port))
        sftp = None
        try:
            connection.connect(username=self.username, password=self.password)
            sftp = paramiko.SFTPClient.from_transport(connection)
            sftp.get_channel().settimeout = self.timeout
            if self.progressbar:
                size = int(sftp.stat(parsed_url["path"]).st_size)
                use_ascii = bool(sys.platform == "win32")
                progress = tqdm(
                    total=size,
                    ncols=79,
                    ascii=use_ascii,
                    unit="B",
                    unit_scale=True,
                    leave=True,
                )
            if self.progressbar:
                with progress:

                    def callback(current, total):
                        "Update the progress bar and write to output"
                        progress.total = int(total)
                        progress.update(int(current - progress.n))

                    sftp.get(parsed_url["path"], output_file, callback=callback)
            else:
                sftp.get(parsed_url["path"], output_file)
        finally:
            connection.close()
            if sftp is not None:
                sftp.close()


class DOIDownloader:  # pylint: disable=too-few-public-methods
    """
    Download manager for fetching files from Digital Object Identifiers (DOIs).

    Open-access data repositories often issue Digital Object Identifiers (DOIs)
    for data which provide a stable link and citation point. The trick is
    finding out the download URL for a file given the DOI.

    When called, this downloader uses the repository's public API to find out
    the download URL from the DOI and file name. It then uses
    :class:`pooch.HTTPDownloader` to download the URL into the specified local
    file. Allowing "URL"s  to be specified with the DOI instead of the actual
    HTTP download link. Uses the :mod:`requests` library to manage downloads
    and interact with the APIs.

    The **format of the "URL"** is: ``doi:{DOI}/{file name}``.

    Notice that there are no ``//`` like in HTTP/FTP and you must specify a
    file name after the DOI (separated by a ``/``).

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to be able to
    download files given the DOI instead of an HTTP link.

    Supported repositories:

    * `figshare <https://www.figshare.com>`__
    * `Zenodo <https://www.zenodo.org>`__
    * `Dataverse <https://dataverse.org/>`__ instances

    .. attention::

        DOIs from other repositories **will not work** since we need to access
        their particular APIs to find the download links. We welcome
        suggestions and contributions adding new repositories.

    Parameters
    ----------
    progressbar : bool or an arbitrary progress bar object
        If True, will print a progress bar of the download to standard error
        (stderr). Requires `tqdm <https://github.com/tqdm/tqdm>`__ to be
        installed. Alternatively, an arbitrary progress bar object can be
        passed. See :ref:`custom-progressbar` for details.
    chunk_size : int
        Files are streamed *chunk_size* bytes at a time instead of loading
        everything into memory at one. Usually doesn't need to be changed.
    **kwargs
        All keyword arguments given when creating an instance of this class
        will be passed to :func:`requests.get`.

    Examples
    --------

    Download one of the data files from the figshare archive of Pooch test
    data:

    >>> import os
    >>> downloader = DOIDownloader()
    >>> url = "doi:10.6084/m9.figshare.14763051.v1/tiny-data.txt"
    >>> # Not using with Pooch.fetch so no need to pass an instance of Pooch
    >>> downloader(url=url, output_file="tiny-data.txt", pooch=None)
    >>> os.path.exists("tiny-data.txt")
    True
    >>> with open("tiny-data.txt") as f:
    ...     print(f.read().strip())
    # A tiny data file for test purposes only
    1  2  3  4  5  6
    >>> os.remove("tiny-data.txt")

    Same thing but for our Zenodo archive:

    >>> url = "doi:10.5281/zenodo.4924875/tiny-data.txt"
    >>> downloader(url=url, output_file="tiny-data.txt", pooch=None)
    >>> os.path.exists("tiny-data.txt")
    True
    >>> with open("tiny-data.txt") as f:
    ...     print(f.read().strip())
    # A tiny data file for test purposes only
    1  2  3  4  5  6
    >>> os.remove("tiny-data.txt")

    """

    def __init__(self, progressbar=False, chunk_size=1024, **kwargs):
        self.kwargs = kwargs
        self.progressbar = progressbar
        self.chunk_size = chunk_size

    def __call__(self, url, output_file, pooch):
        """
        Download the given DOI URL over HTTP to the given output file.

        Uses the repository's API to determine the actual HTTP download URL
        from the given DOI.

        Uses :func:`requests.get`.

        Parameters
        ----------
        url : str
            The URL to the file you want to download.
        output_file : str or file-like object
            Path (and file name) to which the file will be downloaded.
        pooch : :class:`~pooch.Pooch`
            The instance of :class:`~pooch.Pooch` that is calling this method.

        """

        parsed_url = parse_url(url)
        data_repository = doi_to_repository(parsed_url["netloc"])

        # Resolve the URL
        file_name = parsed_url["path"]
        # remove the leading slash in the path
        if file_name[0] == "/":
            file_name = file_name[1:]
        download_url = data_repository.download_url(file_name)

        # Instantiate the downloader object
        downloader = HTTPDownloader(
            progressbar=self.progressbar, chunk_size=self.chunk_size, **self.kwargs
        )
        downloader(download_url, output_file, pooch)


def doi_to_url(doi):
    """
    Follow a DOI link to resolve the URL of the archive.

    Parameters
    ----------
    doi : str
        The DOI of the archive.

    Returns
    -------
    url : str
        The URL of the archive in the data repository.

    """
    # Lazy import requests to speed up import time
    import requests  # pylint: disable=C0415

    # Use doi.org to resolve the DOI to the repository website.
    response = requests.get(f"https://doi.org/{doi}", timeout=DEFAULT_TIMEOUT)
    url = response.url
    if 400 <= response.status_code < 600:
        raise ValueError(
            f"Archive with doi:{doi} not found (see {url}). Is the DOI correct?"
        )
    return url


def doi_to_repository(doi):
    """
    Instantiate a data repository instance from a given DOI.

    This function implements the chain of responsibility dispatch
    to the correct data repository class.

    Parameters
    ----------
    doi : str
        The DOI of the archive.

    Returns
    -------
    data_repository : DataRepository
        The data repository object
    """

    # This should go away in a separate issue: DOI handling should
    # not rely on the (non-)existence of trailing slashes. The issue
    # is documented in https://github.com/fatiando/pooch/issues/324
    if doi[-1] == "/":
        doi = doi[:-1]

    repositories = [
        FigshareRepository,
        ZenodoRepository,
        DataverseRepository,
    ]

    # Extract the DOI and the repository information
    archive_url = doi_to_url(doi)

    # Try the converters one by one until one of them returned a URL
    data_repository = None
    for repo in repositories:
        if data_repository is None:
            data_repository = repo.initialize(
                archive_url=archive_url,
                doi=doi,
            )

    if data_repository is None:
        repository = parse_url(archive_url)["netloc"]
        raise ValueError(
            f"Invalid data repository '{repository}'. "
            "To request or contribute support for this repository, "
            "please open an issue at https://github.com/fatiando/pooch/issues"
        )

    return data_repository


class DataRepository:  # pylint: disable=too-few-public-methods, missing-class-docstring
    @classmethod
    def initialize(cls, doi, archive_url):  # pylint: disable=unused-argument
        """
        Initialize the data repository if the given URL points to a
        corresponding repository.

        Initializes a data repository object. This is done as part of
        a chain of responsibility. If the class cannot handle the given
        repository URL, it returns `None`. Otherwise a `DataRepository`
        instance is returned.

        Parameters
        ----------
        doi : str
            The DOI that identifies the repository
        archive_url : str
            The resolved URL for the DOI
        """

        return None  # pragma: no cover

    def download_url(self, file_name):
        """
        Use the repository API to get the download URL for a file given
        the archive URL.

        Parameters
        ----------
        file_name : str
            The name of the file in the archive that will be downloaded.

        Returns
        -------
        download_url : str
            The HTTP URL that can be used to download the file.
        """

        raise NotImplementedError  # pragma: no cover

    def populate_registry(self, pooch):
        """
        Populate the registry using the data repository's API

        Parameters
        ----------
        pooch : Pooch
            The pooch instance that the registry will be added to.
        """

        raise NotImplementedError  # pragma: no cover


class ZenodoRepository(DataRepository):  # pylint: disable=missing-class-docstring
    base_api_url = "https://zenodo.org/api/records"

    def __init__(self, doi, archive_url):
        self.archive_url = archive_url
        self.doi = doi
        self._api_response = None
        self._api_version = None

    @classmethod
    def initialize(cls, doi, archive_url):
        """
        Initialize the data repository if the given URL points to a
        corresponding repository.

        Initializes a data repository object. This is done as part of
        a chain of responsibility. If the class cannot handle the given
        repository URL, it returns `None`. Otherwise a `DataRepository`
        instance is returned.

        Parameters
        ----------
        doi : str
            The DOI that identifies the repository
        archive_url : str
            The resolved URL for the DOI
        """

        # Check whether this is a Zenodo URL
        parsed_archive_url = parse_url(archive_url)
        if parsed_archive_url["netloc"] != "zenodo.org":
            return None

        return cls(doi, archive_url)

    @property
    def api_response(self):
        """Cached API response from Zenodo"""
        if self._api_response is None:
            # Lazy import requests to speed up import time
            import requests  # pylint: disable=C0415

            article_id = self.archive_url.split("/")[-1]
            self._api_response = requests.get(
                f"{self.base_api_url}/{article_id}",
                timeout=DEFAULT_TIMEOUT,
            ).json()

        return self._api_response

    @property
    def api_version(self):
        """
        Version of the Zenodo API we are interacting with

        The versions can either be :

        - ``"legacy"``: corresponds to the Zenodo API that was supported until
          2023-10-12 (before the migration to InvenioRDM).
        - ``"new"``: corresponds to the new API that went online on 2023-10-13
          after the migration to InvenioRDM.

        The ``"new"`` API breaks backward compatibility with the ``"legacy"``
        one and could probably be replaced by an updated version that restores
        the behaviour of the ``"legacy"`` one.

        Returns
        -------
        str
        """
        if self._api_version is None:
            if all("key" in file for file in self.api_response["files"]):
                self._api_version = "legacy"
            elif all("filename" in file for file in self.api_response["files"]):
                self._api_version = "new"
            else:
                raise ValueError(
                    "Couldn't determine the version of the Zenodo API for "
                    f"{self.archive_url} (doi:{self.doi})."
                )
        return self._api_version

    def download_url(self, file_name):
        """
        Use the repository API to get the download URL for a file given
        the archive URL.

        Parameters
        ----------
        file_name : str
            The name of the file in the archive that will be downloaded.

        Returns
        -------
        download_url : str
            The HTTP URL that can be used to download the file.

        Notes
        -----
        After Zenodo migrated to InvenioRDM on Oct 2023, their API changed. The
        link to the desired files that appears in the API response leads to 404
        errors (by 2023-10-17). The files are available in the following url:
        ``https://zenodo.org/records/{article_id}/files/{file_name}?download=1``.

        This method supports both the legacy and the new API.
        """
        # Create list of files in the repository
        if self.api_version == "legacy":
            files = {item["key"]: item for item in self.api_response["files"]}
        else:
            files = [item["filename"] for item in self.api_response["files"]]
        # Check if file exists in the repository
        if file_name not in files:
            raise ValueError(
                f"File '{file_name}' not found in data archive "
                f"{self.archive_url} (doi:{self.doi})."
            )
        # Build download url
        if self.api_version == "legacy":
            download_url = files[file_name]["links"]["self"]
        else:
            article_id = self.api_response["id"]
            download_url = (
                f"https://zenodo.org/records/{article_id}/files/{file_name}?download=1"
            )
        return download_url

    def populate_registry(self, pooch):
        """
        Populate the registry using the data repository's API

        Parameters
        ----------
        pooch : Pooch
            The pooch instance that the registry will be added to.

        Notes
        -----
        After Zenodo migrated to InvenioRDM on Oct 2023, their API changed. The
        checksums for each file listed in the API reference is now an md5 sum.

        This method supports both the legacy and the new API.
        """
        for filedata in self.api_response["files"]:
            checksum = filedata["checksum"]
            if self.api_version == "legacy":
                key = "key"
            else:
                key = "filename"
                checksum = f"md5:{checksum}"
            pooch.registry[filedata[key]] = checksum


class FigshareRepository(DataRepository):  # pylint: disable=missing-class-docstring
    def __init__(self, doi, archive_url):
        self.archive_url = archive_url
        self.doi = doi
        self._api_response = None

    @classmethod
    def initialize(cls, doi, archive_url):
        """
        Initialize the data repository if the given URL points to a
        corresponding repository.

        Initializes a data repository object. This is done as part of
        a chain of responsibility. If the class cannot handle the given
        repository URL, it returns `None`. Otherwise a `DataRepository`
        instance is returned.

        Parameters
        ----------
        doi : str
            The DOI that identifies the repository
        archive_url : str
            The resolved URL for the DOI
        """

        # Check whether this is a Figshare URL
        parsed_archive_url = parse_url(archive_url)
        if parsed_archive_url["netloc"] != "figshare.com":
            return None

        return cls(doi, archive_url)

    def _parse_version_from_doi(self):
        """
        Parse version from the doi

        Return None if version is not available in the doi.
        """
        # Get suffix of the doi
        _, suffix = self.doi.split("/")
        # Split the suffix by dots and keep the last part
        last_part = suffix.split(".")[-1]
        # Parse the version from the last part
        if last_part[0] != "v":
            return None
        version = int(last_part[1:])
        return version

    @property
    def api_response(self):
        """Cached API response from Figshare"""
        if self._api_response is None:
            # Lazy import requests to speed up import time
            import requests  # pylint: disable=C0415

            # Use the figshare API to find the article ID from the DOI
            article = requests.get(
                f"https://api.figshare.com/v2/articles?doi={self.doi}",
                timeout=DEFAULT_TIMEOUT,
            ).json()[0]
            article_id = article["id"]
            # Parse desired version from the doi
            version = self._parse_version_from_doi()
            # With the ID and version, we can get a list of files and their
            # download links
            if version is None:
                # Figshare returns the latest version available when no version
                # is specified through the DOI.
                warnings.warn(
                    f"The Figshare DOI '{self.doi}' doesn't specify which version of "
                    "the repository should be used. "
                    "Figshare will point to the latest version available.",
                    UserWarning,
                )
                # Define API url using only the article id
                # (figshare will resolve the latest version)
                api_url = f"https://api.figshare.com/v2/articles/{article_id}"
            else:
                # Define API url using article id and the desired version
                # Get list of files using article id and the version
                api_url = (
                    "https://api.figshare.com/v2/articles/"
                    f"{article_id}/versions/{version}"
                )
            # Make the request and return the files in the figshare repository
            response = requests.get(api_url, timeout=DEFAULT_TIMEOUT)
            response.raise_for_status()
            self._api_response = response.json()["files"]

        return self._api_response

    def download_url(self, file_name):
        """
        Use the repository API to get the download URL for a file given
        the archive URL.

        Parameters
        ----------
        file_name : str
            The name of the file in the archive that will be downloaded.

        Returns
        -------
        download_url : str
            The HTTP URL that can be used to download the file.
        """
        files = {item["name"]: item for item in self.api_response}
        if file_name not in files:
            raise ValueError(
                f"File '{file_name}' not found in data archive {self.archive_url} (doi:{self.doi})."
            )
        download_url = files[file_name]["download_url"]
        return download_url

    def populate_registry(self, pooch):
        """
        Populate the registry using the data repository's API

        Parameters
        ----------
        pooch : Pooch
            The pooch instance that the registry will be added to.
        """

        for filedata in self.api_response:
            pooch.registry[filedata["name"]] = f"md5:{filedata['computed_md5']}"


class DataverseRepository(DataRepository):  # pylint: disable=missing-class-docstring
    def __init__(self, doi, archive_url):
        self.archive_url = archive_url
        self.doi = doi
        self._api_response = None

    @classmethod
    def initialize(cls, doi, archive_url):
        """
        Initialize the data repository if the given URL points to a
        corresponding repository.

        Initializes a data repository object. This is done as part of
        a chain of responsibility. If the class cannot handle the given
        repository URL, it returns `None`. Otherwise a `DataRepository`
        instance is returned.

        Parameters
        ----------
        doi : str
            The DOI that identifies the repository
        archive_url : str
            The resolved URL for the DOI
        """
        # Access the DOI as if this was a DataVerse instance
        response = cls._get_api_response(doi, archive_url)

        # If we failed, this is probably not a DataVerse instance
        if 400 <= response.status_code < 600:
            return None

        # Initialize the repository and overwrite the api response
        repository = cls(doi, archive_url)
        repository.api_response = response
        return repository

    @classmethod
    def _get_api_response(cls, doi, archive_url):
        """
        Perform the actual API request

        This has been separated into a separate ``classmethod``, as it can be
        used prior and after the initialization.
        """
        # Lazy import requests to speed up import time
        import requests  # pylint: disable=C0415

        parsed = parse_url(archive_url)
        response = requests.get(
            f"{parsed['protocol']}://{parsed['netloc']}/api/datasets/"
            f":persistentId?persistentId=doi:{doi}",
            timeout=DEFAULT_TIMEOUT,
        )
        return response

    @property
    def api_response(self):
        """Cached API response from a DataVerse instance"""

        if self._api_response is None:
            self._api_response = self._get_api_response(
                self.doi, self.archive_url
            )  # pragma: no cover

        return self._api_response

    @api_response.setter
    def api_response(self, response):
        """Update the cached API response"""

        self._api_response = response

    def download_url(self, file_name):
        """
        Use the repository API to get the download URL for a file given
        the archive URL.

        Parameters
        ----------
        file_name : str
            The name of the file in the archive that will be downloaded.

        Returns
        -------
        download_url : str
            The HTTP URL that can be used to download the file.
        """
        parsed = parse_url(self.archive_url)
        response = self.api_response.json()
        files = {
            file["dataFile"]["filename"]: file["dataFile"]
            for file in response["data"]["latestVersion"]["files"]
        }
        if file_name not in files:
            raise ValueError(
                f"File '{file_name}' not found in data archive "
                f"{self.archive_url} (doi:{self.doi})."
            )
        # Generate download_url using the file id
        download_url = (
            f"{parsed['protocol']}://{parsed['netloc']}/api/access/datafile/"
            f"{files[file_name]['id']}"
        )
        return download_url

    def populate_registry(self, pooch):
        """
        Populate the registry using the data repository's API

        Parameters
        ----------
        pooch : Pooch
            The pooch instance that the registry will be added to.
        """

        for filedata in self.api_response.json()["data"]["latestVersion"]["files"]:
            pooch.registry[filedata["dataFile"]["filename"]] = (
                f"md5:{filedata['dataFile']['md5']}"
            )
