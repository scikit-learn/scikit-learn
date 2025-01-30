# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
# pylint: disable=redefined-outer-name
"""
Test the core class and factory function.
"""
import hashlib
import os
from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from ..core import create, Pooch, retrieve, download_action, stream_download
from ..utils import get_logger, temporary_file, os_cache
from ..hashes import file_hash, hash_matches

# Import the core module so that we can monkeypatch some functions
from .. import core
from ..downloaders import HTTPDownloader, FTPDownloader

from .utils import (
    pooch_test_url,
    data_over_ftp,
    pooch_test_figshare_url,
    pooch_test_zenodo_url,
    pooch_test_zenodo_with_slash_url,
    pooch_test_dataverse_url,
    pooch_test_registry,
    check_tiny_data,
    check_large_data,
    capture_log,
    mirror_directory,
)

DATA_DIR = str(Path(__file__).parent / "data")
REGISTRY = pooch_test_registry()
BASEURL = pooch_test_url()
FIGSHAREURL = pooch_test_figshare_url()
ZENODOURL = pooch_test_zenodo_url()
ZENODOURL_W_SLASH = pooch_test_zenodo_with_slash_url()
DATAVERSEURL = pooch_test_dataverse_url()
REGISTRY_CORRUPTED = {
    # The same data file but I changed the hash manually to a wrong one
    "tiny-data.txt": "098h0894dba14b12085eacb204284b97e362f4f3e5a5807693cc90ef415c1b2d"
}


@pytest.fixture
def data_dir_mirror(tmp_path):
    """
    Mirror the test data folder on a temporary directory. Needed to avoid
    permission errors when pooch is installed on a non-writable path.
    """
    return mirror_directory(DATA_DIR, tmp_path)


@pytest.mark.network
def test_retrieve():
    "Try downloading some data with retrieve"
    with TemporaryDirectory() as local_store:
        data_file = "tiny-data.txt"
        url = BASEURL + data_file
        # Check that the logs say that the file is being downloaded
        with capture_log() as log_file:
            fname = retrieve(url, known_hash=None, path=local_store)
            logs = log_file.getvalue()
            assert logs.split()[0] == "Downloading"
            assert "SHA256 hash of downloaded file:" in logs
            assert REGISTRY[data_file] in logs
        # Check that the downloaded file has the right content
        assert data_file == fname[-len(data_file) :]
        check_tiny_data(fname)
        assert file_hash(fname) == REGISTRY[data_file]
        # Check that no logging happens when not downloading
        with capture_log() as log_file:
            fname = retrieve(url, known_hash=None, path=local_store)
            assert log_file.getvalue() == ""
        with capture_log() as log_file:
            fname = retrieve(url, known_hash=REGISTRY[data_file], path=local_store)
            assert log_file.getvalue() == ""


@pytest.mark.network
def test_retrieve_fname():
    "Try downloading some data with retrieve and setting the file name"
    with TemporaryDirectory() as local_store:
        data_file = "tiny-data.txt"
        url = BASEURL + data_file
        # Check that the logs say that the file is being downloaded
        with capture_log() as log_file:
            fname = retrieve(url, known_hash=None, path=local_store, fname=data_file)
            logs = log_file.getvalue()
            assert logs.split()[0] == "Downloading"
            assert "SHA256 hash of downloaded file:" in logs
            assert REGISTRY[data_file] in logs
        # Check that the downloaded file has the right name and content
        assert data_file == os.path.split(fname)[1]
        check_tiny_data(fname)
        assert file_hash(fname) == REGISTRY[data_file]


@pytest.mark.network
def test_retrieve_default_path():
    "Try downloading some data with retrieve to the default cache location"
    data_file = "tiny-data.txt"
    url = BASEURL + data_file
    expected_location = os_cache("pooch") / data_file
    try:
        # Check that the logs say that the file is being downloaded
        with capture_log() as log_file:
            fname = retrieve(url, known_hash=None, fname=data_file)
            logs = log_file.getvalue()
            assert logs.split()[0] == "Downloading"
            assert str(os_cache("pooch").resolve()) in logs
            assert "SHA256 hash of downloaded file" in logs
            assert REGISTRY[data_file] in logs
        # Check that the downloaded file has the right content
        assert fname == str(expected_location.resolve())
        check_tiny_data(fname)
        assert file_hash(fname) == REGISTRY[data_file]
    finally:
        if os.path.exists(str(expected_location)):
            os.remove(str(expected_location))


def test_pooch_local(data_dir_mirror):
    "Setup a pooch that already has the local data and test the fetch."
    pup = Pooch(path=data_dir_mirror, base_url="some bogus URL", registry=REGISTRY)
    true = str(data_dir_mirror / "tiny-data.txt")
    fname = pup.fetch("tiny-data.txt")
    assert true == fname
    check_tiny_data(fname)


@pytest.mark.network
@pytest.mark.parametrize(
    "url",
    [BASEURL, FIGSHAREURL, ZENODOURL, DATAVERSEURL],
    ids=["https", "figshare", "zenodo", "dataverse"],
)
def test_pooch_custom_url(url):
    "Have pooch download the file from URL that is not base_url"
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        urls = {"tiny-data.txt": url + "tiny-data.txt"}
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url="", registry=REGISTRY, urls=urls)
        # Check that the logs say that the file is being downloaded
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            logs = log_file.getvalue()
            assert logs.split()[0] == "Downloading"
            assert logs.split()[-1] == f"'{path}'."
        check_tiny_data(fname)
        # Check that no logging happens when there are no events
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            assert log_file.getvalue() == ""


@pytest.mark.network
@pytest.mark.parametrize(
    "url",
    [BASEURL, FIGSHAREURL, ZENODOURL, DATAVERSEURL],
    ids=["https", "figshare", "zenodo", "dataverse"],
)
def test_pooch_download(url):
    "Setup a pooch that has no local data and needs to download"
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        true_path = str(path / "tiny-data.txt")
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url=url, registry=REGISTRY)
        # Check that the logs say that the file is being downloaded
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            logs = log_file.getvalue()
            assert logs.split()[0] == "Downloading"
            assert logs.split()[-1] == f"'{path}'."
        # Check that the downloaded file has the right content
        assert true_path == fname
        check_tiny_data(fname)
        assert file_hash(fname) == REGISTRY["tiny-data.txt"]
        # Check that no logging happens when not downloading
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            assert log_file.getvalue() == ""


class FakeHashMatches:  # pylint: disable=too-few-public-methods
    "Create a fake version of hash_matches that fails n times"

    def __init__(self, nfailures):
        self.nfailures = nfailures
        self.failed = 0

    def hash_matches(self, *args, **kwargs):
        "Fail n times before finally passing"
        if self.failed < self.nfailures:
            self.failed += 1
            # Give it an invalid hash to force a failure
            return hash_matches(args[0], "bla", **kwargs)
        return hash_matches(*args, **kwargs)


@pytest.mark.network
def test_pooch_download_retry_off_by_default(monkeypatch):
    "Check that retrying the download is off by default"
    with TemporaryDirectory() as local_store:
        monkeypatch.setattr(core, "hash_matches", FakeHashMatches(3).hash_matches)
        # Setup a pooch without download retrying
        path = Path(local_store)
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY)
        # Make sure it fails with no retries
        with pytest.raises(ValueError) as error:
            with capture_log() as log_file:
                pup.fetch("tiny-data.txt")
        assert "does not match the known hash" in str(error)
        # Check that the log doesn't have the download retry message
        logs = log_file.getvalue().strip().split("\n")
        assert len(logs) == 1
        assert logs[0].startswith("Downloading")
        assert logs[0].endswith(f"'{path}'.")


class FakeSleep:  # pylint: disable=too-few-public-methods
    "Create a fake version of sleep that logs the specified times"

    def __init__(self):
        self.times = []

    def sleep(self, secs):
        "Store the time and doesn't sleep"
        self.times.append(secs)


@pytest.mark.network
def test_pooch_download_retry(monkeypatch):
    "Check that retrying the download works if the hash is different"
    with TemporaryDirectory() as local_store:
        monkeypatch.setattr(core, "hash_matches", FakeHashMatches(11).hash_matches)
        fakesleep = FakeSleep()
        monkeypatch.setattr(core.time, "sleep", fakesleep.sleep)
        # Setup a pooch with download retrying
        path = Path(local_store)
        true_path = str(path / "tiny-data.txt")
        retries = 11
        pup = Pooch(
            path=path, base_url=BASEURL, registry=REGISTRY, retry_if_failed=retries
        )
        # Check that the logs say that the download failed n times
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            logs = log_file.getvalue().strip().split("\n")
            assert len(logs) == 1 + retries
            assert logs[0].startswith("Downloading")
            assert logs[0].endswith(f"'{path}'.")
            for i, line in zip(range(retries, 0, -1), logs[1:]):
                assert "Failed to download" in line
                plural = "s" if i > 1 else ""
                assert f"download again {i} more time{plural}." in line
            # Check that the sleep time increases but stops at 10s
            assert fakesleep.times == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10]
        # Check that the downloaded file has the right content
        assert true_path == fname
        check_tiny_data(fname)
        assert file_hash(fname) == REGISTRY["tiny-data.txt"]


@pytest.mark.network
def test_pooch_download_retry_fails_eventually(monkeypatch):
    "Check that retrying the download fails after the set amount of retries"
    with TemporaryDirectory() as local_store:
        monkeypatch.setattr(core, "hash_matches", FakeHashMatches(3).hash_matches)
        # Setup a pooch with insufficient retry attempts
        path = Path(local_store)
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY, retry_if_failed=1)
        # Make sure it fails with no retries
        with pytest.raises(ValueError) as error:
            # Check that the logs say that the download failed n times
            with capture_log() as log_file:
                pup.fetch("tiny-data.txt")
        logs = log_file.getvalue().strip().split("\n")
        assert len(logs) == 2
        assert logs[0].startswith("Downloading")
        assert logs[0].endswith(f"'{path}'.")
        assert "Failed to download" in logs[1]
        assert "download again 1 more time." in logs[1]
        assert "does not match the known hash" in str(error)


@pytest.mark.network
def test_pooch_logging_level():
    "Setup a pooch and check that no logging happens when the level is raised"
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        urls = {"tiny-data.txt": BASEURL + "tiny-data.txt"}
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url="", registry=REGISTRY, urls=urls)
        # Capture only critical logging events
        with capture_log("CRITICAL") as log_file:
            fname = pup.fetch("tiny-data.txt")
            assert log_file.getvalue() == ""
        check_tiny_data(fname)


@pytest.mark.network
def test_pooch_update():
    "Setup a pooch that already has the local data but the file is outdated"
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        # Create a dummy version of tiny-data.txt that is different from the
        # one in the remote storage
        true_path = str(path / "tiny-data.txt")
        with open(true_path, "w", encoding="utf-8") as fin:
            fin.write("different data")
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY)
        # Check that the logs say that the file is being updated
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            logs = log_file.getvalue()
            assert logs.split()[0] == "Updating"
            assert logs.split()[-1] == f"'{path}'."
        # Check that the updated file has the right content
        assert true_path == fname
        check_tiny_data(fname)
        assert file_hash(fname) == REGISTRY["tiny-data.txt"]
        # Check that no logging happens when not downloading
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt")
            assert log_file.getvalue() == ""


def test_pooch_update_disallowed():
    "Test that disallowing updates works."
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        # Create a dummy version of tiny-data.txt that is different from the
        # one in the remote storage
        true_path = str(path / "tiny-data.txt")
        with open(true_path, "w", encoding="utf-8") as fin:
            fin.write("different data")
        # Setup a pooch in a temp dir
        pup = Pooch(
            path=path,
            base_url=BASEURL,
            registry=REGISTRY,
            allow_updates=False,
        )
        with pytest.raises(ValueError):
            pup.fetch("tiny-data.txt")


def test_pooch_update_disallowed_environment():
    "Test that disallowing updates works through an environment variable."
    variable_name = "MYPROJECT_DISALLOW_UPDATES"
    try:
        os.environ[variable_name] = "False"
        with TemporaryDirectory() as local_store:
            path = Path(local_store)
            # Create a dummy version of tiny-data.txt that is different from
            # the one in the remote storage
            true_path = str(path / "tiny-data.txt")
            with open(true_path, "w", encoding="utf-8") as fin:
                fin.write("different data")
            # Setup a pooch in a temp dir
            pup = create(
                path=path,
                base_url=BASEURL,
                registry=REGISTRY,
                allow_updates=variable_name,
            )
            with pytest.raises(ValueError):
                pup.fetch("tiny-data.txt")
    finally:
        os.environ.pop(variable_name)


def test_pooch_create_base_url_no_trailing_slash():
    """
    Test if pooch.create appends a trailing slash to the base url if missing
    """
    base_url = "https://mybase.url"
    pup = create(base_url=base_url, registry=None, path=DATA_DIR)
    assert pup.base_url == base_url + "/"


@pytest.mark.network
def test_pooch_corrupted(data_dir_mirror):
    "Raise an exception if the file hash doesn't match the registry"
    # Test the case where the file wasn't in the directory
    with TemporaryDirectory() as local_store:
        path = os.path.abspath(local_store)
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY_CORRUPTED)
        with capture_log() as log_file:
            with pytest.raises(ValueError) as error:
                pup.fetch("tiny-data.txt")
            assert "(tiny-data.txt)" in str(error.value)
            logs = log_file.getvalue()
            assert logs.split()[0] == "Downloading"
            assert logs.split()[-1] == f"'{path}'."
    # and the case where the file exists but hash doesn't match
    pup = Pooch(path=data_dir_mirror, base_url=BASEURL, registry=REGISTRY_CORRUPTED)
    with capture_log() as log_file:
        with pytest.raises(ValueError) as error:
            pup.fetch("tiny-data.txt")
        assert "(tiny-data.txt)" in str(error.value)
        logs = log_file.getvalue()
        assert logs.split()[0] == "Updating"
        assert logs.split()[-1] == f"'{data_dir_mirror}'."


def test_pooch_file_not_in_registry():
    "Should raise an exception if the file is not in the registry."
    pup = Pooch(
        path="it shouldn't matter", base_url="this shouldn't either", registry=REGISTRY
    )
    with pytest.raises(ValueError):
        pup.fetch("this-file-does-not-exit.csv")


def test_pooch_load_registry():
    "Loading the registry from a file should work"
    pup = Pooch(path="", base_url="")
    pup.load_registry(os.path.join(DATA_DIR, "registry.txt"))
    assert pup.registry == REGISTRY
    assert pup.registry_files.sort() == list(REGISTRY).sort()


def test_pooch_load_registry_comments():
    "Loading the registry from a file and strip line comments"
    pup = Pooch(path="", base_url="")
    pup.load_registry(os.path.join(DATA_DIR, "registry_comments.txt"))
    assert pup.registry == REGISTRY
    assert pup.registry_files.sort() == list(REGISTRY).sort()


def test_pooch_load_registry_fileobj():
    "Loading the registry from a file object"
    path = os.path.join(DATA_DIR, "registry.txt")

    # Binary mode
    pup = Pooch(path="", base_url="")
    with open(path, "rb") as fin:
        pup.load_registry(fin)
    assert pup.registry == REGISTRY
    assert pup.registry_files.sort() == list(REGISTRY).sort()

    # Text mode
    pup = Pooch(path="", base_url="")
    with open(path, "r", encoding="utf-8") as fin:
        pup.load_registry(fin)
    assert pup.registry == REGISTRY
    assert pup.registry_files.sort() == list(REGISTRY).sort()


def test_pooch_load_registry_custom_url():
    "Load the registry from a file with a custom URL inserted"
    pup = Pooch(path="", base_url="")
    pup.load_registry(os.path.join(DATA_DIR, "registry-custom-url.txt"))
    assert pup.registry == REGISTRY
    assert pup.urls == {"tiny-data.txt": "https://some-site/tiny-data.txt"}


def test_pooch_load_registry_invalid_line():
    "Should raise an exception when a line doesn't have two elements"
    pup = Pooch(path="", base_url="", registry={})
    with pytest.raises(IOError):
        pup.load_registry(os.path.join(DATA_DIR, "registry-invalid.txt"))


def test_pooch_load_registry_with_spaces():
    "Should check that spaces in filenames are allowed in registry files"
    pup = Pooch(path="", base_url="")
    pup.load_registry(os.path.join(DATA_DIR, "registry-spaces.txt"))
    assert "file with spaces.txt" in pup.registry
    assert "other with spaces.txt" in pup.registry


@pytest.mark.network
def test_check_availability():
    "Should correctly check availability of existing and non existing files"
    # Check available remote file
    pup = Pooch(path=DATA_DIR, base_url=BASEURL, registry=REGISTRY)
    assert pup.is_available("tiny-data.txt")
    # Check non available remote file
    pup = Pooch(path=DATA_DIR, base_url=BASEURL + "wrong-url/", registry=REGISTRY)
    assert not pup.is_available("tiny-data.txt")
    # Wrong file name
    registry = {"not-a-real-data-file.txt": "notarealhash"}
    registry.update(REGISTRY)
    pup = Pooch(path=DATA_DIR, base_url=BASEURL, registry=registry)
    assert not pup.is_available("not-a-real-data-file.txt")


def test_check_availability_on_ftp(ftpserver):
    "Should correctly check availability of existing and non existing files"
    with data_over_ftp(ftpserver, "tiny-data.txt") as url:
        # Check available remote file on FTP server
        pup = Pooch(
            path=DATA_DIR,
            base_url=url.replace("tiny-data.txt", ""),
            registry={
                "tiny-data.txt": "baee0894dba14b12085eacb204284b97e362f4f3e5a5807693cc90ef415c1b2d",
                "doesnot_exist.zip": "jdjdjdjdflld",
            },
        )
        downloader = FTPDownloader(port=ftpserver.server_port)
        assert pup.is_available("tiny-data.txt", downloader=downloader)
        # Check non available remote file
        assert not pup.is_available("doesnot_exist.zip", downloader=downloader)


def test_check_availability_invalid_downloader():
    "Should raise an exception if the downloader doesn't support this"

    def downloader(url, output, pooch):  # pylint: disable=unused-argument
        "A downloader that doesn't support check_only"
        return None

    pup = Pooch(path=DATA_DIR, base_url=BASEURL, registry=REGISTRY)
    msg = "does not support availability checks."
    with pytest.raises(NotImplementedError, match=msg):
        pup.is_available("tiny-data.txt", downloader=downloader)


@pytest.mark.network
def test_fetch_with_downloader(capsys):
    "Setup a downloader function for fetch"

    def download(url, output_file, pup):  # pylint: disable=unused-argument
        "Download through HTTP and warn that we're doing it"
        get_logger().info("downloader executed")
        HTTPDownloader()(url, output_file, pup)

    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY)
        # Check that the logs say that the file is being downloaded
        with capture_log() as log_file:
            fname = pup.fetch("large-data.txt", downloader=download)
            logs = log_file.getvalue()
            lines = logs.splitlines()
            assert len(lines) == 2
            assert lines[0].split()[0] == "Downloading"
            assert lines[1] == "downloader executed"
        # Read stderr and make sure no progress bar was printed by default
        assert not capsys.readouterr().err
        # Check that the downloaded file has the right content
        check_large_data(fname)
        # Check that no logging happens when not downloading
        with capture_log() as log_file:
            fname = pup.fetch("large-data.txt")
            assert log_file.getvalue() == ""


def test_invalid_hash_alg(data_dir_mirror):
    "Test an invalid hashing algorithm"
    pup = Pooch(
        path=data_dir_mirror, base_url=BASEURL, registry={"tiny-data.txt": "blah:1234"}
    )
    with pytest.raises(ValueError) as exc:
        pup.fetch("tiny-data.txt")

    assert "'blah'" in str(exc.value)


def test_alternative_hashing_algorithms(data_dir_mirror):
    "Test different hashing algorithms using local data"
    fname = str(data_dir_mirror / "tiny-data.txt")
    check_tiny_data(fname)
    with open(fname, "rb") as fin:
        data = fin.read()
    for alg in ("sha512", "md5"):
        hasher = hashlib.new(alg)
        hasher.update(data)
        registry = {"tiny-data.txt": f"{alg}:{hasher.hexdigest()}"}
        pup = Pooch(path=data_dir_mirror, base_url="some bogus URL", registry=registry)
        assert fname == pup.fetch("tiny-data.txt")
        check_tiny_data(fname)


def test_download_action():
    "Test that the right action is performed based on file existing"
    action, verb = download_action(
        Path("this_file_does_not_exist.txt"), known_hash=None
    )
    assert action == "download"
    assert verb == "Downloading"

    with temporary_file() as tmp:
        action, verb = download_action(Path(tmp), known_hash="not the correct hash")
    assert action == "update"
    assert verb == "Updating"

    with temporary_file() as tmp:
        with open(tmp, "w", encoding="utf-8") as output:
            output.write("some data")
        action, verb = download_action(Path(tmp), known_hash=file_hash(tmp))
    assert action == "fetch"
    assert verb == "Fetching"


@pytest.mark.network
@pytest.mark.parametrize("fname", ["tiny-data.txt", "subdir/tiny-data.txt"])
def test_stream_download(fname):
    "Check that downloading a file over HTTP works as expected"
    # Use the data in store/ because the subdir is in there for some reason
    url = BASEURL + "store/" + fname
    known_hash = REGISTRY[fname]
    downloader = HTTPDownloader()
    with TemporaryDirectory() as local_store:
        destination = Path(local_store) / fname
        assert not destination.exists()
        stream_download(url, destination, known_hash, downloader, pooch=None)
        assert destination.exists()
        check_tiny_data(str(destination))


@pytest.mark.network
@pytest.mark.parametrize(
    "url",
    [FIGSHAREURL, ZENODOURL, DATAVERSEURL],
    ids=["figshare", "zenodo", "dataverse"],
)
def test_load_registry_from_doi(url):
    """Check that the registry is correctly populated from the API"""
    with TemporaryDirectory() as local_store:
        path = os.path.abspath(local_store)
        pup = Pooch(path=path, base_url=url)
        pup.load_registry_from_doi()

        # Check the existence of all files in the registry
        assert len(pup.registry) == 2
        assert "tiny-data.txt" in pup.registry
        assert "store.zip" in pup.registry

        # Ensure that all files have correct checksums by fetching them
        for filename in pup.registry:
            pup.fetch(filename)


@pytest.mark.network
def test_load_registry_from_doi_zenodo_with_slash():
    """
    Check that the registry is correctly populated from the Zenodo API when
    the filename contains a slash
    """
    url = ZENODOURL_W_SLASH
    with TemporaryDirectory() as local_store:
        path = os.path.abspath(local_store)
        pup = Pooch(path=path, base_url=url)
        pup.load_registry_from_doi()

        # Check the existence of all files in the registry
        assert len(pup.registry) == 1
        assert "santisoler/pooch-test-data-v1.zip" in pup.registry

        # Ensure that all files have correct checksums by fetching them
        for filename in pup.registry:
            pup.fetch(filename)


def test_wrong_load_registry_from_doi():
    """Check that non-DOI URLs produce an error"""

    pup = Pooch(path="", base_url=BASEURL)

    with pytest.raises(ValueError) as exc:
        pup.load_registry_from_doi()

    assert "only implemented for DOIs" in str(exc.value)
