# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
Test the downloader classes and functions separately from the Pooch core.
"""
import os
import sys
from tempfile import TemporaryDirectory

import pytest

try:
    import tqdm
except ImportError:
    tqdm = None

try:
    import paramiko
except ImportError:
    paramiko = None

from .. import Pooch
from ..downloaders import (
    HTTPDownloader,
    FTPDownloader,
    SFTPDownloader,
    DOIDownloader,
    choose_downloader,
    FigshareRepository,
    ZenodoRepository,
    DataverseRepository,
    doi_to_url,
)
from ..processors import Unzip
from .utils import (
    pooch_test_url,
    check_large_data,
    check_tiny_data,
    data_over_ftp,
    pooch_test_figshare_url,
    pooch_test_zenodo_url,
    pooch_test_zenodo_with_slash_url,
    pooch_test_dataverse_url,
)


BASEURL = pooch_test_url()
FIGSHAREURL = pooch_test_figshare_url()
ZENODOURL = pooch_test_zenodo_url()
ZENODOURL_W_SLASH = pooch_test_zenodo_with_slash_url()
DATAVERSEURL = pooch_test_dataverse_url()


@pytest.mark.skipif(tqdm is None, reason="requires tqdm")
@pytest.mark.parametrize(
    "url",
    [
        BASEURL + "tiny-data.txt",  # HTTPDownloader
        FIGSHAREURL,  # DOIDownloader
    ],
)
def test_progressbar_kwarg_passed(url):
    """The progressbar keyword argument must pass through choose_downloader"""
    downloader = choose_downloader(url, progressbar=True)
    assert downloader.progressbar is True


@pytest.mark.skipif(paramiko is None, reason="requires paramiko")
def test_progressbar_kwarg_passed_sftp():
    """The progressbar keyword argument must pass through choose_downloader"""
    url = "sftp://test.rebex.net/pub/example/pocketftp.png"
    downloader = choose_downloader(url, progressbar=True)
    assert downloader.progressbar is True


def test_unsupported_protocol():
    "Should raise ValueError when protocol is not supported"
    with pytest.raises(ValueError):
        choose_downloader("httpup://some-invalid-url.com")
    # Simulate the DOI format
    with pytest.raises(ValueError):
        choose_downloader("doii:XXX/XXX/file")


@pytest.mark.network
def test_invalid_doi_repository():
    "Should fail if data repository is not supported"
    with pytest.raises(ValueError) as exc:
        # Use the DOI of the Pooch paper in JOSS (not a data repository)
        DOIDownloader()(
            url="doi:10.21105/joss.01943/file_name.txt", output_file=None, pooch=None
        )
    assert "Invalid data repository 'joss.theoj.org'" in str(exc.value)


@pytest.mark.network
def test_doi_url_not_found():
    "Should fail if the DOI is not found"
    with pytest.raises(ValueError) as exc:
        doi_to_url(doi="NOTAREALDOI")
    assert "Is the DOI correct?" in str(exc.value)


@pytest.mark.network
@pytest.mark.parametrize(
    "repository,doi",
    [
        (FigshareRepository, "10.6084/m9.figshare.14763051.v1"),
        (ZenodoRepository, "10.5281/zenodo.4924875"),
        (DataverseRepository, "10.11588/data/TKCFEF"),
    ],
    ids=["figshare", "zenodo", "dataverse"],
)
def test_figshare_url_file_not_found(repository, doi):
    "Should fail if the file is not found in the archive"
    with pytest.raises(ValueError) as exc:
        url = doi_to_url(doi)
        repo = repository.initialize(doi, url)
        repo.download_url(file_name="bla.txt")
    assert "File 'bla.txt' not found" in str(exc.value)


@pytest.mark.network
@pytest.mark.parametrize(
    "url",
    [FIGSHAREURL, ZENODOURL, DATAVERSEURL],
    ids=["figshare", "zenodo", "dataverse"],
)
def test_doi_downloader(url):
    "Test the DOI downloader"
    # Use the test data we have on the repository
    with TemporaryDirectory() as local_store:
        downloader = DOIDownloader()
        outfile = os.path.join(local_store, "tiny-data.txt")
        downloader(url + "tiny-data.txt", outfile, None)
        check_tiny_data(outfile)


@pytest.mark.network
def test_zenodo_downloader_with_slash_in_fname():
    """
    Test the Zenodo downloader when the path contains a forward slash

    Related to issue #336
    """
    # Use the test data we have on the repository
    with TemporaryDirectory() as local_store:
        base_url = ZENODOURL_W_SLASH + "santisoler/pooch-test-data-v1.zip"
        downloader = DOIDownloader()
        outfile = os.path.join(local_store, "test-data.zip")
        downloader(base_url, outfile, None)
        # unpack the downloaded zip file so we can check the integrity of
        # tiny-data.txt
        fnames = Unzip()(outfile, action="download", pooch=None)
        (fname,) = [f for f in fnames if "tiny-data.txt" in f]
        check_tiny_data(fname)


@pytest.mark.network
def test_figshare_unspecified_version():
    """
    Test if passing a Figshare url without a version warns about it, but still
    downloads it.
    """
    url = FIGSHAREURL
    # Remove the last bits of the doi, where the version is specified and
    url = url[: url.rindex(".")] + "/"
    # Create expected warning message
    doi = url[4:-1]
    warning_msg = f"The Figshare DOI '{doi}' doesn't specify which version of "
    with TemporaryDirectory() as local_store:
        downloader = DOIDownloader()
        outfile = os.path.join(local_store, "tiny-data.txt")
        with pytest.warns(UserWarning, match=warning_msg):
            downloader(url + "tiny-data.txt", outfile, None)


@pytest.mark.network
@pytest.mark.parametrize(
    "version, missing, present",
    [
        (
            1,
            "LC08_L2SP_218074_20190114_20200829_02_T1-cropped.tar.gz",
            "cropped-before.tar.gz",
        ),
        (
            2,
            "cropped-before.tar.gz",
            "LC08_L2SP_218074_20190114_20200829_02_T1-cropped.tar.gz",
        ),
    ],
)
def test_figshare_data_repository_versions(version, missing, present):
    """
    Test if setting the version in Figshare DOI works as expected
    """
    # Use a Figshare repo as example (we won't download files from it since
    # they are too big)
    doi = f"10.6084/m9.figshare.21665630.v{version}"
    url = f"https://doi.org/{doi}/"
    figshare = FigshareRepository(doi, url)
    filenames = [item["name"] for item in figshare.api_response]
    assert present in filenames
    assert missing not in filenames


@pytest.mark.network
def test_ftp_downloader(ftpserver):
    "Test ftp downloader"
    with data_over_ftp(ftpserver, "tiny-data.txt") as url:
        with TemporaryDirectory() as local_store:
            downloader = FTPDownloader(port=ftpserver.server_port)
            outfile = os.path.join(local_store, "tiny-data.txt")
            downloader(url, outfile, None)
            check_tiny_data(outfile)


@pytest.mark.network
@pytest.mark.skipif(paramiko is None, reason="requires paramiko to run SFTP")
def test_sftp_downloader():
    "Test sftp downloader"
    with TemporaryDirectory() as local_store:
        downloader = SFTPDownloader(username="demo", password="password")
        url = "sftp://test.rebex.net/pub/example/pocketftp.png"
        outfile = os.path.join(local_store, "pocketftp.png")
        downloader(url, outfile, None)
        assert os.path.exists(outfile)


@pytest.mark.network
@pytest.mark.skipif(paramiko is None, reason="requires paramiko to run SFTP")
def test_sftp_downloader_fail_if_file_object():
    "Downloader should fail when a file object rather than string is passed"
    with TemporaryDirectory() as local_store:
        downloader = SFTPDownloader(username="demo", password="password")
        url = "sftp://test.rebex.net/pub/example/pocketftp.png"
        outfile = os.path.join(local_store, "pocketftp.png")
        with open(outfile, "wb") as outfile_obj:
            with pytest.raises(TypeError):
                downloader(url, outfile_obj, None)


@pytest.mark.skipif(paramiko is not None, reason="paramiko must be missing")
def test_sftp_downloader_fail_if_paramiko_missing():
    "test must fail if paramiko is not installed"
    with pytest.raises(ValueError) as exc:
        SFTPDownloader()
    assert "'paramiko'" in str(exc.value)


@pytest.mark.skipif(tqdm is not None, reason="tqdm must be missing")
@pytest.mark.parametrize("downloader", [HTTPDownloader, FTPDownloader, SFTPDownloader])
def test_downloader_progressbar_fails(downloader):
    "Make sure an error is raised if trying to use progressbar without tqdm"
    with pytest.raises(ValueError) as exc:
        downloader(progressbar=True)
    assert "'tqdm'" in str(exc.value)


@pytest.mark.network
@pytest.mark.skipif(tqdm is None, reason="requires tqdm")
@pytest.mark.parametrize(
    "url,downloader",
    [(BASEURL, HTTPDownloader), (FIGSHAREURL, DOIDownloader)],
    ids=["http", "figshare"],
)
def test_downloader_progressbar(url, downloader, capsys):
    "Setup a downloader function that prints a progress bar for fetch"
    download = downloader(progressbar=True)
    with TemporaryDirectory() as local_store:
        fname = "tiny-data.txt"
        url = url + fname
        outfile = os.path.join(local_store, fname)
        download(url, outfile, None)
        # Read stderr and make sure the progress bar is printed only when told
        captured = capsys.readouterr()
        printed = captured.err.split("\r")[-1].strip()
        assert len(printed) == 79
        if sys.platform == "win32":
            progress = "100%|####################"
        else:
            progress = "100%|████████████████████"
        # Bar size is not always the same so can't reliably test the whole bar.
        assert printed[:25] == progress
        # Check that the downloaded file has the right content
        check_tiny_data(outfile)


@pytest.mark.network
@pytest.mark.skipif(tqdm is None, reason="requires tqdm")
def test_downloader_progressbar_ftp(capsys, ftpserver):
    "Setup an FTP downloader function that prints a progress bar for fetch"
    with data_over_ftp(ftpserver, "tiny-data.txt") as url:
        download = FTPDownloader(progressbar=True, port=ftpserver.server_port)
        with TemporaryDirectory() as local_store:
            outfile = os.path.join(local_store, "tiny-data.txt")
            download(url, outfile, None)
            # Read stderr and make sure the progress bar is printed only when
            # told
            captured = capsys.readouterr()
            printed = captured.err.split("\r")[-1].strip()
            assert len(printed) == 79
            if sys.platform == "win32":
                progress = "100%|####################"
            else:
                progress = "100%|████████████████████"
            # Bar size is not always the same so can't reliably test the whole
            # bar.
            assert printed[:25] == progress
            # Check that the file was actually downloaded
            check_tiny_data(outfile)


@pytest.mark.network
@pytest.mark.skipif(tqdm is None, reason="requires tqdm")
@pytest.mark.skipif(paramiko is None, reason="requires paramiko")
def test_downloader_progressbar_sftp(capsys):
    "Setup an SFTP downloader function that prints a progress bar for fetch"
    downloader = SFTPDownloader(progressbar=True, username="demo", password="password")
    with TemporaryDirectory() as local_store:
        url = "sftp://test.rebex.net/pub/example/pocketftp.png"
        outfile = os.path.join(local_store, "pocketftp.png")
        downloader(url, outfile, None)
        # Read stderr and make sure the progress bar is printed only when told
        captured = capsys.readouterr()
        printed = captured.err.split("\r")[-1].strip()
        assert len(printed) == 79
        if sys.platform == "win32":
            progress = "100%|####################"
        else:
            progress = "100%|████████████████████"
        # Bar size is not always the same so can't reliably test the whole bar.
        assert printed[:25] == progress
        # Check that the file was actually downloaded
        assert os.path.exists(outfile)


@pytest.mark.network
def test_downloader_arbitrary_progressbar(capsys):
    "Setup a downloader function with an arbitrary progress bar class."

    class MinimalProgressDisplay:
        """A minimalist replacement for tqdm.tqdm"""

        def __init__(self, total):
            self.count = 0
            self.total = total

        def __repr__(self):
            """represent current completion"""
            return str(self.count) + "/" + str(self.total)

        def render(self):
            """print self.__repr__ to stderr"""
            print(f"\r{self}", file=sys.stderr, end="")

        def update(self, i):
            """modify completion and render"""
            self.count = i
            self.render()

        def reset(self):
            """set counter to 0"""
            self.count = 0

        @staticmethod
        def close():
            """print a new empty line"""
            print("", file=sys.stderr)

    pbar = MinimalProgressDisplay(total=None)
    download = HTTPDownloader(progressbar=pbar)
    with TemporaryDirectory() as local_store:
        fname = "large-data.txt"
        url = BASEURL + fname
        outfile = os.path.join(local_store, "large-data.txt")
        download(url, outfile, None)
        # Read stderr and make sure the progress bar is printed only when told
        captured = capsys.readouterr()
        printed = captured.err.split("\r")[-1].strip()

        progress = "336/336"
        assert printed == progress

        # Check that the downloaded file has the right content
        check_large_data(outfile)


class TestZenodoAPISupport:
    """
    Test support for different Zenodo APIs
    """

    article_id = 123456
    doi = f"10.0001/zenodo.{article_id}"
    doi_url = f"https://doi.org/{doi}"
    file_name = "my-file.zip"
    file_url = (
        "https://zenodo.org/api/files/513d7033-93a2-4eeb-821c-2fb0bbab0012/my-file.zip"
    )
    file_checksum = "2942bfabb3d05332b66eb128e0842cff"

    legacy_api_response = {
        "created": "2021-20-19T08:00:00.000000+00:00",
        "modified": "2021-20-19T08:00:00.000000+00:00",
        "id": article_id,
        "doi": doi,
        "doi_url": doi_url,
        "files": [
            {
                "id": "513d7033-93a2-4eeb-821c-2fb0bbab0012",
                "key": file_name,
                "checksum": f"md5:{file_checksum}",
                "links": {
                    "self": file_url,
                },
            }
        ],
    }

    new_api_response = {
        "created": "2021-20-19T08:00:00.000000+00:00",
        "modified": "2021-20-19T08:00:00.000000+00:00",
        "id": article_id,
        "doi": doi,
        "doi_url": doi_url,
        "files": [
            {
                "id": "513d7033-93a2-4eeb-821c-2fb0bbab0012",
                "filename": file_name,
                "checksum": file_checksum,
                "links": {
                    "self": file_url,
                },
            }
        ],
    }

    invalid_api_response = {
        "created": "2021-20-19T08:00:00.000000+00:00",
        "modified": "2021-20-19T08:00:00.000000+00:00",
        "id": article_id,
        "doi": doi,
        "doi_url": doi_url,
        "files": [
            {
                "id": "513d7033-93a2-4eeb-821c-2fb0bbab0012",
                "filename": file_name,
                "checksum": file_checksum,
                "links": {
                    "self": file_url,
                },
            },
            {
                "id": "513d7033-93a2-4eeb-821c-2fb0bbab0012",
                "key": file_name,
                "checksum": f"md5:{file_checksum}",
                "links": {
                    "self": file_url,
                },
            },
        ],
    }

    @pytest.mark.parametrize(
        "api_version, api_response",
        [
            ("legacy", legacy_api_response),
            ("new", new_api_response),
            ("invalid", invalid_api_response),
        ],
    )
    def test_api_version(self, httpserver, api_version, api_response):
        """
        Test if the API version is correctly detected.
        """
        # Create a local http server
        httpserver.expect_request(f"/zenodo.{self.article_id}").respond_with_json(
            api_response
        )
        # Create Zenodo downloader
        downloader = ZenodoRepository(doi=self.doi, archive_url=self.doi_url)
        # Override base url for the API of the downloader
        downloader.base_api_url = httpserver.url_for("")
        # Check if the API version is correctly identified
        if api_version != "invalid":
            assert downloader.api_version == api_version
        else:
            msg = "Couldn't determine the version of the Zenodo API"
            with pytest.raises(ValueError, match=msg):
                api_version = downloader.api_version

    @pytest.mark.parametrize(
        "api_version, api_response",
        [("legacy", legacy_api_response), ("new", new_api_response)],
    )
    def test_download_url(self, httpserver, api_version, api_response):
        """
        Test if the download url is correct for each API version.
        """
        # Create a local http server
        httpserver.expect_request(f"/zenodo.{self.article_id}").respond_with_json(
            api_response
        )
        # Create Zenodo downloader
        downloader = ZenodoRepository(doi=self.doi, archive_url=self.doi_url)
        # Override base url for the API of the downloader
        downloader.base_api_url = httpserver.url_for("")
        # Check if the download url is correct
        download_url = downloader.download_url(file_name=self.file_name)
        if api_version == "legacy":
            assert download_url == self.file_url
        else:
            expected_url = (
                "https://zenodo.org/records/"
                f"{self.article_id}/files/{self.file_name}?download=1"
            )
            assert download_url == expected_url

    @pytest.mark.parametrize(
        "api_response",
        [legacy_api_response, new_api_response],
    )
    def test_populate_registry(self, httpserver, tmp_path, api_response):
        """
        Test if population of registry is correctly done for each API version.
        """
        # Create a local http server
        httpserver.expect_request(f"/zenodo.{self.article_id}").respond_with_json(
            api_response
        )
        # Create sample pooch object
        puppy = Pooch(base_url="", path=tmp_path)
        # Create Zenodo downloader
        downloader = ZenodoRepository(doi=self.doi, archive_url=self.doi_url)
        # Override base url for the API of the downloader
        downloader.base_api_url = httpserver.url_for("")
        # Populate registry
        downloader.populate_registry(puppy)
        assert puppy.registry == {self.file_name: f"md5:{self.file_checksum}"}
