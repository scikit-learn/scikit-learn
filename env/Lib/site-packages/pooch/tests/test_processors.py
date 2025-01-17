# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
"""
Test the processor hooks
"""
from pathlib import Path
from tempfile import TemporaryDirectory
import warnings

import pytest

from .. import Pooch
from ..processors import Unzip, Untar, Decompress

from .utils import pooch_test_url, pooch_test_registry, check_tiny_data, capture_log


REGISTRY = pooch_test_registry()
BASEURL = pooch_test_url()


@pytest.mark.network
@pytest.mark.parametrize(
    "method,ext,name",
    [
        ("auto", "xz", None),
        ("lzma", "xz", None),
        ("xz", "xz", None),
        ("bzip2", "bz2", None),
        ("gzip", "gz", None),
        ("gzip", "gz", "different-name.txt"),
    ],
    ids=["auto", "lzma", "xz", "bz2", "gz", "name"],
)
def test_decompress(method, ext, name):
    "Check that decompression after download works for all formats"
    processor = Decompress(method=method, name=name)
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        if name is None:
            true_path = str(path / ".".join(["tiny-data.txt", ext, "decomp"]))
        else:
            true_path = str(path / name)
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY)
        # Check the logs when downloading and from the processor
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt." + ext, processor=processor)
            logs = log_file.getvalue()
            lines = logs.splitlines()
            assert len(lines) == 2
            assert lines[0].split()[0] == "Downloading"
            assert lines[-1].startswith("Decompressing")
            assert method in lines[-1]
        assert fname == true_path
        check_tiny_data(fname)
        # Check that processor doesn't execute when not downloading
        with capture_log() as log_file:
            fname = pup.fetch("tiny-data.txt." + ext, processor=processor)
            assert log_file.getvalue() == ""
        assert fname == true_path
        check_tiny_data(fname)


@pytest.mark.network
def test_decompress_fails():
    "Should fail if method='auto' and no extension is given in the file name"
    with TemporaryDirectory() as local_store:
        path = Path(local_store)
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY)
        # Invalid extension
        with pytest.raises(ValueError) as exception:
            with warnings.catch_warnings():
                pup.fetch("tiny-data.txt", processor=Decompress(method="auto"))
        assert exception.value.args[0].startswith("Unrecognized file extension '.txt'")
        assert "pooch.Unzip/Untar" not in exception.value.args[0]
        # Should also fail for a bad method name
        with pytest.raises(ValueError) as exception:
            with warnings.catch_warnings():
                pup.fetch("tiny-data.txt", processor=Decompress(method="bla"))
        assert exception.value.args[0].startswith("Invalid compression method 'bla'")
        assert "pooch.Unzip/Untar" not in exception.value.args[0]
        # Point people to Untar and Unzip
        with pytest.raises(ValueError) as exception:
            with warnings.catch_warnings():
                pup.fetch("tiny-data.txt", processor=Decompress(method="zip"))
        assert exception.value.args[0].startswith("Invalid compression method 'zip'")
        assert "pooch.Unzip/Untar" in exception.value.args[0]
        with pytest.raises(ValueError) as exception:
            with warnings.catch_warnings():
                pup.fetch("store.zip", processor=Decompress(method="auto"))
        assert exception.value.args[0].startswith("Unrecognized file extension '.zip'")
        assert "pooch.Unzip/Untar" in exception.value.args[0]


@pytest.mark.network
@pytest.mark.parametrize(
    "target_path", [None, "some_custom_path"], ids=["default_path", "custom_path"]
)
@pytest.mark.parametrize(
    "archive,members",
    [
        ("tiny-data", ["tiny-data.txt"]),
        ("store", None),
        ("store", ["store/tiny-data.txt"]),
        ("store", ["store/subdir/tiny-data.txt"]),
        ("store", ["store/subdir"]),
        ("store", ["store/tiny-data.txt", "store/subdir"]),
    ],
    ids=[
        "single_file",
        "archive_all",
        "archive_file",
        "archive_subdir_file",
        "archive_subdir",
        "archive_multiple",
    ],
)
@pytest.mark.parametrize(
    "processor_class,extension",
    [(Unzip, ".zip"), (Untar, ".tar.gz")],
    ids=["Unzip", "Untar"],
)
def test_unpacking(processor_class, extension, target_path, archive, members):
    "Tests the behaviour of processors for unpacking archives (Untar, Unzip)"
    processor = processor_class(members=members, extract_dir=target_path)
    if target_path is None:
        target_path = archive + extension + processor.suffix
    with TemporaryDirectory() as path:
        path = Path(path)
        true_paths, expected_log = _unpacking_expected_paths_and_logs(
            archive, members, path / target_path, processor_class.__name__
        )
        # Setup a pooch in a temp dir
        pup = Pooch(path=path, base_url=BASEURL, registry=REGISTRY)
        # Capture logs and check for the right processor message
        with capture_log() as log_file:
            fnames = pup.fetch(archive + extension, processor=processor)
            assert set(fnames) == true_paths
            _check_logs(log_file, expected_log)
        for fname in fnames:
            check_tiny_data(fname)
        # Check that processor doesn't execute when not downloading
        with capture_log() as log_file:
            fnames = pup.fetch(archive + extension, processor=processor)
            assert set(fnames) == true_paths
            _check_logs(log_file, [])
        for fname in fnames:
            check_tiny_data(fname)


@pytest.mark.network
@pytest.mark.parametrize(
    "processor_class,extension",
    [(Unzip, ".zip"), (Untar, ".tar.gz")],
)
def test_multiple_unpacking(processor_class, extension):
    "Test that multiple subsequent calls to a processor yield correct results"

    with TemporaryDirectory() as local_store:
        pup = Pooch(path=Path(local_store), base_url=BASEURL, registry=REGISTRY)

        # Do a first fetch with the one member only
        processor1 = processor_class(members=["store/tiny-data.txt"])
        filenames1 = pup.fetch("store" + extension, processor=processor1)
        assert len(filenames1) == 1
        check_tiny_data(filenames1[0])

        # Do a second fetch with the other member
        processor2 = processor_class(
            members=["store/tiny-data.txt", "store/subdir/tiny-data.txt"]
        )
        filenames2 = pup.fetch("store" + extension, processor=processor2)
        assert len(filenames2) == 2
        check_tiny_data(filenames2[0])
        check_tiny_data(filenames2[1])

        # Do a third fetch, again with one member and assert
        # that only this member was returned
        filenames3 = pup.fetch("store" + extension, processor=processor1)
        assert len(filenames3) == 1
        check_tiny_data(filenames3[0])


@pytest.mark.network
@pytest.mark.parametrize(
    "processor_class,extension",
    [(Unzip, ".zip"), (Untar, ".tar.gz")],
)
def test_unpack_members_with_leading_dot(processor_class, extension):
    "Test that unpack members can also be specifed both with a leading ./"

    with TemporaryDirectory() as local_store:
        pup = Pooch(path=Path(local_store), base_url=BASEURL, registry=REGISTRY)

        # Do a first fetch with the one member only
        processor1 = processor_class(members=["./store/tiny-data.txt"])
        filenames1 = pup.fetch("store" + extension, processor=processor1)
        assert len(filenames1) == 1
        check_tiny_data(filenames1[0])


def _check_logs(log_file, expected_lines):
    """
    Assert that the lines in the log match the expected ones.
    """
    lines = log_file.getvalue().splitlines()
    assert len(lines) == len(expected_lines)
    for line, expected_line in zip(lines, expected_lines):
        assert line.startswith(expected_line)


def _unpacking_expected_paths_and_logs(archive, members, path, name):
    """
    Generate the appropriate expected paths and log message depending on the
    parameters for the test.
    """
    log_lines = ["Downloading"]
    if archive == "tiny-data":
        true_paths = {str(path / "tiny-data.txt")}
        log_lines.append("Extracting 'tiny-data.txt'")
    elif archive == "store" and members is None:
        true_paths = {
            str(path / "store" / "tiny-data.txt"),
            str(path / "store" / "subdir" / "tiny-data.txt"),
        }
        log_lines.append(f"{name}{name[-1]}ing contents")
    elif archive == "store" and members is not None:
        true_paths = []
        for member in members:
            true_path = path / Path(*member.split("/"))
            if not str(true_path).endswith("tiny-data.txt"):
                true_path = true_path / "tiny-data.txt"
            true_paths.append(str(true_path))
            log_lines.append(f"Extracting '{member}'")
        true_paths = set(true_paths)
    return true_paths, log_lines


@pytest.mark.network
@pytest.mark.parametrize(
    "processor_class,extension",
    [(Unzip, ".zip"), (Untar, ".tar.gz")],
)
def test_unpacking_members_then_no_members(processor_class, extension):
    """
    Test that calling with valid members then without them works.
    https://github.com/fatiando/pooch/issues/364
    """
    with TemporaryDirectory() as local_store:
        pup = Pooch(path=Path(local_store), base_url=BASEURL, registry=REGISTRY)

        # Do a first fetch with an existing member
        processor1 = processor_class(members=["store/tiny-data.txt"])
        filenames1 = pup.fetch("store" + extension, processor=processor1)
        assert len(filenames1) == 1

        # Do a second fetch with no members
        processor2 = processor_class()
        filenames2 = pup.fetch("store" + extension, processor=processor2)
        assert len(filenames2) > 1


@pytest.mark.network
@pytest.mark.parametrize(
    "processor_class,extension",
    [(Unzip, ".zip"), (Untar, ".tar.gz")],
)
def test_unpacking_wrong_members_then_no_members(processor_class, extension):
    """
    Test that calling with invalid members then without them works.
    https://github.com/fatiando/pooch/issues/364
    """
    with TemporaryDirectory() as local_store:
        pup = Pooch(path=Path(local_store), base_url=BASEURL, registry=REGISTRY)

        # Do a first fetch with incorrect member
        processor1 = processor_class(members=["not-a-valid-file.csv"])
        filenames1 = pup.fetch("store" + extension, processor=processor1)
        assert len(filenames1) == 0

        # Do a second fetch with no members
        processor2 = processor_class()
        filenames2 = pup.fetch("store" + extension, processor=processor2)
        assert len(filenames2) > 0
