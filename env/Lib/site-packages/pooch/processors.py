# Copyright (c) 2018 The Pooch Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
#
# This code is part of the Fatiando a Terra project (https://www.fatiando.org)
#
# pylint: disable=line-too-long
"""
Post-processing hooks
"""
import abc
import os
import bz2
import gzip
import lzma
import shutil
from zipfile import ZipFile
from tarfile import TarFile

from .utils import get_logger


class ExtractorProcessor(abc.ABC):  # pylint: disable=too-few-public-methods
    """
    Abstract base class for extractions from compressed archives.

    Subclasses can be used with :meth:`pooch.Pooch.fetch` and
    :func:`pooch.retrieve` to unzip a downloaded data file into a folder in the
    local data store. :meth:`~pooch.Pooch.fetch` will return a list with the
    names of the extracted files instead of the archive.

    Parameters
    ----------
    members : list or None
        If None, will unpack all files in the archive. Otherwise, *members*
        must be a list of file names to unpack from the archive. Only these
        files will be unpacked.
    extract_dir : str or None
        If None, files will be unpacked to the default location (a folder in
        the same location as the downloaded zip file, with a suffix added).
        Otherwise, files will be unpacked to ``extract_dir``, which is
        interpreted as a *relative path* (relative to the cache location
        provided by :func:`pooch.retrieve` or :meth:`pooch.Pooch.fetch`).

    """

    def __init__(self, members=None, extract_dir=None):
        self.members = members
        self.extract_dir = extract_dir

    @property
    @abc.abstractmethod
    def suffix(self):
        """
        String appended to unpacked archive folder name.
        Only used if extract_dir is None.
        MUST BE IMPLEMENTED BY CHILD CLASSES.
        """

    @abc.abstractmethod
    def _all_members(self, fname):
        """
        Return all the members in the archive.
        MUST BE IMPLEMENTED BY CHILD CLASSES.
        """

    @abc.abstractmethod
    def _extract_file(self, fname, extract_dir):
        """
        This method receives an argument for the archive to extract and the
        destination path.
        MUST BE IMPLEMENTED BY CHILD CLASSES.
        """

    def __call__(self, fname, action, pooch):
        """
        Extract all files from the given archive.

        Parameters
        ----------
        fname : str
            Full path of the zipped file in local storage.
        action : str
            Indicates what action was taken by :meth:`pooch.Pooch.fetch` or
            :func:`pooch.retrieve`:

            * ``"download"``: File didn't exist locally and was downloaded
            * ``"update"``: Local file was outdated and was re-download
            * ``"fetch"``: File exists and is updated so it wasn't downloaded

        pooch : :class:`pooch.Pooch`
            The instance of :class:`pooch.Pooch` that is calling this.

        Returns
        -------
        fnames : list of str
            A list of the full path to all files in the extracted archive.

        """
        if self.extract_dir is None:
            self.extract_dir = fname + self.suffix
        else:
            archive_dir = fname.rsplit(os.path.sep, maxsplit=1)[0]
            self.extract_dir = os.path.join(archive_dir, self.extract_dir)
        # Get a list of everyone who is supposed to be in the unpacked folder
        # so we can check if they are all there or if we need to extract new
        # files.
        if self.members is None or not self.members:
            members = self._all_members(fname)
        else:
            members = self.members
        if (
            (action in ("update", "download"))
            or (not os.path.exists(self.extract_dir))
            or not all(
                os.path.exists(os.path.join(self.extract_dir, m)) for m in members
            )
        ):
            # Make sure that the folder with the extracted files exists
            os.makedirs(self.extract_dir, exist_ok=True)
            self._extract_file(fname, self.extract_dir)

        # Get a list of all file names (including subdirectories) in our folder
        # of unzipped files, filtered by the given members list
        fnames = []
        for path, _, files in os.walk(self.extract_dir):
            for filename in files:
                relpath = os.path.normpath(
                    os.path.join(os.path.relpath(path, self.extract_dir), filename)
                )
                if self.members is None or any(
                    relpath.startswith(os.path.normpath(m)) for m in self.members
                ):
                    fnames.append(os.path.join(path, filename))

        return fnames


class Unzip(ExtractorProcessor):  # pylint: disable=too-few-public-methods
    """
    Processor that unpacks a zip archive and returns a list of all files.

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to unzip a
    downloaded data file into a folder in the local data store. The
    method/function will return a list with the names of the unzipped files
    instead of the zip archive.

    The output folder is ``{fname}.unzip``.

    Parameters
    ----------
    members : list or None
        If None, will unpack all files in the zip archive. Otherwise, *members*
        must be a list of file names to unpack from the archive. Only these
        files will be unpacked.
    extract_dir : str or None
        If None, files will be unpacked to the default location (a folder in
        the same location as the downloaded zip file, with the suffix
        ``.unzip`` added). Otherwise, files will be unpacked to
        ``extract_dir``, which is interpreted as a *relative path* (relative to
        the cache location provided by :func:`pooch.retrieve` or
        :meth:`pooch.Pooch.fetch`).

    """

    @property
    def suffix(self):
        """
        String appended to unpacked archive folder name.
        Only used if extract_dir is None.
        """
        return ".unzip"

    def _all_members(self, fname):
        """Return all members from a given archive."""
        with ZipFile(fname, "r") as zip_file:
            return zip_file.namelist()

    def _extract_file(self, fname, extract_dir):
        """
        This method receives an argument for the archive to extract and the
        destination path.
        """
        with ZipFile(fname, "r") as zip_file:
            if self.members is None:
                get_logger().info(
                    "Unzipping contents of '%s' to '%s'", fname, extract_dir
                )
                # Unpack all files from the archive into our new folder
                zip_file.extractall(path=extract_dir)
            else:
                for member in self.members:
                    get_logger().info(
                        "Extracting '%s' from '%s' to '%s'", member, fname, extract_dir
                    )
                    # If the member is a dir, we need to get the names of the
                    # elements it contains for extraction (ZipFile does not
                    # support dirs on .extract). If it's not a dir, this will
                    # only include the member itself.
                    # Based on:
                    # https://stackoverflow.com/questions/8008829/extract-only-a-single-directory-from-tar
                    subdir_members = [
                        name
                        for name in zip_file.namelist()
                        if os.path.normpath(name).startswith(os.path.normpath(member))
                    ]
                    # Extract the data file from within the archive
                    zip_file.extractall(members=subdir_members, path=extract_dir)


class Untar(ExtractorProcessor):  # pylint: disable=too-few-public-methods
    """
    Processor that unpacks a tar archive and returns a list of all files.

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to untar a
    downloaded data file into a folder in the local data store. The
    method/function will return a list with the names of the extracted files
    instead of the archive.

    The output folder is ``{fname}.untar``.


    Parameters
    ----------
    members : list or None
        If None, will unpack all files in the archive. Otherwise, *members*
        must be a list of file names to unpack from the archive. Only these
        files will be unpacked.
    extract_dir : str or None
        If None, files will be unpacked to the default location (a folder in
        the same location as the downloaded tar file, with the suffix
        ``.untar`` added). Otherwise, files will be unpacked to
        ``extract_dir``, which is interpreted as a *relative path* (relative to
        the cache location  provided by :func:`pooch.retrieve` or
        :meth:`pooch.Pooch.fetch`).
    """

    @property
    def suffix(self):
        """
        String appended to unpacked archive folder name.
        Only used if extract_dir is None.
        """
        return ".untar"

    def _all_members(self, fname):
        """Return all members from a given archive."""
        with TarFile.open(fname, "r") as tar_file:
            return [info.name for info in tar_file.getmembers()]

    def _extract_file(self, fname, extract_dir):
        """
        This method receives an argument for the archive to extract and the
        destination path.
        """
        with TarFile.open(fname, "r") as tar_file:
            if self.members is None:
                get_logger().info(
                    "Untarring contents of '%s' to '%s'", fname, extract_dir
                )
                # Unpack all files from the archive into our new folder
                tar_file.extractall(path=extract_dir)
            else:
                for member in self.members:
                    get_logger().info(
                        "Extracting '%s' from '%s' to '%s'", member, fname, extract_dir
                    )
                    # If the member is a dir, we need to get the names of the
                    # elements it contains for extraction (TarFile does not
                    # support dirs on .extract). If it's not a dir, this will
                    # only include the member itself.
                    # Based on:
                    # https://stackoverflow.com/questions/8008829/extract-only-a-single-directory-from-tar
                    # Can't use .getnames because extractall expects TarInfo
                    # objects.
                    subdir_members = [
                        info
                        for info in tar_file.getmembers()
                        if os.path.normpath(info.name).startswith(
                            os.path.normpath(member)
                        )
                    ]
                    # Extract the data file from within the archive
                    tar_file.extractall(members=subdir_members, path=extract_dir)


class Decompress:  # pylint: disable=too-few-public-methods
    """
    Processor that decompress a file and returns the decompressed version.

    Use with :meth:`pooch.Pooch.fetch` or :func:`pooch.retrieve` to decompress
    a downloaded data file so that it can be easily opened. Useful for data
    files that take a long time to decompress (exchanging disk space for
    speed).

    Supported decompression methods are LZMA (``.xz``), bzip2 (``.bz2``), and
    gzip (``.gz``).

    File names with the standard extensions (see above) can use
    ``method="auto"`` to automatically determine the compression method. This
    can be overwritten by setting the *method* argument.

    .. note::

        To unpack zip and tar archives with one or more files, use
        :class:`pooch.Unzip` and :class:`pooch.Untar` instead.

    The output file is ``{fname}.decomp`` by default but it can be changed by
    setting the ``name`` parameter.

    .. warning::

        Passing in ``name`` can cause existing data to be lost! For example, if
        a file already exists with the specified name it will be overwritten
        with the new decompressed file content. **Use this option with
        caution.**

    Parameters
    ----------
    method : str
        Name of the compression method. Can be "auto", "lzma", "xz", "bzip2",
        or "gzip".
    name : None or str
        Defines the decompressed file name. The file name will be
        ``{fname}.decomp`` if ``None`` (default) or the given name otherwise.
        Note that the name should **not** include the full (or relative) path,
        it should be just the file name itself.

    """

    modules = {"auto": None, "lzma": lzma, "xz": lzma, "gzip": gzip, "bzip2": bz2}
    extensions = {".xz": "lzma", ".gz": "gzip", ".bz2": "bzip2"}

    def __init__(self, method="auto", name=None):
        self.method = method
        self.name = name

    def __call__(self, fname, action, pooch):
        """
        Decompress the given file.

        The output file will be either ``{fname}.decomp`` or the given *name*
        class attribute.

        Parameters
        ----------
        fname : str
            Full path of the compressed file in local storage.
        action : str
            Indicates what action was taken by :meth:`pooch.Pooch.fetch` or
            :func:`pooch.retrieve`:

            - ``"download"``: File didn't exist locally and was downloaded
            - ``"update"``: Local file was outdated and was re-download
            - ``"fetch"``: File exists and is updated so it wasn't downloaded

        pooch : :class:`pooch.Pooch`
            The instance of :class:`pooch.Pooch` that is calling this.

        Returns
        -------
        fname : str
            The full path to the decompressed file.
        """
        if self.name is None:
            decompressed = fname + ".decomp"
        else:
            decompressed = os.path.join(os.path.dirname(fname), self.name)
        if action in ("update", "download") or not os.path.exists(decompressed):
            get_logger().info(
                "Decompressing '%s' to '%s' using method '%s'.",
                fname,
                decompressed,
                self.method,
            )
            module = self._compression_module(fname)
            with open(decompressed, "w+b") as output:
                with module.open(fname) as compressed:
                    shutil.copyfileobj(compressed, output)
        return decompressed

    def _compression_module(self, fname):
        """
        Get the Python module compatible with fname and the chosen method.

        If the *method* attribute is "auto", will select a method based on the
        extension. If no recognized extension is in the file name, will raise a
        ValueError.
        """
        error_archives = "To unpack zip/tar archives, use pooch.Unzip/Untar instead."
        if self.method not in self.modules:
            message = (
                f"Invalid compression method '{self.method}'. "
                f"Must be one of '{list(self.modules.keys())}'."
            )
            if self.method in {"zip", "tar"}:
                message = " ".join([message, error_archives])
            raise ValueError(message)
        if self.method == "auto":
            ext = os.path.splitext(fname)[-1]
            if ext not in self.extensions:
                message = (
                    f"Unrecognized file extension '{ext}'. "
                    f"Must be one of '{list(self.extensions.keys())}'."
                )
                if ext in {".zip", ".tar"}:
                    message = " ".join([message, error_archives])
                raise ValueError(message)
            return self.modules[self.extensions[ext]]
        return self.modules[self.method]
