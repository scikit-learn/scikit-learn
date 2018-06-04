from __future__ import print_function, division, absolute_import

from glob import glob
import os

from . import core
from ..base import tokenize


class LocalFileSystem(core.FileSystem):
    """API spec for the methods a filesystem

    A filesystem must provide these methods, if it is to be registered as
    a backend for dask.

    Implementation for local disc"""
    sep = os.sep

    def __init__(self, **storage_options):
        """
        Parameters
        ----------
        storage_options: key-value
            May be credentials, or other configuration specific to the backend.
        """
        self.cwd = os.getcwd()

    def _normalize_path(self, path):
        """Ensure paths are absolute and normalized"""
        if not os.path.isabs(path):
            return os.path.join(self.cwd, path)
        return os.path.normpath(path)

    def glob(self, path):
        """For a template path, return matching files"""
        return sorted(glob(self._normalize_path(path)))

    def mkdirs(self, path):
        """Make any intermediate directories to make path writable"""
        path = self._normalize_path(path)
        try:
            os.makedirs(path)
        except OSError:
            assert os.path.isdir(path)

    def open(self, path, mode='rb', **kwargs):
        """Make a file-like object

        Parameters
        ----------
        mode: string
            normally "rb", "wb" or "ab" or other.
        kwargs: key-value
            Any other parameters, such as buffer size. May be better to set
            these on the filesystem instance, to apply to all files created by
            it. Not used for local.
        """
        return open(self._normalize_path(path), mode=mode)

    def ukey(self, path):
        """Unique identifier, so we can tell if a file changed"""
        path = self._normalize_path(path)
        return tokenize(path, os.stat(path).st_mtime)

    def size(self, path):
        """Size in bytes of the file at path"""
        return os.stat(self._normalize_path(path)).st_size

    def _get_pyarrow_filesystem(self):
        """Get an equivalent pyarrow filesystem"""
        import pyarrow as pa
        return pa.filesystem.LocalFileSystem.get_instance()


core._filesystems['file'] = LocalFileSystem
