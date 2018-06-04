from __future__ import print_function, division, absolute_import

import posixpath

from .glob import generic_glob
from ..base import tokenize

import pyarrow as pa


class HDFS3Wrapper(pa.filesystem.DaskFileSystem):
    """Wrapper around `hdfs3.HDFileSystem` that allows it to be passed to
    pyarrow methods"""
    def isdir(self, path):
        return self.fs.isdir(path)

    def isfile(self, path):
        return self.fs.isfile(path)


_MIN_PYARROW_VERSION_SUPPORTED = '0.8.1.dev81'


class PyArrowHadoopFileSystem(object):
    sep = "/"

    def __init__(self, **kwargs):
        self.fs = pa.hdfs.HadoopFileSystem(**kwargs)

    @classmethod
    def from_pyarrow(cls, fs):
        out = object.__new__(cls)
        out.fs = fs
        return out

    def open(self, path, mode='rb', **kwargs):
        return self.fs.open(path, mode=mode, **kwargs)

    def glob(self, path):
        return sorted(generic_glob(self.fs, posixpath, path))

    def mkdirs(self, path):
        return self.fs.mkdir(path, create_parents=True)

    def ukey(self, path):
        return tokenize(path, self.fs.info(path)['last_modified'])

    def size(self, path):
        return self.fs.info(path)['size']

    def _get_pyarrow_filesystem(self):
        return self.fs
