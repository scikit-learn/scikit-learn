from __future__ import print_function, division, absolute_import

import posixpath

from .glob import generic_glob
from ..base import tokenize

import hdfs3


class HDFS3HadoopFileSystem(object):
    sep = "/"

    def __init__(self, **kwargs):
        self.fs = hdfs3.HDFileSystem(**kwargs)

    @classmethod
    def from_hdfs3(cls, fs):
        out = object.__new__(cls)
        out.fs = fs
        return out

    def open(self, path, mode='rb', **kwargs):
        return self.fs.open(path, mode=mode, **kwargs)

    def glob(self, path):
        return sorted(generic_glob(self.fs, posixpath, path))

    def mkdirs(self, path):
        return self.fs.makedirs(path)

    def ukey(self, path):
        return tokenize(path, self.fs.info(path)['last_mod'])

    def size(self, path):
        return self.fs.info(path)['size']

    def _get_pyarrow_filesystem(self):
        from .pyarrow import HDFS3Wrapper
        return HDFS3Wrapper(self.fs)
