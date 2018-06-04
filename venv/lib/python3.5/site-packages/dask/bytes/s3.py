from __future__ import print_function, division, absolute_import

from s3fs import S3FileSystem

from . import core


class DaskS3FileSystem(S3FileSystem, core.FileSystem):

    sep = '/'

    def __init__(self, key=None, username=None, secret=None, password=None,
                 path=None, host=None, s3=None, **kwargs):
        if username is not None:
            if key is not None:
                raise KeyError("S3 storage options got secrets argument "
                               "collision. Please, use either `key` "
                               "storage option or password field in URLpath, "
                               "not both options together.")
            key = username
        if key is not None:
            kwargs['key'] = key
        if password is not None:
            if secret is not None:
                raise KeyError("S3 storage options got secrets argument "
                               "collision. Please, use either `secret` "
                               "storage option or password field in URLpath, "
                               "not both options together.")
            secret = password
        if secret is not None:
            kwargs['secret'] = secret
        super(DaskS3FileSystem, self).__init__(**kwargs)

    def mkdirs(self, path):
        pass  # no need to pre-make paths on S3

    def ukey(self, path):
        return self.info(path)['ETag']

    def size(self, path):
        return self.info(path)['Size']

    def _get_pyarrow_filesystem(self):
        """Get an equivalent pyarrow fileystem"""
        import pyarrow as pa
        return pa.filesystem.S3FSWrapper(self)


core._filesystems['s3'] = DaskS3FileSystem
