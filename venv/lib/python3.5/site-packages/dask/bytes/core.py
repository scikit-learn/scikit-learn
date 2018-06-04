from __future__ import print_function, division, absolute_import

import io
import os
from distutils.version import LooseVersion
from warnings import warn

from toolz import merge

from .compression import seekable_files, files as compress_files
from .utils import (SeekableFile, read_block, infer_compression,
                    infer_storage_options, build_name_function,
                    update_storage_options)
from ..compatibility import unicode
from ..context import _globals
from ..base import tokenize
from ..delayed import delayed
from ..utils import import_required, is_integer


def read_bytes(urlpath, delimiter=None, not_zero=False, blocksize=2**27,
               sample=True, compression=None, **kwargs):
    """Given a path or paths, return delayed objects that read from those paths.

    The path may be a filename like ``'2015-01-01.csv'`` or a globstring
    like ``'2015-*-*.csv'``.

    The path may be preceded by a protocol, like ``s3://`` or ``hdfs://`` if
    those libraries are installed.

    This cleanly breaks data by a delimiter if given, so that block boundaries
    start directly after a delimiter and end on the delimiter.

    Parameters
    ----------
    urlpath : string or list
        Absolute or relative filepath(s). Prefix with a protocol like ``s3://``
        to read from alternative filesystems. To read from multiple files you
        can pass a globstring or a list of paths, with the caveat that they
        must all have the same protocol.
    delimiter : bytes
        An optional delimiter, like ``b'\\n'`` on which to split blocks of
        bytes.
    not_zero : bool
        Force seek of start-of-file delimiter, discarding header.
    blocksize : int (=128MB)
        Chunk size in bytes
    compression : string or None
        String like 'gzip' or 'xz'.  Must support efficient random access.
    sample : bool or int
        Whether or not to return a header sample. If an integer is given it is
        used as sample size, otherwise the default sample size is 10kB.
    **kwargs : dict
        Extra options that make sense to a particular storage connection, e.g.
        host, port, username, password, etc.

    Examples
    --------
    >>> sample, blocks = read_bytes('2015-*-*.csv', delimiter=b'\\n')  # doctest: +SKIP
    >>> sample, blocks = read_bytes('s3://bucket/2015-*-*.csv', delimiter=b'\\n')  # doctest: +SKIP

    Returns
    -------
    sample : bytes
        The sample header
    blocks : list of lists of ``dask.Delayed``
        Each list corresponds to a file, and each delayed object computes to a
        block of bytes from that file.
    """
    fs, fs_token, paths = get_fs_token_paths(urlpath, mode='rb',
                                             storage_options=kwargs)

    if len(paths) == 0:
        raise IOError("%s resolved to no files" % urlpath)

    if blocksize is not None:
        if not is_integer(blocksize):
            raise TypeError("blocksize must be an integer")
        blocksize = int(blocksize)

    if blocksize is None:
        offsets = [[0]] * len(paths)
        lengths = [[None]] * len(paths)
    else:
        offsets = []
        lengths = []
        for path in paths:
            try:
                size = logical_size(fs, path, compression)
            except ValueError:
                raise ValueError("Cannot do chunked reads on files compressed "
                                 "with compression=%r. To read, set "
                                 "blocksize=None" % get_compression(path, compression))
            off = list(range(0, size, blocksize))
            length = [blocksize] * len(off)
            if not_zero:
                off[0] = 1
                length[0] -= 1
            offsets.append(off)
            lengths.append(length)

    delayed_read = delayed(read_block_from_file)

    out = []
    for path, offset, length in zip(paths, offsets, lengths):
        token = tokenize(fs_token, delimiter, path, fs.ukey(path),
                         compression, offset)
        keys = ['read-block-%s-%s' % (o, token) for o in offset]
        out.append([delayed_read(OpenFile(fs, path, compression=compression),
                                 o, l, delimiter, dask_key_name=key)
                    for o, key, l in zip(offset, keys, length)])

    if sample:
        with OpenFile(fs, paths[0], compression=compression) as f:
            nbytes = 10000 if sample is True else sample
            sample = read_block(f, 0, nbytes, delimiter)

    return sample, out


def read_block_from_file(lazy_file, off, bs, delimiter):
    with lazy_file as f:
        return read_block(f, off, bs, delimiter)


class OpenFile(object):
    """
    File-like object to be used in a context

    These instances are safe to serialize, as the low-level file object
    is not created until invoked using `with`.

    Parameters
    ----------
    fs : FileSystem
        The file system to use for opening the file. Should match the interface
        of ``dask.bytes.local.LocalFileSystem``.
    path : str
        Location to open
    mode : str like 'rb', optional
        Mode of the opened file
    compression : str or None, optional
        Compression to apply
    encoding : str or None, optional
        The encoding to use if opened in text mode.
    errors : str or None, optional
        How to handle encoding errors if opened in text mode.
    """
    def __init__(self, fs, path, mode='rb', compression=None, encoding=None,
                 errors=None):
        self.fs = fs
        self.path = path
        self.mode = mode
        self.compression = get_compression(path, compression)
        self.encoding = encoding
        self.errors = errors
        self.fobjects = []

    def __reduce__(self):
        return (OpenFile, (self.fs, self.path, self.mode, self.compression,
                           self.encoding, self.errors))

    def __enter__(self):
        mode = self.mode.replace('t', '').replace('b', '') + 'b'

        f = SeekableFile(self.fs.open(self.path, mode=mode))

        fobjects = [f]

        if self.compression is not None:
            compress = merge(seekable_files, compress_files)[self.compression]
            f = compress(f, mode=mode)
            fobjects.append(f)

        if 't' in self.mode:
            f = io.TextIOWrapper(f, encoding=self.encoding,
                                 errors=self.errors)
            fobjects.append(f)

        self.fobjects = fobjects
        return f

    def __exit__(self, *args):
        self.close()

    def close(self):
        """Close all encapsulated file objects"""
        for f in reversed(self.fobjects):
            f.close()
        self.fobjects = []


def open_files(urlpath, mode='rb', compression=None, encoding='utf8',
               errors=None, name_function=None, num=1, **kwargs):
    """ Given a path or paths, return a list of ``OpenFile`` objects.

    Parameters
    ----------
    urlpath : string or list
        Absolute or relative filepath(s). Prefix with a protocol like ``s3://``
        to read from alternative filesystems. To read from multiple files you
        can pass a globstring or a list of paths, with the caveat that they
        must all have the same protocol.
    mode : 'rb', 'wt', etc.
    compression : string
        Compression to use.  See ``dask.bytes.compression.files`` for options.
    encoding : str
        For text mode only
    errors : None or str
        Passed to TextIOWrapper in text mode
    name_function : function or None
        if opening a set of files for writing, those files do not yet exist,
        so we need to generate their names by formatting the urlpath for
        each sequence number
    num : int [1]
        if writing mode, number of files we expect to create (passed to
        name+function)
    **kwargs : dict
        Extra options that make sense to a particular storage connection, e.g.
        host, port, username, password, etc.

    Examples
    --------
    >>> files = open_files('2015-*-*.csv')  # doctest: +SKIP
    >>> files = open_files('s3://bucket/2015-*-*.csv.gz', compression='gzip')  # doctest: +SKIP

    Returns
    -------
    List of ``OpenFile`` objects.
    """
    fs, fs_token, paths = get_fs_token_paths(urlpath, mode, num=num,
                                             name_function=name_function,
                                             storage_options=kwargs)

    return [OpenFile(fs, path, mode=mode, compression=compression,
                     encoding=encoding, errors=errors)
            for path in paths]


def get_compression(urlpath, compression):
    if compression == 'infer':
        compression = infer_compression(urlpath)
    if compression is not None and compression not in compress_files:
        raise ValueError("Compression type %s not supported" % compression)
    return compression


def infer_options(urlpath):
    if hasattr(urlpath, 'name'):
        # deal with pathlib.Path objects - must be local
        urlpath = str(urlpath)
        ispath = True
    else:
        ispath = False

    options = infer_storage_options(urlpath)
    protocol = options.pop('protocol')
    urlpath = options.pop('path')

    if ispath and protocol != 'file':
        raise ValueError("Only use pathlib.Path with local files.")

    return urlpath, protocol, options


def get_fs_token_paths(urlpath, mode='rb', num=1, name_function=None,
                       storage_options=None):
    """Filesystem, deterministic token, and paths from a urlpath and options.

    Parameters
    ----------
    urlpath : string
        Absolute or relative filepath, URL (may include protocols like
        ``s3://``), or globstring pointing to data.
    mode : str, optional
        Mode in which to open files.
    num : int, optional
        If opening in writing mode, number of files we expect to create.
    name_function : callable, optional
        If opening in writing mode, this callable is used to generate path
        names. Names are generated for each partition by
        ``urlpath.replace('*', name_function(partition_index))``.
    storage_options : dict, optional
        Additional keywords to pass to the filesystem class.
    """
    if isinstance(urlpath, (list, tuple)):
        if not urlpath:
            raise ValueError("empty urlpath sequence")
        paths, protocols, options_list = zip(*map(infer_options, urlpath))
        protocol = protocols[0]
        options = options_list[0]
        if not (all(p == protocol for p in protocols) and
                all(o == options for o in options_list)):
            raise ValueError("When specifying a list of paths, all paths must "
                             "share the same protocol and options")
        update_storage_options(options, storage_options)
        paths = list(paths)

        fs, fs_token = get_fs(protocol, options)

    elif isinstance(urlpath, (str, unicode)) or hasattr(urlpath, 'name'):
        urlpath, protocol, options = infer_options(urlpath)
        update_storage_options(options, storage_options)

        fs, fs_token = get_fs(protocol, options)

        if 'w' in mode:
            paths = _expand_paths(urlpath, name_function, num)
        elif "*" in urlpath:
            paths = sorted(fs.glob(urlpath))
        else:
            paths = [urlpath]

    else:
        raise TypeError('url type not understood: %s' % urlpath)

    return fs, fs_token, paths


def open_text_files(urlpath, compression=None, mode='rt', encoding='utf8',
                    errors='strict', **kwargs):
    """ Given path return dask.delayed file-like objects in text mode

    This function is deprecated, use ``open_files(path, mode='rt', ...)``.

    Parameters
    ----------
    urlpath: string
        Absolute or relative filepath, URL (may include protocols like
        ``s3://``), or globstring pointing to data.
    encoding: string
    errors: string
    compression: string
        Compression to use.  See ``dask.bytes.compression.files`` for options.
    **kwargs: dict
        Extra options that make sense to a particular storage connection, e.g.
        host, port, username, password, etc.

    Examples
    --------
    >>> files = open_text_files('2015-*-*.csv', encoding='utf-8')  # doctest: +SKIP
    >>> files = open_text_files('s3://bucket/2015-*-*.csv')  # doctest: +SKIP

    Returns
    -------
    List of ``dask.delayed`` objects that compute to text file-like objects
    """
    warn("DeprecationWarning: open_text_files is deprecated, use `open_files` "
         "with mode='rt' or mode='wt'")
    return open_files(urlpath, mode=mode.replace('b', 't'),
                      compression=compression, encoding=encoding,
                      errors=errors, **kwargs)


def _expand_paths(path, name_function, num):
    if isinstance(path, (str, unicode)):
        if path.count('*') > 1:
            raise ValueError("Output path spec must contain at most one '*'.")
        elif '*' not in path:
            path = os.path.join(path, '*.part')

        if name_function is None:
            name_function = build_name_function(num - 1)

        paths = [path.replace('*', name_function(i)) for i in range(num)]
        if paths != sorted(paths):
            warn("In order to preserve order between partitions paths created "
                 "with ``name_function`` should sort to partition order")
    elif isinstance(path, (tuple, list)):
        assert len(path) == num
        paths = list(path)
    else:
        raise ValueError("Path should be either\n"
                         "1. A list of paths: ['foo.json', 'bar.json', ...]\n"
                         "2. A directory: 'foo/\n"
                         "3. A path with a '*' in it: 'foo.*.json'")
    return paths


def get_hdfs_driver(driver="auto"):
    """Get the hdfs driver implementation.

    Parameters
    ----------
    driver : {'auto', 'hdfs3', 'pyarrow'}, default 'auto'
        HDFS library to use. Default is first installed in this list.

    Returns
    -------
    A filesystem class
    """
    if driver == 'auto':
        for d in ['hdfs3', 'pyarrow']:
            try:
                return get_hdfs_driver(d)
            except RuntimeError:
                pass
        else:
            raise RuntimeError("Please install either `hdfs3` or `pyarrow`")

    elif driver == 'hdfs3':
        import_required('hdfs3', "`hdfs3` not installed")
        from dask.bytes.hdfs3 import HDFS3HadoopFileSystem as cls
        return cls

    elif driver == 'pyarrow':
        pa = import_required('pyarrow', "`pyarrow` not installed")
        from dask.bytes.pyarrow import (_MIN_PYARROW_VERSION_SUPPORTED,
                                        PyArrowHadoopFileSystem as cls)
        if LooseVersion(pa.__version__) < _MIN_PYARROW_VERSION_SUPPORTED:
            raise RuntimeError("pyarrow version >= %r required for hdfs driver "
                               "support" % _MIN_PYARROW_VERSION_SUPPORTED)
        return cls

    else:
        raise ValueError('Unsupported hdfs driver: {0}'.format(driver))


def get_fs(protocol, storage_options=None):
    """Create a filesystem object from a protocol and options.

    Parameters
    ----------
    protocol : str
        The specific protocol to use.
    storage_options : dict, optional
        Keywords to pass to the filesystem class.
    """
    if protocol in _filesystems:
        cls = _filesystems[protocol]

    elif protocol == 's3':
        import_required('s3fs',
                        "Need to install `s3fs` library for s3 support\n"
                        "    conda install s3fs -c conda-forge\n"
                        "    or\n"
                        "    pip install s3fs")
        import dask.bytes.s3  # noqa, register the s3 backend
        cls = _filesystems[protocol]

    elif protocol in ('gs', 'gcs'):
        import_required('gcsfs',
                        "Need to install `gcsfs` library for Google Cloud Storage support\n"
                        "    conda install gcsfs -c conda-forge\n"
                        "    or\n"
                        "    pip install gcsfs")
        cls = _filesystems[protocol]

    elif protocol == 'hdfs':
        cls = get_hdfs_driver(_globals.get("hdfs_driver", "auto"))

    elif protocol in ['http', 'https']:
        import_required('requests',
                        "Need to install `requests` for HTTP support\n"
                        "   conda install requests\n"
                        "    or\n"
                        "   pip install requests")
        import dask.bytes.http  # noqa, registers HTTP backend
        cls = _filesystems[protocol]

    else:
        raise ValueError("Unknown protocol %s" % protocol)

    if storage_options is None:
        storage_options = {}
    fs = cls(**storage_options)
    fs_token = tokenize(cls, protocol, storage_options)
    return fs, fs_token


_filesystems = dict()


class FileSystem(object):
    """Deprecated, do not use. Implement filesystems by matching the interface
    of `dask.bytes.local.LocalFileSystem` instead of subclassing."""
    pass


def logical_size(fs, path, compression='infer'):
    if compression == 'infer':
        compression = infer_compression(path)

    if compression is None:
        return fs.size(path)
    elif compression in seekable_files:
        with OpenFile(fs, path, compression=compression) as f:
            f.seek(0, 2)
            return f.tell()
    else:
        raise ValueError("Cannot infer logical size from file compressed with "
                         "compression=%r" % compression)


def get_pyarrow_filesystem(fs):
    """Get an equivalent pyarrow filesystem.

    Not for public use, will be removed once a consistent filesystem api
    is defined."""
    try:
        return fs._get_pyarrow_filesystem()
    except AttributeError:
        raise NotImplementedError("Using pyarrow with a %r "
                                  "filesystem object" % type(fs).__name__)
