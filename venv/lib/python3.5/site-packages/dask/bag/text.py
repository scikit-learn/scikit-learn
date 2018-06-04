from __future__ import print_function, division, absolute_import

import io
import os

from toolz import concat

from ..utils import system_encoding
from ..delayed import delayed
from ..bytes import open_files, read_bytes
from .core import from_delayed

delayed = delayed(pure=True)


def read_text(urlpath, blocksize=None, compression='infer',
              encoding=system_encoding, errors='strict',
              linedelimiter=os.linesep, collection=True,
              storage_options=None):
    """ Read lines from text files

    Parameters
    ----------
    urlpath : string or list
        Absolute or relative filepath(s). Prefix with a protocol like ``s3://``
        to read from alternative filesystems. To read from multiple files you
        can pass a globstring or a list of paths, with the caveat that they
        must all have the same protocol.
    blocksize: None or int
        Size (in bytes) to cut up larger files.  Streams by default.
    compression: string
        Compression format like 'gzip' or 'xz'.  Defaults to 'infer'
    encoding: string
    errors: string
    linedelimiter: string
    collection: bool, optional
        Return dask.bag if True, or list of delayed values if false
    storage_options: dict
        Extra options that make sense to a particular storage connection, e.g.
        host, port, username, password, etc.

    Examples
    --------
    >>> b = read_text('myfiles.1.txt')  # doctest: +SKIP
    >>> b = read_text('myfiles.*.txt')  # doctest: +SKIP
    >>> b = read_text('myfiles.*.txt.gz')  # doctest: +SKIP
    >>> b = read_text('s3://bucket/myfiles.*.txt')  # doctest: +SKIP
    >>> b = read_text('s3://key:secret@bucket/myfiles.*.txt')  # doctest: +SKIP
    >>> b = read_text('hdfs://namenode.example.com/myfiles.*.txt')  # doctest: +SKIP

    Parallelize a large file by providing the number of uncompressed bytes to
    load into each partition.

    >>> b = read_text('largefile.txt', blocksize=1e7)  # doctest: +SKIP

    Returns
    -------
    dask.bag.Bag if collection is True or list of Delayed lists otherwise

    See Also
    --------
    from_sequence: Build bag from Python sequence
    """
    if blocksize is None:
        files = open_files(urlpath, mode='rt', encoding=encoding,
                           errors=errors, compression=compression,
                           **(storage_options or {}))
        blocks = [delayed(list)(delayed(file_to_blocks)(fil)) for fil in files]

    else:
        _, blocks = read_bytes(urlpath, delimiter=linedelimiter.encode(),
                               blocksize=blocksize, sample=False,
                               compression=compression,
                               **(storage_options or {}))
        blocks = [delayed(decode)(b, encoding, errors) for b in concat(blocks)]

    if not blocks:
        raise ValueError("No files found", urlpath)

    if not collection:
        return blocks
    else:
        return from_delayed(blocks)


def file_to_blocks(lazy_file):
    with lazy_file as f:
        for line in f:
            yield line


def decode(block, encoding, errors):
    text = block.decode(encoding, errors)
    lines = io.StringIO(text)
    return list(lines)
