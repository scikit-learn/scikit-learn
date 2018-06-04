from __future__ import print_function, division, absolute_import

import bz2
import sys
import zlib

from toolz import identity

from ..compatibility import gzip_compress, gzip_decompress, GzipFile
from ..utils import ignoring


def noop_file(file, **kwargs):
    return file


compress = {'gzip': gzip_compress,
            'zlib': zlib.compress,
            'bz2': bz2.compress,
            None: identity}
decompress = {'gzip': gzip_decompress,
              'zlib': zlib.decompress,
              'bz2': bz2.decompress,
              None: identity}
files = {'gzip': lambda f, **kwargs: GzipFile(fileobj=f, **kwargs),
         None: noop_file}
seekable_files = {None: noop_file}


with ignoring(ImportError):
    import snappy
    compress['snappy'] = snappy.compress
    decompress['snappy'] = snappy.decompress


try:
    import lz4.block
    compress['lz4'] = lz4.block.compress
    compress['lz4'] = lz4.block.decompress
except ImportError:
    try:
        import lz4
        compress['lz4'] = lz4.LZ4_compress
        compress['lz4'] = lz4.LZ4_uncompress
    except ImportError:
        pass

with ignoring(ImportError):
    from ..compatibility import LZMAFile, lzma_compress, lzma_decompress
    compress['xz'] = lzma_compress
    decompress['xz'] = lzma_decompress
    files['xz'] = LZMAFile

# Seekable xz files actually tend to scan whole file - see `get_xz_blocks`
# with ignoring(ImportError):
#     import lzma
#     seekable_files['xz'] = lzma.LZMAFile
#
# with ignoring(ImportError):
#     import lzmaffi
#     seekable_files['xz'] = lzmaffi.LZMAFile


if sys.version_info[0] >= 3:
    import bz2
    files['bz2'] = bz2.BZ2File


def get_xz_blocks(fp):
    from lzmaffi import (STREAM_HEADER_SIZE, decode_stream_footer,
                         decode_index, LZMAError)
    fp.seek(0, 2)

    def _peek(f, size):
        data = f.read(size)
        f.seek(-size, 1)
        return data

    if fp.tell() < 2 * STREAM_HEADER_SIZE:
        raise LZMAError("file too small")

    # read stream paddings (4 bytes each)
    fp.seek(-4, 1)
    padding = 0
    while _peek(fp, 4) == b'\x00\x00\x00\x00':
        fp.seek(-4, 1)
        padding += 4

    fp.seek(-STREAM_HEADER_SIZE + 4, 1)

    stream_flags = decode_stream_footer(_peek(fp, STREAM_HEADER_SIZE))
    fp.seek(-stream_flags.backward_size, 1)

    index = decode_index(_peek(fp, stream_flags.backward_size), padding)
    return {'offsets': [b.compressed_file_offset for i, b in index],
            'lengths': [b.unpadded_size for i, b in index],
            'check': stream_flags.check}


def xz_decompress(data, check):
    from lzmaffi import decode_block_header_size, LZMADecompressor, FORMAT_BLOCK
    hsize = decode_block_header_size(data[:1])
    header = data[:hsize]
    dc = LZMADecompressor(format=FORMAT_BLOCK, header=header,
                          unpadded_size=len(data), check=check)
    return dc.decompress(data[len(header):])
