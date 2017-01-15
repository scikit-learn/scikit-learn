"""Utilities for fast persistence of big data, with optional compression."""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import pickle
import sys
import io
import zlib
import gzip
import bz2
import warnings
import contextlib
from contextlib import closing

from ._compat import PY3_OR_LATER, PY26, PY27, _basestring

try:
    from threading import RLock
except ImportError:
    from dummy_threading import RLock

if PY3_OR_LATER:
    Unpickler = pickle._Unpickler
    Pickler = pickle._Pickler
    xrange = range
else:
    Unpickler = pickle.Unpickler
    Pickler = pickle.Pickler

try:
    import numpy as np
except ImportError:
    np = None

try:
    import lzma
except ImportError:
    lzma = None


# Magic numbers of supported compression file formats.        '
_ZFILE_PREFIX = b'ZF'  # used with pickle files created before 0.9.3.
_ZLIB_PREFIX = b'\x78'
_GZIP_PREFIX = b'\x1f\x8b'
_BZ2_PREFIX = b'BZ'
_XZ_PREFIX = b'\xfd\x37\x7a\x58\x5a'
_LZMA_PREFIX = b'\x5d\x00'

# Supported compressors
_COMPRESSORS = ('zlib', 'bz2', 'lzma', 'xz', 'gzip')
_COMPRESSOR_CLASSES = [gzip.GzipFile, bz2.BZ2File]
if lzma is not None:
    _COMPRESSOR_CLASSES.append(lzma.LZMAFile)

# The max magic number length of supported compression file types.
_MAX_PREFIX_LEN = max(len(prefix)
                      for prefix in (_ZFILE_PREFIX, _GZIP_PREFIX, _BZ2_PREFIX,
                                     _XZ_PREFIX, _LZMA_PREFIX))

# Buffer size used in io.BufferedReader and io.BufferedWriter
_IO_BUFFER_SIZE = 1024 ** 2


###############################################################################
# Cache file utilities
def _detect_compressor(fileobj):
    """Return the compressor matching fileobj.

    Parameters
    ----------
    fileobj: file object

    Returns
    -------
    str in {'zlib', 'gzip', 'bz2', 'lzma', 'xz', 'compat', 'not-compressed'}
    """
    # Ensure we read the first bytes.
    fileobj.seek(0)
    first_bytes = fileobj.read(_MAX_PREFIX_LEN)
    fileobj.seek(0)

    if first_bytes.startswith(_ZLIB_PREFIX):
        return "zlib"
    elif first_bytes.startswith(_GZIP_PREFIX):
        return "gzip"
    elif first_bytes.startswith(_BZ2_PREFIX):
        return "bz2"
    elif first_bytes.startswith(_LZMA_PREFIX):
        return "lzma"
    elif first_bytes.startswith(_XZ_PREFIX):
        return "xz"
    elif first_bytes.startswith(_ZFILE_PREFIX):
        return "compat"

    return "not-compressed"


def _buffered_read_file(fobj):
    """Return a buffered version of a read file object."""
    if PY26 or (PY27 and isinstance(fobj, bz2.BZ2File)):
        # Python 2.6 doesn't fully support io.BufferedReader.
        # Python 2.7 doesn't work with BZ2File through a buffer: "no
        # attribute 'readable'" error.
        return fobj
    else:
        return io.BufferedReader(fobj, buffer_size=_IO_BUFFER_SIZE)


def _buffered_write_file(fobj):
    """Return a buffered version of a write file object."""
    if PY26 or (PY27 and isinstance(fobj, bz2.BZ2File)):
        # Python 2.6 doesn't fully support io.BufferedWriter.
        # Python 2.7 doesn't work with BZ2File through a buffer: no attribute
        # 'writable'.
        # BZ2File doesn't implement the file object context manager in python 2
        # so we wrap the fileobj using `closing`.
        return closing(fobj)
    else:
        return io.BufferedWriter(fobj, buffer_size=_IO_BUFFER_SIZE)


@contextlib.contextmanager
def _read_fileobject(fileobj, filename, mmap_mode=None):
    """Utility function opening the right fileobject from a filename.

    The magic number is used to choose between the type of file object to open:
    * regular file object (default)
    * zlib file object
    * gzip file object
    * bz2 file object
    * lzma file object (for xz and lzma compressor)

    Parameters
    ----------
    fileobj: file object
    compressor: str in {'zlib', 'gzip', 'bz2', 'lzma', 'xz', 'compat',
                        'not-compressed'}
    filename: str
        filename path corresponding to the fileobj parameter.
    mmap_mode: str
        memory map mode that should be used to open the pickle file. This
        parameter is useful to verify that the user is not trying to one with
        compression. Default: None.

    Returns
    -------
        a file like object

    """
    # Detect if the fileobj contains compressed data.
    compressor = _detect_compressor(fileobj)
    if isinstance(fileobj, tuple(_COMPRESSOR_CLASSES)):
        compressor = fileobj.__class__.__name__
    if compressor == 'compat':
        # Compatibility with old pickle mode: simply return the input
        # filename "as-is" and let the compatibility function be called by the
        # caller.
        warnings.warn("The file '%s' has been generated with a joblib "
                      "version less than 0.10. "
                      "Please regenerate this pickle file." % filename,
                      DeprecationWarning, stacklevel=2)
        yield filename
    else:
        # Checking if incompatible load parameters with the type of file:
        # mmap_mode cannot be used with compressed file or in memory buffers
        # such as io.BytesIO.
        if ((compressor in _COMPRESSORS or
                isinstance(fileobj, tuple(_COMPRESSOR_CLASSES))) and
                mmap_mode is not None):
            warnings.warn('File "%(filename)s" is compressed using '
                          '"%(compressor)s" which is not compatible with '
                          'mmap_mode "%(mmap_mode)s" flag passed. mmap_mode '
                          'option will be ignored.'
                          % locals(), stacklevel=2)
        if isinstance(fileobj, io.BytesIO) and mmap_mode is not None:
            warnings.warn('In memory persistence is not compatible with '
                          'mmap_mode "%(mmap_mode)s" flag passed. mmap_mode '
                          'option will be ignored.'
                          % locals(), stacklevel=2)

        # if the passed fileobj is in the supported list of decompressor
        # objects (GzipFile, BZ2File, LzmaFile), we simply return it.
        if isinstance(fileobj, tuple(_COMPRESSOR_CLASSES)):
            yield fileobj
        # otherwise, based on the compressor detected in the file, we open the
        # correct decompressor file object, wrapped in a buffer.
        elif compressor == 'zlib':
            yield _buffered_read_file(BinaryZlibFile(fileobj, 'rb'))
        elif compressor == 'gzip':
            yield _buffered_read_file(BinaryGzipFile(fileobj, 'rb'))
        elif compressor == 'bz2':
            if PY3_OR_LATER:
                yield _buffered_read_file(bz2.BZ2File(fileobj, 'rb'))
            else:
                # In python 2, BZ2File doesn't support a fileobj opened in
                # binary mode. In this case, we pass the filename.
                yield _buffered_read_file(bz2.BZ2File(fileobj.name, 'rb'))
        elif (compressor == 'lzma' or compressor == 'xz'):
            if lzma is not None:
                yield _buffered_read_file(lzma.LZMAFile(fileobj, 'rb'))
            else:
                raise NotImplementedError("Lzma decompression is not "
                                          "available for this version of "
                                          "python ({0}.{1})"
                                          .format(sys.version_info[0],
                                                  sys.version_info[1]))
        # No compression detected => returning the input file object (open)
        else:
            yield fileobj


def _write_fileobject(filename, compress=("zlib", 3)):
    """Return the right compressor file object in write mode."""
    compressmethod = compress[0]
    compresslevel = compress[1]
    if compressmethod == "gzip":
        return _buffered_write_file(BinaryGzipFile(filename, 'wb',
                                    compresslevel=compresslevel))
    elif compressmethod == "bz2":
        return _buffered_write_file(bz2.BZ2File(filename, 'wb',
                                                compresslevel=compresslevel))
    elif lzma is not None and compressmethod == "xz":
        return _buffered_write_file(lzma.LZMAFile(filename, 'wb',
                                                  check=lzma.CHECK_NONE,
                                                  preset=compresslevel))
    elif lzma is not None and compressmethod == "lzma":
        return _buffered_write_file(lzma.LZMAFile(filename, 'wb',
                                                  preset=compresslevel,
                                                  format=lzma.FORMAT_ALONE))
    else:
        return _buffered_write_file(BinaryZlibFile(filename, 'wb',
                                    compresslevel=compresslevel))


###############################################################################
#  Joblib zlib compression file object definition

_MODE_CLOSED = 0
_MODE_READ = 1
_MODE_READ_EOF = 2
_MODE_WRITE = 3
_BUFFER_SIZE = 8192


class BinaryZlibFile(io.BufferedIOBase):
    """A file object providing transparent zlib (de)compression.

    A BinaryZlibFile can act as a wrapper for an existing file object, or refer
    directly to a named file on disk.

    Note that BinaryZlibFile provides only a *binary* file interface: data read
    is returned as bytes, and data to be written should be given as bytes.

    This object is an adaptation of the BZ2File object and is compatible with
    versions of python >= 2.6.

    If filename is a str or bytes object, it gives the name
    of the file to be opened. Otherwise, it should be a file object,
    which will be used to read or write the compressed data.

    mode can be 'rb' for reading (default) or 'wb' for (over)writing

    If mode is 'wb', compresslevel can be a number between 1
    and 9 specifying the level of compression: 1 produces the least
    compression, and 9 (default) produces the most compression.
    """

    wbits = zlib.MAX_WBITS

    def __init__(self, filename, mode="rb", compresslevel=9):
        # This lock must be recursive, so that BufferedIOBase's
        # readline(), readlines() and writelines() don't deadlock.
        self._lock = RLock()
        self._fp = None
        self._closefp = False
        self._mode = _MODE_CLOSED
        self._pos = 0
        self._size = -1

        if not isinstance(compresslevel, int) or not (1 <= compresslevel <= 9):
            raise ValueError("compresslevel must be between an integer "
                             "between 1 and 9, you gave {0}"
                             .format(compresslevel))

        if mode == "rb":
            mode_code = _MODE_READ
            self._decompressor = zlib.decompressobj(self.wbits)
            self._buffer = b""
            self._buffer_offset = 0
        elif mode == "wb":
            mode_code = _MODE_WRITE
            self._compressor = zlib.compressobj(compresslevel,
                                                zlib.DEFLATED,
                                                self.wbits,
                                                zlib.DEF_MEM_LEVEL,
                                                0)
        else:
            raise ValueError("Invalid mode: %r" % (mode,))

        if isinstance(filename, _basestring):
            self._fp = open(filename, mode)
            self._closefp = True
            self._mode = mode_code
        elif hasattr(filename, "read") or hasattr(filename, "write"):
            self._fp = filename
            self._mode = mode_code
        else:
            raise TypeError("filename must be a str or bytes object, "
                            "or a file")

    def close(self):
        """Flush and close the file.

        May be called more than once without error. Once the file is
        closed, any other operation on it will raise a ValueError.
        """
        with self._lock:
            if self._mode == _MODE_CLOSED:
                return
            try:
                if self._mode in (_MODE_READ, _MODE_READ_EOF):
                    self._decompressor = None
                elif self._mode == _MODE_WRITE:
                    self._fp.write(self._compressor.flush())
                    self._compressor = None
            finally:
                try:
                    if self._closefp:
                        self._fp.close()
                finally:
                    self._fp = None
                    self._closefp = False
                    self._mode = _MODE_CLOSED
                    self._buffer = b""
                    self._buffer_offset = 0

    @property
    def closed(self):
        """True if this file is closed."""
        return self._mode == _MODE_CLOSED

    def fileno(self):
        """Return the file descriptor for the underlying file."""
        self._check_not_closed()
        return self._fp.fileno()

    def seekable(self):
        """Return whether the file supports seeking."""
        return self.readable() and self._fp.seekable()

    def readable(self):
        """Return whether the file was opened for reading."""
        self._check_not_closed()
        return self._mode in (_MODE_READ, _MODE_READ_EOF)

    def writable(self):
        """Return whether the file was opened for writing."""
        self._check_not_closed()
        return self._mode == _MODE_WRITE

    # Mode-checking helper functions.

    def _check_not_closed(self):
        if self.closed:
            fname = getattr(self._fp, 'name', None)
            msg = "I/O operation on closed file"
            if fname is not None:
                msg += " {0}".format(fname)
            msg += "."
            raise ValueError(msg)

    def _check_can_read(self):
        if self._mode not in (_MODE_READ, _MODE_READ_EOF):
            self._check_not_closed()
            raise io.UnsupportedOperation("File not open for reading")

    def _check_can_write(self):
        if self._mode != _MODE_WRITE:
            self._check_not_closed()
            raise io.UnsupportedOperation("File not open for writing")

    def _check_can_seek(self):
        if self._mode not in (_MODE_READ, _MODE_READ_EOF):
            self._check_not_closed()
            raise io.UnsupportedOperation("Seeking is only supported "
                                          "on files open for reading")
        if not self._fp.seekable():
            raise io.UnsupportedOperation("The underlying file object "
                                          "does not support seeking")

    # Fill the readahead buffer if it is empty. Returns False on EOF.
    def _fill_buffer(self):
        if self._mode == _MODE_READ_EOF:
            return False
        # Depending on the input data, our call to the decompressor may not
        # return any data. In this case, try again after reading another block.
        while self._buffer_offset == len(self._buffer):
            try:
                rawblock = (self._decompressor.unused_data or
                            self._fp.read(_BUFFER_SIZE))

                if not rawblock:
                    raise EOFError
            except EOFError:
                # End-of-stream marker and end of file. We're good.
                self._mode = _MODE_READ_EOF
                self._size = self._pos
                return False
            else:
                self._buffer = self._decompressor.decompress(rawblock)
            self._buffer_offset = 0
        return True

    # Read data until EOF.
    # If return_data is false, consume the data without returning it.
    def _read_all(self, return_data=True):
        # The loop assumes that _buffer_offset is 0. Ensure that this is true.
        self._buffer = self._buffer[self._buffer_offset:]
        self._buffer_offset = 0

        blocks = []
        while self._fill_buffer():
            if return_data:
                blocks.append(self._buffer)
            self._pos += len(self._buffer)
            self._buffer = b""
        if return_data:
            return b"".join(blocks)

    # Read a block of up to n bytes.
    # If return_data is false, consume the data without returning it.
    def _read_block(self, n_bytes, return_data=True):
        # If we have enough data buffered, return immediately.
        end = self._buffer_offset + n_bytes
        if end <= len(self._buffer):
            data = self._buffer[self._buffer_offset: end]
            self._buffer_offset = end
            self._pos += len(data)
            return data if return_data else None

        # The loop assumes that _buffer_offset is 0. Ensure that this is true.
        self._buffer = self._buffer[self._buffer_offset:]
        self._buffer_offset = 0

        blocks = []
        while n_bytes > 0 and self._fill_buffer():
            if n_bytes < len(self._buffer):
                data = self._buffer[:n_bytes]
                self._buffer_offset = n_bytes
            else:
                data = self._buffer
                self._buffer = b""
            if return_data:
                blocks.append(data)
            self._pos += len(data)
            n_bytes -= len(data)
        if return_data:
            return b"".join(blocks)

    def read(self, size=-1):
        """Read up to size uncompressed bytes from the file.

        If size is negative or omitted, read until EOF is reached.
        Returns b'' if the file is already at EOF.
        """
        with self._lock:
            self._check_can_read()
            if size == 0:
                return b""
            elif size < 0:
                return self._read_all()
            else:
                return self._read_block(size)

    def readinto(self, b):
        """Read up to len(b) bytes into b.

        Returns the number of bytes read (0 for EOF).
        """
        with self._lock:
            return io.BufferedIOBase.readinto(self, b)

    def write(self, data):
        """Write a byte string to the file.

        Returns the number of uncompressed bytes written, which is
        always len(data). Note that due to buffering, the file on disk
        may not reflect the data written until close() is called.
        """
        with self._lock:
            self._check_can_write()
            # Convert data type if called by io.BufferedWriter.
            if not PY26 and isinstance(data, memoryview):
                data = data.tobytes()

            compressed = self._compressor.compress(data)
            self._fp.write(compressed)
            self._pos += len(data)
            return len(data)

    # Rewind the file to the beginning of the data stream.
    def _rewind(self):
        self._fp.seek(0, 0)
        self._mode = _MODE_READ
        self._pos = 0
        self._decompressor = zlib.decompressobj(self.wbits)
        self._buffer = b""
        self._buffer_offset = 0

    def seek(self, offset, whence=0):
        """Change the file position.

        The new position is specified by offset, relative to the
        position indicated by whence. Values for whence are:

            0: start of stream (default); offset must not be negative
            1: current stream position
            2: end of stream; offset must not be positive

        Returns the new file position.

        Note that seeking is emulated, so depending on the parameters,
        this operation may be extremely slow.
        """
        with self._lock:
            self._check_can_seek()

            # Recalculate offset as an absolute file position.
            if whence == 0:
                pass
            elif whence == 1:
                offset = self._pos + offset
            elif whence == 2:
                # Seeking relative to EOF - we need to know the file's size.
                if self._size < 0:
                    self._read_all(return_data=False)
                offset = self._size + offset
            else:
                raise ValueError("Invalid value for whence: %s" % (whence,))

            # Make it so that offset is the number of bytes to skip forward.
            if offset < self._pos:
                self._rewind()
            else:
                offset -= self._pos

            # Read and discard data until we reach the desired position.
            self._read_block(offset, return_data=False)

            return self._pos

    def tell(self):
        """Return the current file position."""
        with self._lock:
            self._check_not_closed()
            return self._pos


class BinaryGzipFile(BinaryZlibFile):
    """A file object providing transparent gzip (de)compression.

    If filename is a str or bytes object, it gives the name
    of the file to be opened. Otherwise, it should be a file object,
    which will be used to read or write the compressed data.

    mode can be 'rb' for reading (default) or 'wb' for (over)writing

    If mode is 'wb', compresslevel can be a number between 1
    and 9 specifying the level of compression: 1 produces the least
    compression, and 9 (default) produces the most compression.
    """

    wbits = 31  # zlib compressor/decompressor wbits value for gzip format.


# Utility functions/variables from numpy required for writing arrays.
# We need at least the functions introduced in version 1.9 of numpy. Here,
# we use the ones from numpy 1.10.2.
BUFFER_SIZE = 2 ** 18  # size of buffer for reading npz files in bytes


def _read_bytes(fp, size, error_template="ran out of data"):
    """Read from file-like object until size bytes are read.

    Raises ValueError if not EOF is encountered before size bytes are read.
    Non-blocking objects only supported if they derive from io objects.

    Required as e.g. ZipExtFile in python 2.6 can return less data than
    requested.

    This function was taken from numpy/lib/format.py in version 1.10.2.

    Parameters
    ----------
    fp: file-like object
    size: int
    error_template: str

    Returns
    -------
    a bytes object
        The data read in bytes.

    """
    data = bytes()
    while True:
        # io files (default in python3) return None or raise on
        # would-block, python2 file will truncate, probably nothing can be
        # done about that.  note that regular files can't be non-blocking
        try:
            r = fp.read(size - len(data))
            data += r
            if len(r) == 0 or len(data) == size:
                break
        except io.BlockingIOError:
            pass
    if len(data) != size:
        msg = "EOF: reading %s, expected %d bytes got %d"
        raise ValueError(msg % (error_template, size, len(data)))
    else:
        return data
