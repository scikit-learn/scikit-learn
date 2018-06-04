from __future__ import print_function, division, absolute_import

import requests
import uuid

from . import core

DEFAULT_BLOCK_SIZE = 5 * 2 ** 20


class HTTPFileSystem(core.FileSystem):
    """
    Simple File-System for fetching data via HTTP(S)

    Unlike other file-systems, HTTP is limited in that it does not provide glob
    or write capability.
    """
    sep = '/'

    def __init__(self, **storage_options):
        """
        Parameters
        ----------
        storage_options: key-value
            May be credentials, e.g., `{'auth': ('username', 'pword')}` or any
            other parameters for requests
        """
        self.block_size = storage_options.pop('block_size', DEFAULT_BLOCK_SIZE)
        self.kwargs = storage_options
        self.session = requests.Session()

    def glob(self, url):
        """For a template path, return matching files"""
        raise NotImplementedError

    def mkdirs(self, url):
        """Make any intermediate directories to make path writable"""
        raise NotImplementedError

    def open(self, url, mode='rb', block_size=None, **kwargs):
        """Make a file-like object

        Parameters
        ----------
        url: str
            Full URL with protocol
        mode: string
            must be "rb"
        kwargs: key-value
            Any other parameters, passed to requests calls
        """
        if mode != 'rb':
            raise NotImplementedError
        block_size = block_size if block_size is not None else self.block_size
        return HTTPFile(url, self.session, block_size, **self.kwargs)

    def ukey(self, url):
        """Unique identifier, implied file might have changed every time"""
        return uuid.uuid1().hex

    def size(self, url):
        """Size in bytes of the file at path"""
        return file_size(url, session=self.session, **self.kwargs)


core._filesystems['http'] = HTTPFileSystem
core._filesystems['https'] = HTTPFileSystem


class HTTPFile(object):
    """
    A file-like object pointing to a remove HTTP(S) resource

    Supports only reading, with read-ahead of a predermined block-size.

    In the case that the server does not supply the filesize, only reading of
    the complete file in one go is supported.

    Parameters
    ----------
    url: str
        Full URL of the remote resource, including the protocol
    session: requests.Session or None
        All calls will be made within this session, to avoid restarting
        connections where the server allows this
    block_size: int or None
        The amount of read-ahead to do, in bytes. Default is 5MB, or the value
        configured for the FileSystem creating this file
    kwargs: all other key-values are passed to reqeuests calls.
    """

    def __init__(self, url, session=None, block_size=None, **kwargs):
        self.url = url
        self.kwargs = kwargs
        self.loc = 0
        self.session = session if session is not None else requests.Session()
        self.blocksize = (block_size if block_size is not None
                          else DEFAULT_BLOCK_SIZE)
        try:
            self.size = file_size(url, self.session, allow_redirects=True,
                                  **self.kwargs)
        except ValueError:
            # No size information - only allow read() and no seek()
            self.size = None
        self.cache = None
        self.closed = False
        self.start = None
        self.end = None

    def seek(self, where, whence=0):
        """Set file position

        Parameters
        ----------
        where: int
            Location to set
        whence: int (default 0)
            If zero, set from start of file (value should be positive); if 1,
            set relative to current position; if 2, set relative to end of file
            (value shoulf be negative)

        Returns the position.
        """
        if self.size is None:
            raise ValueError('Cannot seek since size of file is not known')
        if whence == 0:
            nloc = where
        elif whence == 1:
            nloc += where
        elif whence == 2:
            nloc = self.size + where
        else:
            raise ValueError('Whence must be in [1, 2, 3], but got %s' % whence)
        if nloc < 0:
            raise ValueError('Seek before start of file')
        self.loc = nloc
        return nloc

    def tell(self):
        """Get current file byte position"""
        return self.loc

    def read(self, length=-1):
        """Read bytes from file

        Parameters
        ----------
        length: int
            Read up to this many bytes. If negative, read all content to end of
            file. If the server has not supplied the filesize, attempting to
            read only part of the data will raise a ValueError.
        """
        if self.size is None:
            if length >= 0:
                raise ValueError('File size is unknown, must read all data')
            else:
                return self._fetch_all()
        if length < 0 or self.loc + length > self.size:
            end = self.size
        else:
            end = self.loc + length
        if self.loc >= self.size:
            return b''
        self. _fetch(self.loc, end)
        data = self.cache[self.loc - self.start:end - self.start]
        self.loc = end
        return data

    def _fetch(self, start, end):
        """Set new bounds for data cache and fetch data, if required"""
        if self.start is None and self.end is None:
            # First read
            self.start = start
            self.end = end + self.blocksize
            self.cache = self._fetch_range(start, self.end)
        elif start < self.start:
            if self.end - end > self.blocksize:
                self.start = start
                self.end = end + self.blocksize
                self.cache = self._fetch_range(self.start, self.end)
            else:
                new = self._fetch_range(start, self.start)
                self.start = start
                self.cache = new + self.cache
        elif end > self.end:
            if self.end > self.size:
                return
            if end - self.end > self.blocksize:
                self.start = start
                self.end = end + self.blocksize
                self.cache = self._fetch_range(self.start, self.end)
            else:
                new = self._fetch_range(self.end, end + self.blocksize)
                self.end = end + self.blocksize
                self.cache = self.cache + new

    def _fetch_all(self):
        """Read whole file in one shot, without caching

        This is only called when size is None and read() is called without a
        byte-count.
        """
        r = self.session.get(self.url, **self.kwargs)
        r.raise_for_status()
        return r.content

    def _fetch_range(self, start, end):
        """Download a block of data

        The expectation is that the server returns only the requested bytes,
        with HTTP code 206. If this is not the case, we first check the headers,
        and then stream the output - if the data size is bigger than we
        requested, an exception is raised.
        """
        kwargs = self.kwargs.copy()
        headers = self.kwargs.pop('headers', {})
        headers['Range'] = 'bytes=%i-%i' % (start, end - 1)
        r = self.session.get(self.url, headers=headers, stream=True, **kwargs)
        r.raise_for_status()
        if r.status_code == 206:
            # partial content, as expected
            return r.content
        if 'Content-Length' in r.headers:
            cl = int(r.headers['Content-Length'])
            if cl <= end - start:
                # data size OK
                return r.content
            else:
                raise ValueError('Got more bytes (%i) than requested (%i)' % (
                    cl, end - start))
        cl = 0
        out = []
        for chunk in r.iter_content(chunk_size=2 ** 20):
            # data size unknown, let's see if it goes too big
            if chunk:
                out.append(chunk)
                cl += len(chunk)
                if cl > end - start:
                    raise ValueError(
                        'Got more bytes so far (>%i) than requested (%i)' % (
                            cl, end - start))
            else:
                break
        return b''.join(out)

    def __enter__(self):
        self.loc = 0
        return self

    def __exit__(self, *args):
        self.close()

    def __iter__(self):
        # no text lines here, use TextIOWrapper
        raise NotImplementedError

    def write(self):
        raise NotImplementedError

    def flush(self):
        pass

    def close(self):
        self.closed = True

    def seekable(self):
        return True

    def writable(self):
        return False

    def readable(self):
        return True


def file_size(url, session, **kwargs):
    """Call HEAD on the server to get file size"""
    r = session.head(url, **kwargs)
    r.raise_for_status()
    if 'Content-Length' in r.headers:
        return int(r.headers['Content-Length'])
    else:
        raise ValueError("Server did not supply size of %s" % url)
