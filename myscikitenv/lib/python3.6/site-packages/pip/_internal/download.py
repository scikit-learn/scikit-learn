from __future__ import absolute_import

import cgi
import email.utils
import getpass
import json
import logging
import mimetypes
import os
import platform
import re
import shutil
import sys

from pip._vendor import requests, six, urllib3
from pip._vendor.cachecontrol import CacheControlAdapter
from pip._vendor.cachecontrol.caches import FileCache
from pip._vendor.lockfile import LockError
from pip._vendor.requests.adapters import BaseAdapter, HTTPAdapter
from pip._vendor.requests.auth import AuthBase, HTTPBasicAuth
from pip._vendor.requests.models import CONTENT_CHUNK_SIZE, Response
from pip._vendor.requests.structures import CaseInsensitiveDict
from pip._vendor.requests.utils import get_netrc_auth
# NOTE: XMLRPC Client is not annotated in typeshed as on 2017-07-17, which is
#       why we ignore the type on this import
from pip._vendor.six.moves import xmlrpc_client  # type: ignore
from pip._vendor.six.moves.urllib import parse as urllib_parse
from pip._vendor.six.moves.urllib import request as urllib_request
from pip._vendor.urllib3.util import IS_PYOPENSSL

import pip
from pip._internal.exceptions import HashMismatch, InstallationError
from pip._internal.locations import write_delete_marker_file
from pip._internal.models.index import PyPI
from pip._internal.utils.encoding import auto_decode
from pip._internal.utils.filesystem import check_path_owner
from pip._internal.utils.glibc import libc_ver
from pip._internal.utils.misc import (
    ARCHIVE_EXTENSIONS, ask_path_exists, backup_dir, consume, display_path,
    format_size, get_installed_version, rmtree, split_auth_from_netloc,
    splitext, unpack_file,
)
from pip._internal.utils.temp_dir import TempDirectory
from pip._internal.utils.typing import MYPY_CHECK_RUNNING
from pip._internal.utils.ui import DownloadProgressProvider
from pip._internal.vcs import vcs

if MYPY_CHECK_RUNNING:
    from typing import (
        Optional, Tuple, Dict, IO, Text, Union
    )
    from pip._internal.models.link import Link
    from pip._internal.utils.hashes import Hashes
    from pip._internal.vcs import AuthInfo

try:
    import ssl  # noqa
except ImportError:
    ssl = None

HAS_TLS = (ssl is not None) or IS_PYOPENSSL

__all__ = ['get_file_content',
           'is_url', 'url_to_path', 'path_to_url',
           'is_archive_file', 'unpack_vcs_link',
           'unpack_file_url', 'is_vcs_url', 'is_file_url',
           'unpack_http_url', 'unpack_url']


logger = logging.getLogger(__name__)


# These are environment variables present when running under various
# CI systems.  For each variable, some CI systems that use the variable
# are indicated.  The collection was chosen so that for each of a number
# of popular systems, at least one of the environment variables is used.
# This list is used to provide some indication of and lower bound for
# CI traffic to PyPI.  Thus, it is okay if the list is not comprehensive.
# For more background, see: https://github.com/pypa/pip/issues/5499
CI_ENVIRONMENT_VARIABLES = (
    # Azure Pipelines
    'BUILD_BUILDID',
    # Jenkins
    'BUILD_ID',
    # AppVeyor, CircleCI, Codeship, Gitlab CI, Shippable, Travis CI
    'CI',
)


def looks_like_ci():
    # type: () -> bool
    """
    Return whether it looks like pip is running under CI.
    """
    # We don't use the method of checking for a tty (e.g. using isatty())
    # because some CI systems mimic a tty (e.g. Travis CI).  Thus that
    # method doesn't provide definitive information in either direction.
    return any(name in os.environ for name in CI_ENVIRONMENT_VARIABLES)


def user_agent():
    """
    Return a string representing the user agent.
    """
    data = {
        "installer": {"name": "pip", "version": pip.__version__},
        "python": platform.python_version(),
        "implementation": {
            "name": platform.python_implementation(),
        },
    }

    if data["implementation"]["name"] == 'CPython':
        data["implementation"]["version"] = platform.python_version()
    elif data["implementation"]["name"] == 'PyPy':
        if sys.pypy_version_info.releaselevel == 'final':
            pypy_version_info = sys.pypy_version_info[:3]
        else:
            pypy_version_info = sys.pypy_version_info
        data["implementation"]["version"] = ".".join(
            [str(x) for x in pypy_version_info]
        )
    elif data["implementation"]["name"] == 'Jython':
        # Complete Guess
        data["implementation"]["version"] = platform.python_version()
    elif data["implementation"]["name"] == 'IronPython':
        # Complete Guess
        data["implementation"]["version"] = platform.python_version()

    if sys.platform.startswith("linux"):
        from pip._vendor import distro
        distro_infos = dict(filter(
            lambda x: x[1],
            zip(["name", "version", "id"], distro.linux_distribution()),
        ))
        libc = dict(filter(
            lambda x: x[1],
            zip(["lib", "version"], libc_ver()),
        ))
        if libc:
            distro_infos["libc"] = libc
        if distro_infos:
            data["distro"] = distro_infos

    if sys.platform.startswith("darwin") and platform.mac_ver()[0]:
        data["distro"] = {"name": "macOS", "version": platform.mac_ver()[0]}

    if platform.system():
        data.setdefault("system", {})["name"] = platform.system()

    if platform.release():
        data.setdefault("system", {})["release"] = platform.release()

    if platform.machine():
        data["cpu"] = platform.machine()

    if HAS_TLS:
        data["openssl_version"] = ssl.OPENSSL_VERSION

    setuptools_version = get_installed_version("setuptools")
    if setuptools_version is not None:
        data["setuptools_version"] = setuptools_version

    # Use None rather than False so as not to give the impression that
    # pip knows it is not being run under CI.  Rather, it is a null or
    # inconclusive result.  Also, we include some value rather than no
    # value to make it easier to know that the check has been run.
    data["ci"] = True if looks_like_ci() else None

    user_data = os.environ.get("PIP_USER_AGENT_USER_DATA")
    if user_data is not None:
        data["user_data"] = user_data

    return "{data[installer][name]}/{data[installer][version]} {json}".format(
        data=data,
        json=json.dumps(data, separators=(",", ":"), sort_keys=True),
    )


class MultiDomainBasicAuth(AuthBase):

    def __init__(self, prompting=True):
        # type: (bool) -> None
        self.prompting = prompting
        self.passwords = {}  # type: Dict[str, AuthInfo]

    def __call__(self, req):
        parsed = urllib_parse.urlparse(req.url)

        # Split the credentials from the netloc.
        netloc, url_user_password = split_auth_from_netloc(parsed.netloc)

        # Set the url of the request to the url without any credentials
        req.url = urllib_parse.urlunparse(parsed[:1] + (netloc,) + parsed[2:])

        # Use any stored credentials that we have for this netloc
        username, password = self.passwords.get(netloc, (None, None))

        # Use the credentials embedded in the url if we have none stored
        if username is None:
            username, password = url_user_password

        # Get creds from netrc if we still don't have them
        if username is None and password is None:
            netrc_auth = get_netrc_auth(req.url)
            username, password = netrc_auth if netrc_auth else (None, None)

        if username or password:
            # Store the username and password
            self.passwords[netloc] = (username, password)

            # Send the basic auth with this request
            req = HTTPBasicAuth(username or "", password or "")(req)

        # Attach a hook to handle 401 responses
        req.register_hook("response", self.handle_401)

        return req

    def handle_401(self, resp, **kwargs):
        # We only care about 401 responses, anything else we want to just
        #   pass through the actual response
        if resp.status_code != 401:
            return resp

        # We are not able to prompt the user so simply return the response
        if not self.prompting:
            return resp

        parsed = urllib_parse.urlparse(resp.url)

        # Prompt the user for a new username and password
        username = six.moves.input("User for %s: " % parsed.netloc)
        password = getpass.getpass("Password: ")

        # Store the new username and password to use for future requests
        if username or password:
            self.passwords[parsed.netloc] = (username, password)

        # Consume content and release the original connection to allow our new
        #   request to reuse the same one.
        resp.content
        resp.raw.release_conn()

        # Add our new username and password to the request
        req = HTTPBasicAuth(username or "", password or "")(resp.request)
        req.register_hook("response", self.warn_on_401)

        # Send our new request
        new_resp = resp.connection.send(req, **kwargs)
        new_resp.history.append(resp)

        return new_resp

    def warn_on_401(self, resp, **kwargs):
        # warn user that they provided incorrect credentials
        if resp.status_code == 401:
            logger.warning('401 Error, Credentials not correct for %s',
                           resp.request.url)


class LocalFSAdapter(BaseAdapter):

    def send(self, request, stream=None, timeout=None, verify=None, cert=None,
             proxies=None):
        pathname = url_to_path(request.url)

        resp = Response()
        resp.status_code = 200
        resp.url = request.url

        try:
            stats = os.stat(pathname)
        except OSError as exc:
            resp.status_code = 404
            resp.raw = exc
        else:
            modified = email.utils.formatdate(stats.st_mtime, usegmt=True)
            content_type = mimetypes.guess_type(pathname)[0] or "text/plain"
            resp.headers = CaseInsensitiveDict({
                "Content-Type": content_type,
                "Content-Length": stats.st_size,
                "Last-Modified": modified,
            })

            resp.raw = open(pathname, "rb")
            resp.close = resp.raw.close

        return resp

    def close(self):
        pass


class SafeFileCache(FileCache):
    """
    A file based cache which is safe to use even when the target directory may
    not be accessible or writable.
    """

    def __init__(self, *args, **kwargs):
        super(SafeFileCache, self).__init__(*args, **kwargs)

        # Check to ensure that the directory containing our cache directory
        # is owned by the user current executing pip. If it does not exist
        # we will check the parent directory until we find one that does exist.
        # If it is not owned by the user executing pip then we will disable
        # the cache and log a warning.
        if not check_path_owner(self.directory):
            logger.warning(
                "The directory '%s' or its parent directory is not owned by "
                "the current user and the cache has been disabled. Please "
                "check the permissions and owner of that directory. If "
                "executing pip with sudo, you may want sudo's -H flag.",
                self.directory,
            )

            # Set our directory to None to disable the Cache
            self.directory = None

    def get(self, *args, **kwargs):
        # If we don't have a directory, then the cache should be a no-op.
        if self.directory is None:
            return

        try:
            return super(SafeFileCache, self).get(*args, **kwargs)
        except (LockError, OSError, IOError):
            # We intentionally silence this error, if we can't access the cache
            # then we can just skip caching and process the request as if
            # caching wasn't enabled.
            pass

    def set(self, *args, **kwargs):
        # If we don't have a directory, then the cache should be a no-op.
        if self.directory is None:
            return

        try:
            return super(SafeFileCache, self).set(*args, **kwargs)
        except (LockError, OSError, IOError):
            # We intentionally silence this error, if we can't access the cache
            # then we can just skip caching and process the request as if
            # caching wasn't enabled.
            pass

    def delete(self, *args, **kwargs):
        # If we don't have a directory, then the cache should be a no-op.
        if self.directory is None:
            return

        try:
            return super(SafeFileCache, self).delete(*args, **kwargs)
        except (LockError, OSError, IOError):
            # We intentionally silence this error, if we can't access the cache
            # then we can just skip caching and process the request as if
            # caching wasn't enabled.
            pass


class InsecureHTTPAdapter(HTTPAdapter):

    def cert_verify(self, conn, url, verify, cert):
        conn.cert_reqs = 'CERT_NONE'
        conn.ca_certs = None


class PipSession(requests.Session):

    timeout = None  # type: Optional[int]

    def __init__(self, *args, **kwargs):
        retries = kwargs.pop("retries", 0)
        cache = kwargs.pop("cache", None)
        insecure_hosts = kwargs.pop("insecure_hosts", [])

        super(PipSession, self).__init__(*args, **kwargs)

        # Attach our User Agent to the request
        self.headers["User-Agent"] = user_agent()

        # Attach our Authentication handler to the session
        self.auth = MultiDomainBasicAuth()

        # Create our urllib3.Retry instance which will allow us to customize
        # how we handle retries.
        retries = urllib3.Retry(
            # Set the total number of retries that a particular request can
            # have.
            total=retries,

            # A 503 error from PyPI typically means that the Fastly -> Origin
            # connection got interrupted in some way. A 503 error in general
            # is typically considered a transient error so we'll go ahead and
            # retry it.
            # A 500 may indicate transient error in Amazon S3
            # A 520 or 527 - may indicate transient error in CloudFlare
            status_forcelist=[500, 503, 520, 527],

            # Add a small amount of back off between failed requests in
            # order to prevent hammering the service.
            backoff_factor=0.25,
        )

        # We want to _only_ cache responses on securely fetched origins. We do
        # this because we can't validate the response of an insecurely fetched
        # origin, and we don't want someone to be able to poison the cache and
        # require manual eviction from the cache to fix it.
        if cache:
            secure_adapter = CacheControlAdapter(
                cache=SafeFileCache(cache, use_dir_lock=True),
                max_retries=retries,
            )
        else:
            secure_adapter = HTTPAdapter(max_retries=retries)

        # Our Insecure HTTPAdapter disables HTTPS validation. It does not
        # support caching (see above) so we'll use it for all http:// URLs as
        # well as any https:// host that we've marked as ignoring TLS errors
        # for.
        insecure_adapter = InsecureHTTPAdapter(max_retries=retries)

        self.mount("https://", secure_adapter)
        self.mount("http://", insecure_adapter)

        # Enable file:// urls
        self.mount("file://", LocalFSAdapter())

        # We want to use a non-validating adapter for any requests which are
        # deemed insecure.
        for host in insecure_hosts:
            self.mount("https://{}/".format(host), insecure_adapter)

    def request(self, method, url, *args, **kwargs):
        # Allow setting a default timeout on a session
        kwargs.setdefault("timeout", self.timeout)

        # Dispatch the actual request
        return super(PipSession, self).request(method, url, *args, **kwargs)


def get_file_content(url, comes_from=None, session=None):
    # type: (str, Optional[str], Optional[PipSession]) -> Tuple[str, Text]
    """Gets the content of a file; it may be a filename, file: URL, or
    http: URL.  Returns (location, content).  Content is unicode.

    :param url:         File path or url.
    :param comes_from:  Origin description of requirements.
    :param session:     Instance of pip.download.PipSession.
    """
    if session is None:
        raise TypeError(
            "get_file_content() missing 1 required keyword argument: 'session'"
        )

    match = _scheme_re.search(url)
    if match:
        scheme = match.group(1).lower()
        if (scheme == 'file' and comes_from and
                comes_from.startswith('http')):
            raise InstallationError(
                'Requirements file %s references URL %s, which is local'
                % (comes_from, url))
        if scheme == 'file':
            path = url.split(':', 1)[1]
            path = path.replace('\\', '/')
            match = _url_slash_drive_re.match(path)
            if match:
                path = match.group(1) + ':' + path.split('|', 1)[1]
            path = urllib_parse.unquote(path)
            if path.startswith('/'):
                path = '/' + path.lstrip('/')
            url = path
        else:
            # FIXME: catch some errors
            resp = session.get(url)
            resp.raise_for_status()
            return resp.url, resp.text
    try:
        with open(url, 'rb') as f:
            content = auto_decode(f.read())
    except IOError as exc:
        raise InstallationError(
            'Could not open requirements file: %s' % str(exc)
        )
    return url, content


_scheme_re = re.compile(r'^(http|https|file):', re.I)
_url_slash_drive_re = re.compile(r'/*([a-z])\|', re.I)


def is_url(name):
    # type: (Union[str, Text]) -> bool
    """Returns true if the name looks like a URL"""
    if ':' not in name:
        return False
    scheme = name.split(':', 1)[0].lower()
    return scheme in ['http', 'https', 'file', 'ftp'] + vcs.all_schemes


def url_to_path(url):
    # type: (str) -> str
    """
    Convert a file: URL to a path.
    """
    assert url.startswith('file:'), (
        "You can only turn file: urls into filenames (not %r)" % url)

    _, netloc, path, _, _ = urllib_parse.urlsplit(url)

    if not netloc or netloc == 'localhost':
        # According to RFC 8089, same as empty authority.
        netloc = ''
    elif sys.platform == 'win32':
        # If we have a UNC path, prepend UNC share notation.
        netloc = '\\\\' + netloc
    else:
        raise ValueError(
            'non-local file URIs are not supported on this platform: %r'
            % url
        )

    path = urllib_request.url2pathname(netloc + path)
    return path


def path_to_url(path):
    # type: (Union[str, Text]) -> str
    """
    Convert a path to a file: URL.  The path will be made absolute and have
    quoted path parts.
    """
    path = os.path.normpath(os.path.abspath(path))
    url = urllib_parse.urljoin('file:', urllib_request.pathname2url(path))
    return url


def is_archive_file(name):
    # type: (str) -> bool
    """Return True if `name` is a considered as an archive file."""
    ext = splitext(name)[1].lower()
    if ext in ARCHIVE_EXTENSIONS:
        return True
    return False


def unpack_vcs_link(link, location):
    vcs_backend = _get_used_vcs_backend(link)
    vcs_backend.unpack(location)


def _get_used_vcs_backend(link):
    for backend in vcs.backends:
        if link.scheme in backend.schemes:
            vcs_backend = backend(link.url)
            return vcs_backend


def is_vcs_url(link):
    # type: (Link) -> bool
    return bool(_get_used_vcs_backend(link))


def is_file_url(link):
    # type: (Link) -> bool
    return link.url.lower().startswith('file:')


def is_dir_url(link):
    # type: (Link) -> bool
    """Return whether a file:// Link points to a directory.

    ``link`` must not have any other scheme but file://. Call is_file_url()
    first.

    """
    link_path = url_to_path(link.url_without_fragment)
    return os.path.isdir(link_path)


def _progress_indicator(iterable, *args, **kwargs):
    return iterable


def _download_url(
    resp,  # type: Response
    link,  # type: Link
    content_file,  # type: IO
    hashes,  # type: Hashes
    progress_bar  # type: str
):
    # type: (...) -> None
    try:
        total_length = int(resp.headers['content-length'])
    except (ValueError, KeyError, TypeError):
        total_length = 0

    cached_resp = getattr(resp, "from_cache", False)
    if logger.getEffectiveLevel() > logging.INFO:
        show_progress = False
    elif cached_resp:
        show_progress = False
    elif total_length > (40 * 1000):
        show_progress = True
    elif not total_length:
        show_progress = True
    else:
        show_progress = False

    show_url = link.show_url

    def resp_read(chunk_size):
        try:
            # Special case for urllib3.
            for chunk in resp.raw.stream(
                    chunk_size,
                    # We use decode_content=False here because we don't
                    # want urllib3 to mess with the raw bytes we get
                    # from the server. If we decompress inside of
                    # urllib3 then we cannot verify the checksum
                    # because the checksum will be of the compressed
                    # file. This breakage will only occur if the
                    # server adds a Content-Encoding header, which
                    # depends on how the server was configured:
                    # - Some servers will notice that the file isn't a
                    #   compressible file and will leave the file alone
                    #   and with an empty Content-Encoding
                    # - Some servers will notice that the file is
                    #   already compressed and will leave the file
                    #   alone and will add a Content-Encoding: gzip
                    #   header
                    # - Some servers won't notice anything at all and
                    #   will take a file that's already been compressed
                    #   and compress it again and set the
                    #   Content-Encoding: gzip header
                    #
                    # By setting this not to decode automatically we
                    # hope to eliminate problems with the second case.
                    decode_content=False):
                yield chunk
        except AttributeError:
            # Standard file-like object.
            while True:
                chunk = resp.raw.read(chunk_size)
                if not chunk:
                    break
                yield chunk

    def written_chunks(chunks):
        for chunk in chunks:
            content_file.write(chunk)
            yield chunk

    progress_indicator = _progress_indicator

    if link.netloc == PyPI.netloc:
        url = show_url
    else:
        url = link.url_without_fragment

    if show_progress:  # We don't show progress on cached responses
        progress_indicator = DownloadProgressProvider(progress_bar,
                                                      max=total_length)
        if total_length:
            logger.info("Downloading %s (%s)", url, format_size(total_length))
        else:
            logger.info("Downloading %s", url)
    elif cached_resp:
        logger.info("Using cached %s", url)
    else:
        logger.info("Downloading %s", url)

    logger.debug('Downloading from URL %s', link)

    downloaded_chunks = written_chunks(
        progress_indicator(
            resp_read(CONTENT_CHUNK_SIZE),
            CONTENT_CHUNK_SIZE
        )
    )
    if hashes:
        hashes.check_against_chunks(downloaded_chunks)
    else:
        consume(downloaded_chunks)


def _copy_file(filename, location, link):
    copy = True
    download_location = os.path.join(location, link.filename)
    if os.path.exists(download_location):
        response = ask_path_exists(
            'The file %s exists. (i)gnore, (w)ipe, (b)ackup, (a)abort' %
            display_path(download_location), ('i', 'w', 'b', 'a'))
        if response == 'i':
            copy = False
        elif response == 'w':
            logger.warning('Deleting %s', display_path(download_location))
            os.remove(download_location)
        elif response == 'b':
            dest_file = backup_dir(download_location)
            logger.warning(
                'Backing up %s to %s',
                display_path(download_location),
                display_path(dest_file),
            )
            shutil.move(download_location, dest_file)
        elif response == 'a':
            sys.exit(-1)
    if copy:
        shutil.copy(filename, download_location)
        logger.info('Saved %s', display_path(download_location))


def unpack_http_url(
    link,  # type: Link
    location,  # type: str
    download_dir=None,  # type: Optional[str]
    session=None,  # type: Optional[PipSession]
    hashes=None,  # type: Optional[Hashes]
    progress_bar="on"  # type: str
):
    # type: (...) -> None
    if session is None:
        raise TypeError(
            "unpack_http_url() missing 1 required keyword argument: 'session'"
        )

    with TempDirectory(kind="unpack") as temp_dir:
        # If a download dir is specified, is the file already downloaded there?
        already_downloaded_path = None
        if download_dir:
            already_downloaded_path = _check_download_dir(link,
                                                          download_dir,
                                                          hashes)

        if already_downloaded_path:
            from_path = already_downloaded_path
            content_type = mimetypes.guess_type(from_path)[0]
        else:
            # let's download to a tmp dir
            from_path, content_type = _download_http_url(link,
                                                         session,
                                                         temp_dir.path,
                                                         hashes,
                                                         progress_bar)

        # unpack the archive to the build dir location. even when only
        # downloading archives, they have to be unpacked to parse dependencies
        unpack_file(from_path, location, content_type, link)

        # a download dir is specified; let's copy the archive there
        if download_dir and not already_downloaded_path:
            _copy_file(from_path, download_dir, link)

        if not already_downloaded_path:
            os.unlink(from_path)


def unpack_file_url(
    link,  # type: Link
    location,  # type: str
    download_dir=None,  # type: Optional[str]
    hashes=None  # type: Optional[Hashes]
):
    # type: (...) -> None
    """Unpack link into location.

    If download_dir is provided and link points to a file, make a copy
    of the link file inside download_dir.
    """
    link_path = url_to_path(link.url_without_fragment)

    # If it's a url to a local directory
    if is_dir_url(link):
        if os.path.isdir(location):
            rmtree(location)
        shutil.copytree(link_path, location, symlinks=True)
        if download_dir:
            logger.info('Link is a directory, ignoring download_dir')
        return

    # If --require-hashes is off, `hashes` is either empty, the
    # link's embedded hash, or MissingHashes; it is required to
    # match. If --require-hashes is on, we are satisfied by any
    # hash in `hashes` matching: a URL-based or an option-based
    # one; no internet-sourced hash will be in `hashes`.
    if hashes:
        hashes.check_against_path(link_path)

    # If a download dir is specified, is the file already there and valid?
    already_downloaded_path = None
    if download_dir:
        already_downloaded_path = _check_download_dir(link,
                                                      download_dir,
                                                      hashes)

    if already_downloaded_path:
        from_path = already_downloaded_path
    else:
        from_path = link_path

    content_type = mimetypes.guess_type(from_path)[0]

    # unpack the archive to the build dir location. even when only downloading
    # archives, they have to be unpacked to parse dependencies
    unpack_file(from_path, location, content_type, link)

    # a download dir is specified and not already downloaded
    if download_dir and not already_downloaded_path:
        _copy_file(from_path, download_dir, link)


class PipXmlrpcTransport(xmlrpc_client.Transport):
    """Provide a `xmlrpclib.Transport` implementation via a `PipSession`
    object.
    """

    def __init__(self, index_url, session, use_datetime=False):
        xmlrpc_client.Transport.__init__(self, use_datetime)
        index_parts = urllib_parse.urlparse(index_url)
        self._scheme = index_parts.scheme
        self._session = session

    def request(self, host, handler, request_body, verbose=False):
        parts = (self._scheme, host, handler, None, None, None)
        url = urllib_parse.urlunparse(parts)
        try:
            headers = {'Content-Type': 'text/xml'}
            response = self._session.post(url, data=request_body,
                                          headers=headers, stream=True)
            response.raise_for_status()
            self.verbose = verbose
            return self.parse_response(response.raw)
        except requests.HTTPError as exc:
            logger.critical(
                "HTTP error %s while getting %s",
                exc.response.status_code, url,
            )
            raise


def unpack_url(
    link,  # type: Optional[Link]
    location,  # type: Optional[str]
    download_dir=None,  # type: Optional[str]
    only_download=False,  # type: bool
    session=None,  # type: Optional[PipSession]
    hashes=None,  # type: Optional[Hashes]
    progress_bar="on"  # type: str
):
    # type: (...) -> None
    """Unpack link.
       If link is a VCS link:
         if only_download, export into download_dir and ignore location
          else unpack into location
       for other types of link:
         - unpack into location
         - if download_dir, copy the file into download_dir
         - if only_download, mark location for deletion

    :param hashes: A Hashes object, one of whose embedded hashes must match,
        or HashMismatch will be raised. If the Hashes is empty, no matches are
        required, and unhashable types of requirements (like VCS ones, which
        would ordinarily raise HashUnsupported) are allowed.
    """
    # non-editable vcs urls
    if is_vcs_url(link):
        unpack_vcs_link(link, location)

    # file urls
    elif is_file_url(link):
        unpack_file_url(link, location, download_dir, hashes=hashes)

    # http urls
    else:
        if session is None:
            session = PipSession()

        unpack_http_url(
            link,
            location,
            download_dir,
            session,
            hashes=hashes,
            progress_bar=progress_bar
        )
    if only_download:
        write_delete_marker_file(location)


def _download_http_url(
    link,  # type: Link
    session,  # type: PipSession
    temp_dir,  # type: str
    hashes,  # type: Hashes
    progress_bar  # type: str
):
    # type: (...) -> Tuple[str, str]
    """Download link url into temp_dir using provided session"""
    target_url = link.url.split('#', 1)[0]
    try:
        resp = session.get(
            target_url,
            # We use Accept-Encoding: identity here because requests
            # defaults to accepting compressed responses. This breaks in
            # a variety of ways depending on how the server is configured.
            # - Some servers will notice that the file isn't a compressible
            #   file and will leave the file alone and with an empty
            #   Content-Encoding
            # - Some servers will notice that the file is already
            #   compressed and will leave the file alone and will add a
            #   Content-Encoding: gzip header
            # - Some servers won't notice anything at all and will take
            #   a file that's already been compressed and compress it again
            #   and set the Content-Encoding: gzip header
            # By setting this to request only the identity encoding We're
            # hoping to eliminate the third case. Hopefully there does not
            # exist a server which when given a file will notice it is
            # already compressed and that you're not asking for a
            # compressed file and will then decompress it before sending
            # because if that's the case I don't think it'll ever be
            # possible to make this work.
            headers={"Accept-Encoding": "identity"},
            stream=True,
        )
        resp.raise_for_status()
    except requests.HTTPError as exc:
        logger.critical(
            "HTTP error %s while getting %s", exc.response.status_code, link,
        )
        raise

    content_type = resp.headers.get('content-type', '')
    filename = link.filename  # fallback
    # Have a look at the Content-Disposition header for a better guess
    content_disposition = resp.headers.get('content-disposition')
    if content_disposition:
        type, params = cgi.parse_header(content_disposition)
        # We use ``or`` here because we don't want to use an "empty" value
        # from the filename param.
        filename = params.get('filename') or filename
    ext = splitext(filename)[1]
    if not ext:
        ext = mimetypes.guess_extension(content_type)
        if ext:
            filename += ext
    if not ext and link.url != resp.url:
        ext = os.path.splitext(resp.url)[1]
        if ext:
            filename += ext
    file_path = os.path.join(temp_dir, filename)
    with open(file_path, 'wb') as content_file:
        _download_url(resp, link, content_file, hashes, progress_bar)
    return file_path, content_type


def _check_download_dir(link, download_dir, hashes):
    # type: (Link, str, Hashes) -> Optional[str]
    """ Check download_dir for previously downloaded file with correct hash
        If a correct file is found return its path else None
    """
    download_path = os.path.join(download_dir, link.filename)
    if os.path.exists(download_path):
        # If already downloaded, does its hash match?
        logger.info('File was already downloaded %s', download_path)
        if hashes:
            try:
                hashes.check_against_path(download_path)
            except HashMismatch:
                logger.warning(
                    'Previously-downloaded file %s has bad hash. '
                    'Re-downloading.',
                    download_path
                )
                os.unlink(download_path)
                return None
        return download_path
    return None
