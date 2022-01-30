"""
    sphinx.util
    ~~~~~~~~~~~

    Utility functions for Sphinx.

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import functools
import hashlib
import os
import posixpath
import re
import sys
import tempfile
import traceback
import unicodedata
import warnings
from datetime import datetime
from importlib import import_module
from os import path
from time import mktime, strptime
from typing import (IO, TYPE_CHECKING, Any, Callable, Dict, Iterable, Iterator, List, Optional,
                    Pattern, Set, Tuple, Type)
from urllib.parse import parse_qsl, quote_plus, urlencode, urlsplit, urlunsplit

from sphinx.deprecation import RemovedInSphinx50Warning
from sphinx.errors import ExtensionError, FiletypeNotFoundError, SphinxParallelError
from sphinx.locale import __
from sphinx.util import logging
from sphinx.util.console import bold, colorize, strip_colors, term_width_line  # type: ignore
from sphinx.util.matching import patfilter  # noqa
from sphinx.util.nodes import (caption_ref_re, explicit_title_re,  # noqa
                               nested_parse_with_titles, split_explicit_title)
# import other utilities; partly for backwards compatibility, so don't
# prune unused ones indiscriminately
from sphinx.util.osutil import (SEP, copyfile, copytimes, ensuredir, make_filename,  # noqa
                                movefile, mtimes_of_files, os_path, relative_uri)
from sphinx.util.typing import PathMatcher

if TYPE_CHECKING:
    from sphinx.application import Sphinx


logger = logging.getLogger(__name__)

# Generally useful regular expressions.
ws_re: Pattern = re.compile(r'\s+')
url_re: Pattern = re.compile(r'(?P<schema>.+)://.*')


# High-level utility functions.

def docname_join(basedocname: str, docname: str) -> str:
    return posixpath.normpath(
        posixpath.join('/' + basedocname, '..', docname))[1:]


def path_stabilize(filepath: str) -> str:
    "Normalize path separator and unicode string"
    newpath = filepath.replace(os.path.sep, SEP)
    return unicodedata.normalize('NFC', newpath)


def get_matching_files(dirname: str,
                       exclude_matchers: Tuple[PathMatcher, ...] = ()) -> Iterable[str]:  # NOQA
    """Get all file names in a directory, recursively.

    Exclude files and dirs matching some matcher in *exclude_matchers*.
    """
    # dirname is a normalized absolute path.
    dirname = path.normpath(path.abspath(dirname))

    for root, dirs, files in os.walk(dirname, followlinks=True):
        relativeroot = path.relpath(root, dirname)
        if relativeroot == ".":
            relativeroot = ""  # suppress dirname for files on the target dir

        qdirs = enumerate(path_stabilize(path.join(relativeroot, dn))
                          for dn in dirs)  # type: Iterable[Tuple[int, str]]
        qfiles = enumerate(path_stabilize(path.join(relativeroot, fn))
                           for fn in files)  # type: Iterable[Tuple[int, str]]
        for matcher in exclude_matchers:
            qdirs = [entry for entry in qdirs if not matcher(entry[1])]
            qfiles = [entry for entry in qfiles if not matcher(entry[1])]

        dirs[:] = sorted(dirs[i] for (i, _) in qdirs)

        for _i, filename in sorted(qfiles):
            yield filename


def get_filetype(source_suffix: Dict[str, str], filename: str) -> str:
    for suffix, filetype in source_suffix.items():
        if filename.endswith(suffix):
            # If default filetype (None), considered as restructuredtext.
            return filetype or 'restructuredtext'
    else:
        raise FiletypeNotFoundError


class FilenameUniqDict(dict):
    """
    A dictionary that automatically generates unique names for its keys,
    interpreted as filenames, and keeps track of a set of docnames they
    appear in.  Used for images and downloadable files in the environment.
    """
    def __init__(self) -> None:
        self._existing: Set[str] = set()

    def add_file(self, docname: str, newfile: str) -> str:
        if newfile in self:
            self[newfile][0].add(docname)
            return self[newfile][1]
        uniquename = path.basename(newfile)
        base, ext = path.splitext(uniquename)
        i = 0
        while uniquename in self._existing:
            i += 1
            uniquename = '%s%s%s' % (base, i, ext)
        self[newfile] = ({docname}, uniquename)
        self._existing.add(uniquename)
        return uniquename

    def purge_doc(self, docname: str) -> None:
        for filename, (docs, unique) in list(self.items()):
            docs.discard(docname)
            if not docs:
                del self[filename]
                self._existing.discard(unique)

    def merge_other(self, docnames: Set[str], other: Dict[str, Tuple[Set[str], Any]]) -> None:
        for filename, (docs, _unique) in other.items():
            for doc in docs & set(docnames):
                self.add_file(doc, filename)

    def __getstate__(self) -> Set[str]:
        return self._existing

    def __setstate__(self, state: Set[str]) -> None:
        self._existing = state


def md5(data=b'', **kwargs):
    """Wrapper around hashlib.md5

    Attempt call with 'usedforsecurity=False' if we get a ValueError, which happens when
    OpenSSL FIPS mode is enabled:
    ValueError: error:060800A3:digital envelope routines:EVP_DigestInit_ex:disabled for fips

    See: https://github.com/sphinx-doc/sphinx/issues/7611
    """

    try:
        return hashlib.md5(data, **kwargs)
    except ValueError:
        return hashlib.md5(data, **kwargs, usedforsecurity=False)  # type: ignore


def sha1(data=b'', **kwargs):
    """Wrapper around hashlib.sha1

    Attempt call with 'usedforsecurity=False' if we get a ValueError

    See: https://github.com/sphinx-doc/sphinx/issues/7611
    """

    try:
        return hashlib.sha1(data, **kwargs)
    except ValueError:
        return hashlib.sha1(data, **kwargs, usedforsecurity=False)  # type: ignore


class DownloadFiles(dict):
    """A special dictionary for download files.

    .. important:: This class would be refactored in nearly future.
                   Hence don't hack this directly.
    """

    def add_file(self, docname: str, filename: str) -> str:
        if filename not in self:
            digest = md5(filename.encode()).hexdigest()
            dest = '%s/%s' % (digest, os.path.basename(filename))
            self[filename] = (set(), dest)

        self[filename][0].add(docname)
        return self[filename][1]

    def purge_doc(self, docname: str) -> None:
        for filename, (docs, _dest) in list(self.items()):
            docs.discard(docname)
            if not docs:
                del self[filename]

    def merge_other(self, docnames: Set[str], other: Dict[str, Tuple[Set[str], Any]]) -> None:
        for filename, (docs, _dest) in other.items():
            for docname in docs & set(docnames):
                self.add_file(docname, filename)


_DEBUG_HEADER = '''\
# Sphinx version: %s
# Python version: %s (%s)
# Docutils version: %s %s
# Jinja2 version: %s
# Last messages:
%s
# Loaded extensions:
'''


def save_traceback(app: "Sphinx") -> str:
    """Save the current exception's traceback in a temporary file."""
    import platform

    import docutils
    import jinja2

    import sphinx
    exc = sys.exc_info()[1]
    if isinstance(exc, SphinxParallelError):
        exc_format = '(Error in parallel process)\n' + exc.traceback
    else:
        exc_format = traceback.format_exc()
    fd, path = tempfile.mkstemp('.log', 'sphinx-err-')
    last_msgs = ''
    if app is not None:
        last_msgs = '\n'.join(
            '#   %s' % strip_colors(s).strip()
            for s in app.messagelog)
    os.write(fd, (_DEBUG_HEADER %
                  (sphinx.__display_version__,
                   platform.python_version(),
                   platform.python_implementation(),
                   docutils.__version__, docutils.__version_details__,
                   jinja2.__version__,
                   last_msgs)).encode())
    if app is not None:
        for ext in app.extensions.values():
            modfile = getattr(ext.module, '__file__', 'unknown')
            if ext.version != 'builtin':
                os.write(fd, ('#   %s (%s) from %s\n' %
                              (ext.name, ext.version, modfile)).encode())
    os.write(fd, exc_format.encode())
    os.close(fd)
    return path


def get_full_modname(modname: str, attribute: str) -> str:
    if modname is None:
        # Prevents a TypeError: if the last getattr() call will return None
        # then it's better to return it directly
        return None
    module = import_module(modname)

    # Allow an attribute to have multiple parts and incidentally allow
    # repeated .s in the attribute.
    value = module
    for attr in attribute.split('.'):
        if attr:
            value = getattr(value, attr)

    return getattr(value, '__module__', None)


# a regex to recognize coding cookies
_coding_re = re.compile(r'coding[:=]\s*([-\w.]+)')


class UnicodeDecodeErrorHandler:
    """Custom error handler for open() that warns and replaces."""

    def __init__(self, docname: str) -> None:
        self.docname = docname

    def __call__(self, error: UnicodeDecodeError) -> Tuple[str, int]:
        linestart = error.object.rfind(b'\n', 0, error.start)
        lineend = error.object.find(b'\n', error.start)
        if lineend == -1:
            lineend = len(error.object)
        lineno = error.object.count(b'\n', 0, error.start) + 1
        logger.warning(__('undecodable source characters, replacing with "?": %r'),
                       (error.object[linestart + 1:error.start] + b'>>>' +
                        error.object[error.start:error.end] + b'<<<' +
                        error.object[error.end:lineend]),
                       location=(self.docname, lineno))
        return ('?', error.end)


# Low-level utility functions and classes.

class Tee:
    """
    File-like object writing to two streams.
    """
    def __init__(self, stream1: IO, stream2: IO) -> None:
        self.stream1 = stream1
        self.stream2 = stream2

    def write(self, text: str) -> None:
        self.stream1.write(text)
        self.stream2.write(text)

    def flush(self) -> None:
        if hasattr(self.stream1, 'flush'):
            self.stream1.flush()
        if hasattr(self.stream2, 'flush'):
            self.stream2.flush()


def parselinenos(spec: str, total: int) -> List[int]:
    """Parse a line number spec (such as "1,2,4-6") and return a list of
    wanted line numbers.
    """
    items = list()
    parts = spec.split(',')
    for part in parts:
        try:
            begend = part.strip().split('-')
            if ['', ''] == begend:
                raise ValueError
            elif len(begend) == 1:
                items.append(int(begend[0]) - 1)
            elif len(begend) == 2:
                start = int(begend[0] or 1)  # left half open (cf. -10)
                end = int(begend[1] or max(start, total))  # right half open (cf. 10-)
                if start > end:  # invalid range (cf. 10-1)
                    raise ValueError
                items.extend(range(start - 1, end))
            else:
                raise ValueError
        except Exception as exc:
            raise ValueError('invalid line number spec: %r' % spec) from exc

    return items


def force_decode(string: str, encoding: str) -> str:
    """Forcibly get a unicode string out of a bytestring."""
    warnings.warn('force_decode() is deprecated.',
                  RemovedInSphinx50Warning, stacklevel=2)
    if isinstance(string, bytes):
        try:
            if encoding:
                string = string.decode(encoding)
            else:
                # try decoding with utf-8, should only work for real UTF-8
                string = string.decode()
        except UnicodeError:
            # last resort -- can't fail
            string = string.decode('latin1')
    return string


def rpartition(s: str, t: str) -> Tuple[str, str]:
    """Similar to str.rpartition from 2.5, but doesn't return the separator."""
    warnings.warn('rpartition() is now deprecated.', RemovedInSphinx50Warning, stacklevel=2)
    i = s.rfind(t)
    if i != -1:
        return s[:i], s[i + len(t):]
    return '', s


def split_into(n: int, type: str, value: str) -> List[str]:
    """Split an index entry into a given number of parts at semicolons."""
    parts = [x.strip() for x in value.split(';', n - 1)]
    if sum(1 for part in parts if part) < n:
        raise ValueError('invalid %s index entry %r' % (type, value))
    return parts


def split_index_msg(type: str, value: str) -> List[str]:
    # new entry types must be listed in directives/other.py!
    if type == 'single':
        try:
            result = split_into(2, 'single', value)
        except ValueError:
            result = split_into(1, 'single', value)
    elif type == 'pair':
        result = split_into(2, 'pair', value)
    elif type == 'triple':
        result = split_into(3, 'triple', value)
    elif type == 'see':
        result = split_into(2, 'see', value)
    elif type == 'seealso':
        result = split_into(2, 'see', value)
    else:
        raise ValueError('invalid %s index entry %r' % (type, value))

    return result


def format_exception_cut_frames(x: int = 1) -> str:
    """Format an exception with traceback, but only the last x frames."""
    typ, val, tb = sys.exc_info()
    # res = ['Traceback (most recent call last):\n']
    res: List[str] = []
    tbres = traceback.format_tb(tb)
    res += tbres[-x:]
    res += traceback.format_exception_only(typ, val)
    return ''.join(res)


def import_object(objname: str, source: str = None) -> Any:
    """Import python object by qualname."""
    try:
        objpath = objname.split('.')
        modname = objpath.pop(0)
        obj = import_module(modname)
        for name in objpath:
            modname += '.' + name
            try:
                obj = getattr(obj, name)
            except AttributeError:
                obj = import_module(modname)

        return obj
    except (AttributeError, ImportError) as exc:
        if source:
            raise ExtensionError('Could not import %s (needed for %s)' %
                                 (objname, source), exc) from exc
        else:
            raise ExtensionError('Could not import %s' % objname, exc) from exc


def split_full_qualified_name(name: str) -> Tuple[Optional[str], str]:
    """Split full qualified name to a pair of modname and qualname.

    A qualname is an abbreviation for "Qualified name" introduced at PEP-3155
    (https://www.python.org/dev/peps/pep-3155/).  It is a dotted path name
    from the module top-level.

    A "full" qualified name means a string containing both module name and
    qualified name.

    .. note:: This function actually imports the module to check its existence.
              Therefore you need to mock 3rd party modules if needed before
              calling this function.
    """
    parts = name.split('.')
    for i, _part in enumerate(parts, 1):
        try:
            modname = ".".join(parts[:i])
            import_module(modname)
        except ImportError:
            if parts[:i - 1]:
                return ".".join(parts[:i - 1]), ".".join(parts[i - 1:])
            else:
                return None, ".".join(parts)
        except IndexError:
            pass

    return name, ""


def encode_uri(uri: str) -> str:
    split = list(urlsplit(uri))
    split[1] = split[1].encode('idna').decode('ascii')
    split[2] = quote_plus(split[2].encode(), '/')
    query = list((q, v.encode()) for (q, v) in parse_qsl(split[3]))
    split[3] = urlencode(query)
    return urlunsplit(split)


def isurl(url: str) -> bool:
    """Check *url* is URL or not."""
    if url and '://' in url:
        return True
    else:
        return False


def display_chunk(chunk: Any) -> str:
    if isinstance(chunk, (list, tuple)):
        if len(chunk) == 1:
            return str(chunk[0])
        return '%s .. %s' % (chunk[0], chunk[-1])
    return str(chunk)


def old_status_iterator(iterable: Iterable, summary: str, color: str = "darkgreen",
                        stringify_func: Callable[[Any], str] = display_chunk) -> Iterator:
    l = 0
    for item in iterable:
        if l == 0:
            logger.info(bold(summary), nonl=True)
            l = 1
        logger.info(stringify_func(item), color=color, nonl=True)
        logger.info(" ", nonl=True)
        yield item
    if l == 1:
        logger.info('')


# new version with progress info
def status_iterator(iterable: Iterable, summary: str, color: str = "darkgreen",
                    length: int = 0, verbosity: int = 0,
                    stringify_func: Callable[[Any], str] = display_chunk) -> Iterable:
    if length == 0:
        yield from old_status_iterator(iterable, summary, color, stringify_func)
        return
    l = 0
    summary = bold(summary)
    for item in iterable:
        l += 1
        s = '%s[%3d%%] %s' % (summary, 100 * l / length, colorize(color, stringify_func(item)))
        if verbosity:
            s += '\n'
        else:
            s = term_width_line(s)
        logger.info(s, nonl=True)
        yield item
    if l > 0:
        logger.info('')


class SkipProgressMessage(Exception):
    pass


class progress_message:
    def __init__(self, message: str) -> None:
        self.message = message

    def __enter__(self) -> None:
        logger.info(bold(self.message + '... '), nonl=True)

    def __exit__(self, exc_type: Type[Exception], exc_value: Exception, traceback: Any) -> bool:  # NOQA
        if isinstance(exc_value, SkipProgressMessage):
            logger.info(__('skipped'))
            if exc_value.args:
                logger.info(*exc_value.args)
            return True
        elif exc_type:
            logger.info(__('failed'))
        else:
            logger.info(__('done'))

        return False

    def __call__(self, f: Callable) -> Callable:
        @functools.wraps(f)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            with self:
                return f(*args, **kwargs)

        return wrapper


def epoch_to_rfc1123(epoch: float) -> str:
    """Convert datetime format epoch to RFC1123."""
    from babel.dates import format_datetime

    dt = datetime.fromtimestamp(epoch)
    fmt = 'EEE, dd LLL yyyy hh:mm:ss'
    return format_datetime(dt, fmt, locale='en') + ' GMT'


def rfc1123_to_epoch(rfc1123: str) -> float:
    return mktime(strptime(rfc1123, '%a, %d %b %Y %H:%M:%S %Z'))


def xmlname_checker() -> Pattern:
    # https://www.w3.org/TR/REC-xml/#NT-Name
    name_start_chars = [
        ':', ['A', 'Z'], '_', ['a', 'z'], ['\u00C0', '\u00D6'],
        ['\u00D8', '\u00F6'], ['\u00F8', '\u02FF'], ['\u0370', '\u037D'],
        ['\u037F', '\u1FFF'], ['\u200C', '\u200D'], ['\u2070', '\u218F'],
        ['\u2C00', '\u2FEF'], ['\u3001', '\uD7FF'], ['\uF900', '\uFDCF'],
        ['\uFDF0', '\uFFFD'], ['\U00010000', '\U000EFFFF']]

    name_chars = [
        "\\-", "\\.", ['0', '9'], '\u00B7', ['\u0300', '\u036F'],
        ['\u203F', '\u2040']
    ]

    def convert(entries: Any, splitter: str = '|') -> str:
        results = []
        for entry in entries:
            if isinstance(entry, list):
                results.append('[%s]' % convert(entry, '-'))
            else:
                results.append(entry)
        return splitter.join(results)

    start_chars_regex = convert(name_start_chars)
    name_chars_regex = convert(name_chars)
    return re.compile('(%s)(%s|%s)*' % (
        start_chars_regex, start_chars_regex, name_chars_regex))
