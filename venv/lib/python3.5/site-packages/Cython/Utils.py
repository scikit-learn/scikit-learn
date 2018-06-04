#
#   Cython -- Things that don't belong
#            anywhere else in particular
#

from __future__ import absolute_import

try:
    from __builtin__ import basestring
except ImportError:
    basestring = str

import os
import sys
import re
import io
import codecs
import shutil
from contextlib import contextmanager

modification_time = os.path.getmtime


def cached_function(f):
    cache = {}
    uncomputed = object()
    def wrapper(*args):
        res = cache.get(args, uncomputed)
        if res is uncomputed:
            res = cache[args] = f(*args)
        return res
    wrapper.uncached = f
    return wrapper

def cached_method(f):
    cache_name = '__%s_cache' % f.__name__
    def wrapper(self, *args):
        cache = getattr(self, cache_name, None)
        if cache is None:
            cache = {}
            setattr(self, cache_name, cache)
        if args in cache:
            return cache[args]
        res = cache[args] = f(self, *args)
        return res
    return wrapper

def replace_suffix(path, newsuf):
    base, _ = os.path.splitext(path)
    return base + newsuf


def open_new_file(path):
    if os.path.exists(path):
        # Make sure to create a new file here so we can
        # safely hard link the output files.
        os.unlink(path)

    # we use the ISO-8859-1 encoding here because we only write pure
    # ASCII strings or (e.g. for file names) byte encoded strings as
    # Unicode, so we need a direct mapping from the first 256 Unicode
    # characters to a byte sequence, which ISO-8859-1 provides

    # note: can't use io.open() in Py2 as we may be writing str objects
    return codecs.open(path, "w", encoding="ISO-8859-1")


def castrate_file(path, st):
    #  Remove junk contents from an output file after a
    #  failed compilation.
    #  Also sets access and modification times back to
    #  those specified by st (a stat struct).
    try:
        f = open_new_file(path)
    except EnvironmentError:
        pass
    else:
        f.write(
            "#error Do not use this file, it is the result of a failed Cython compilation.\n")
        f.close()
        if st:
            os.utime(path, (st.st_atime, st.st_mtime-1))

def file_newer_than(path, time):
    ftime = modification_time(path)
    return ftime > time


def safe_makedirs(path):
    try:
        os.makedirs(path)
    except OSError:
        if not os.path.isdir(path):
            raise


def copy_file_to_dir_if_newer(sourcefile, destdir):
    """
    Copy file sourcefile to directory destdir (creating it if needed),
    preserving metadata. If the destination file exists and is not
    older than the source file, the copying is skipped.
    """
    destfile = os.path.join(destdir, os.path.basename(sourcefile))
    try:
        desttime = modification_time(destfile)
    except OSError:
        # New file does not exist, destdir may or may not exist
        safe_makedirs(destdir)
    else:
        # New file already exists
        if not file_newer_than(sourcefile, desttime):
            return
    shutil.copy2(sourcefile, destfile)


@cached_function
def search_include_directories(dirs, qualified_name, suffix, pos,
                               include=False, sys_path=False):
    # Search the list of include directories for the given
    # file name. If a source file position is given, first
    # searches the directory containing that file. Returns
    # None if not found, but does not report an error.
    # The 'include' option will disable package dereferencing.
    # If 'sys_path' is True, also search sys.path.
    if sys_path:
        dirs = dirs + tuple(sys.path)
    if pos:
        file_desc = pos[0]
        from Cython.Compiler.Scanning import FileSourceDescriptor
        if not isinstance(file_desc, FileSourceDescriptor):
            raise RuntimeError("Only file sources for code supported")
        if include:
            dirs = (os.path.dirname(file_desc.filename),) + dirs
        else:
            dirs = (find_root_package_dir(file_desc.filename),) + dirs

    dotted_filename = qualified_name
    if suffix:
        dotted_filename += suffix
    if not include:
        names = qualified_name.split('.')
        package_names = tuple(names[:-1])
        module_name = names[-1]
        module_filename = module_name + suffix
        package_filename = "__init__" + suffix

    for dir in dirs:
        path = os.path.join(dir, dotted_filename)
        if path_exists(path):
            return path
        if not include:
            package_dir = check_package_dir(dir, package_names)
            if package_dir is not None:
                path = os.path.join(package_dir, module_filename)
                if path_exists(path):
                    return path
                path = os.path.join(dir, package_dir, module_name,
                                    package_filename)
                if path_exists(path):
                    return path
    return None


@cached_function
def find_root_package_dir(file_path):
    dir = os.path.dirname(file_path)
    if file_path == dir:
        return dir
    elif is_package_dir(dir):
        return find_root_package_dir(dir)
    else:
        return dir

@cached_function
def check_package_dir(dir, package_names):
    for dirname in package_names:
        dir = os.path.join(dir, dirname)
        if not is_package_dir(dir):
            return None
    return dir

@cached_function
def is_package_dir(dir_path):
    for filename in ("__init__.py",
                     "__init__.pyc",
                     "__init__.pyx",
                     "__init__.pxd"):
        path = os.path.join(dir_path, filename)
        if path_exists(path):
            return 1

@cached_function
def path_exists(path):
    # try on the filesystem first
    if os.path.exists(path):
        return True
    # figure out if a PEP 302 loader is around
    try:
        loader = __loader__
        # XXX the code below assumes a 'zipimport.zipimporter' instance
        # XXX should be easy to generalize, but too lazy right now to write it
        archive_path = getattr(loader, 'archive', None)
        if archive_path:
            normpath = os.path.normpath(path)
            if normpath.startswith(archive_path):
                arcname = normpath[len(archive_path)+1:]
                try:
                    loader.get_data(arcname)
                    return True
                except IOError:
                    return False
    except NameError:
        pass
    return False

# file name encodings

def decode_filename(filename):
    if isinstance(filename, bytes):
        try:
            filename_encoding = sys.getfilesystemencoding()
            if filename_encoding is None:
                filename_encoding = sys.getdefaultencoding()
            filename = filename.decode(filename_encoding)
        except UnicodeDecodeError:
            pass
    return filename

# support for source file encoding detection

_match_file_encoding = re.compile(u"coding[:=]\s*([-\w.]+)").search


def detect_file_encoding(source_filename):
    f = open_source_file(source_filename, encoding="UTF-8", error_handling='ignore')
    try:
        return detect_opened_file_encoding(f)
    finally:
        f.close()


def detect_opened_file_encoding(f):
    # PEPs 263 and 3120
    # Most of the time the first two lines fall in the first 250 chars,
    # and this bulk read/split is much faster.
    lines = f.read(250).split(u"\n")
    if len(lines) > 1:
        m = _match_file_encoding(lines[0])
        if m:
            return m.group(1)
        elif len(lines) > 2:
            m = _match_file_encoding(lines[1])
            if m:
                return m.group(1)
            else:
                return "UTF-8"
    # Fallback to one-char-at-a-time detection.
    f.seek(0)
    chars = []
    for i in range(2):
        c = f.read(1)
        while c and c != u'\n':
            chars.append(c)
            c = f.read(1)
        encoding = _match_file_encoding(u''.join(chars))
        if encoding:
            return encoding.group(1)
    return "UTF-8"


def skip_bom(f):
    """
    Read past a BOM at the beginning of a source file.
    This could be added to the scanner, but it's *substantially* easier
    to keep it at this level.
    """
    if f.read(1) != u'\uFEFF':
        f.seek(0)


def open_source_file(source_filename, mode="r",
                     encoding=None, error_handling=None):
    if encoding is None:
        # Most of the time the coding is unspecified, so be optimistic that
        # it's UTF-8.
        f = open_source_file(source_filename, encoding="UTF-8", mode=mode, error_handling='ignore')
        encoding = detect_opened_file_encoding(f)
        if encoding == "UTF-8" and error_handling == 'ignore':
            f.seek(0)
            skip_bom(f)
            return f
        else:
            f.close()

    if not os.path.exists(source_filename):
        try:
            loader = __loader__
            if source_filename.startswith(loader.archive):
                return open_source_from_loader(
                    loader, source_filename,
                    encoding, error_handling)
        except (NameError, AttributeError):
            pass

    stream = io.open(source_filename, mode=mode,
                     encoding=encoding, errors=error_handling)
    skip_bom(stream)
    return stream


def open_source_from_loader(loader,
                            source_filename,
                            encoding=None, error_handling=None):
    nrmpath = os.path.normpath(source_filename)
    arcname = nrmpath[len(loader.archive)+1:]
    data = loader.get_data(arcname)
    return io.TextIOWrapper(io.BytesIO(data),
                            encoding=encoding,
                            errors=error_handling)


def str_to_number(value):
    # note: this expects a string as input that was accepted by the
    # parser already, with an optional "-" sign in front
    is_neg = False
    if value[:1] == '-':
        is_neg = True
        value = value[1:]
    if len(value) < 2:
        value = int(value, 0)
    elif value[0] == '0':
        literal_type = value[1]  # 0'o' - 0'b' - 0'x'
        if literal_type in 'xX':
            # hex notation ('0x1AF')
            value = int(value[2:], 16)
        elif literal_type in 'oO':
            # Py3 octal notation ('0o136')
            value = int(value[2:], 8)
        elif literal_type in 'bB':
            # Py3 binary notation ('0b101')
            value = int(value[2:], 2)
        else:
            # Py2 octal notation ('0136')
            value = int(value, 8)
    else:
        value = int(value, 0)
    return -value if is_neg else value


def long_literal(value):
    if isinstance(value, basestring):
        value = str_to_number(value)
    return not -2**31 <= value < 2**31


@cached_function
def get_cython_cache_dir():
    """get the cython cache dir

    Priority:

    1. CYTHON_CACHE_DIR
    2. (OS X): ~/Library/Caches/Cython
       (posix not OS X): XDG_CACHE_HOME/cython if XDG_CACHE_HOME defined
    3. ~/.cython

    """
    if 'CYTHON_CACHE_DIR' in os.environ:
        return os.environ['CYTHON_CACHE_DIR']

    parent = None
    if os.name == 'posix':
        if sys.platform == 'darwin':
            parent = os.path.expanduser('~/Library/Caches')
        else:
            # this could fallback on ~/.cache
            parent = os.environ.get('XDG_CACHE_HOME')

    if parent and os.path.isdir(parent):
        return os.path.join(parent, 'cython')

    # last fallback: ~/.cython
    return os.path.expanduser(os.path.join('~', '.cython'))


@contextmanager
def captured_fd(stream=2, encoding=None):
    pipe_in = t = None
    orig_stream = os.dup(stream)  # keep copy of original stream
    try:
        pipe_in, pipe_out = os.pipe()
        os.dup2(pipe_out, stream)  # replace stream by copy of pipe
        try:
            os.close(pipe_out)  # close original pipe-out stream
            data = []

            def copy():
                try:
                    while True:
                        d = os.read(pipe_in, 1000)
                        if d:
                            data.append(d)
                        else:
                            break
                finally:
                    os.close(pipe_in)

            def get_output():
                output = b''.join(data)
                if encoding:
                    output = output.decode(encoding)
                return output

            from threading import Thread
            t = Thread(target=copy)
            t.daemon = True  # just in case
            t.start()
            yield get_output
        finally:
            os.dup2(orig_stream, stream)  # restore original stream
            if t is not None:
                t.join()
    finally:
        os.close(orig_stream)


def print_bytes(s, end=b'\n', file=sys.stdout, flush=True):
    file.flush()
    try:
        out = file.buffer  # Py3
    except AttributeError:
        out = file         # Py2
    out.write(s)
    if end:
        out.write(end)
    if flush:
        out.flush()

class LazyStr:
    def __init__(self, callback):
        self.callback = callback
    def __str__(self):
        return self.callback()
    def __repr__(self):
        return self.callback()
    def __add__(self, right):
        return self.callback() + right
    def __radd__(self, left):
        return left + self.callback()


class OrderedSet(object):
  def __init__(self, elements=()):
    self._list = []
    self._set = set()
    self.update(elements)
  def __iter__(self):
    return iter(self._list)
  def update(self, elements):
    for e in elements:
      self.add(e)
  def add(self, e):
    if e not in self._set:
      self._list.append(e)
      self._set.add(e)


# Class decorator that adds a metaclass and recreates the class with it.
# Copied from 'six'.
def add_metaclass(metaclass):
    """Class decorator for creating a class with a metaclass."""
    def wrapper(cls):
        orig_vars = cls.__dict__.copy()
        slots = orig_vars.get('__slots__')
        if slots is not None:
            if isinstance(slots, str):
                slots = [slots]
            for slots_var in slots:
                orig_vars.pop(slots_var)
        orig_vars.pop('__dict__', None)
        orig_vars.pop('__weakref__', None)
        return metaclass(cls.__name__, cls.__bases__, orig_vars)
    return wrapper
