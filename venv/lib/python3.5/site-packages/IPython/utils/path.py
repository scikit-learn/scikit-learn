# encoding: utf-8
"""
Utilities for path handling.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import os
import sys
import errno
import shutil
import random
import glob
from warnings import warn

from IPython.utils.process import system
from IPython.utils import py3compat
from IPython.utils.decorators import undoc

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

fs_encoding = sys.getfilesystemencoding()

def _writable_dir(path):
    """Whether `path` is a directory, to which the user has write access."""
    return os.path.isdir(path) and os.access(path, os.W_OK)

if sys.platform == 'win32':
    def _get_long_path_name(path):
        """Get a long path name (expand ~) on Windows using ctypes.

        Examples
        --------

        >>> get_long_path_name('c:\\docume~1')
        'c:\\\\Documents and Settings'

        """
        try:
            import ctypes
        except ImportError:
            raise ImportError('you need to have ctypes installed for this to work')
        _GetLongPathName = ctypes.windll.kernel32.GetLongPathNameW
        _GetLongPathName.argtypes = [ctypes.c_wchar_p, ctypes.c_wchar_p,
            ctypes.c_uint ]

        buf = ctypes.create_unicode_buffer(260)
        rv = _GetLongPathName(path, buf, 260)
        if rv == 0 or rv > 260:
            return path
        else:
            return buf.value
else:
    def _get_long_path_name(path):
        """Dummy no-op."""
        return path



def get_long_path_name(path):
    """Expand a path into its long form.

    On Windows this expands any ~ in the paths. On other platforms, it is
    a null operation.
    """
    return _get_long_path_name(path)


def unquote_filename(name, win32=(sys.platform=='win32')):
    """ On Windows, remove leading and trailing quotes from filenames.

    This function has been deprecated and should not be used any more:
    unquoting is now taken care of by :func:`IPython.utils.process.arg_split`.
    """
    warn("'unquote_filename' is deprecated since IPython 5.0 and should not "
         "be used anymore", DeprecationWarning, stacklevel=2)
    if win32:
        if name.startswith(("'", '"')) and name.endswith(("'", '"')):
            name = name[1:-1]
    return name


def compress_user(path):
    """Reverse of :func:`os.path.expanduser`
    """
    home = os.path.expanduser('~')
    if path.startswith(home):
        path =  "~" + path[len(home):]
    return path

def get_py_filename(name, force_win32=None):
    """Return a valid python filename in the current directory.

    If the given name is not a file, it adds '.py' and searches again.
    Raises IOError with an informative message if the file isn't found.
    """

    name = os.path.expanduser(name)
    if force_win32 is not None:
        warn("The 'force_win32' argument to 'get_py_filename' is deprecated "
             "since IPython 5.0 and should not be used anymore",
            DeprecationWarning, stacklevel=2)
    if not os.path.isfile(name) and not name.endswith('.py'):
        name += '.py'
    if os.path.isfile(name):
        return name
    else:
        raise IOError('File `%r` not found.' % name)


def filefind(filename, path_dirs=None):
    """Find a file by looking through a sequence of paths.

    This iterates through a sequence of paths looking for a file and returns
    the full, absolute path of the first occurrence of the file.  If no set of
    path dirs is given, the filename is tested as is, after running through
    :func:`expandvars` and :func:`expanduser`.  Thus a simple call::

        filefind('myfile.txt')

    will find the file in the current working dir, but::

        filefind('~/myfile.txt')

    Will find the file in the users home directory.  This function does not
    automatically try any paths, such as the cwd or the user's home directory.

    Parameters
    ----------
    filename : str
        The filename to look for.
    path_dirs : str, None or sequence of str
        The sequence of paths to look for the file in.  If None, the filename
        need to be absolute or be in the cwd.  If a string, the string is
        put into a sequence and the searched.  If a sequence, walk through
        each element and join with ``filename``, calling :func:`expandvars`
        and :func:`expanduser` before testing for existence.

    Returns
    -------
    Raises :exc:`IOError` or returns absolute path to file.
    """

    # If paths are quoted, abspath gets confused, strip them...
    filename = filename.strip('"').strip("'")
    # If the input is an absolute path, just check it exists
    if os.path.isabs(filename) and os.path.isfile(filename):
        return filename

    if path_dirs is None:
        path_dirs = ("",)
    elif isinstance(path_dirs, str):
        path_dirs = (path_dirs,)

    for path in path_dirs:
        if path == '.': path = os.getcwd()
        testname = expand_path(os.path.join(path, filename))
        if os.path.isfile(testname):
            return os.path.abspath(testname)

    raise IOError("File %r does not exist in any of the search paths: %r" %
                  (filename, path_dirs) )


class HomeDirError(Exception):
    pass


def get_home_dir(require_writable=False):
    """Return the 'home' directory, as a unicode string.

    Uses os.path.expanduser('~'), and checks for writability.

    See stdlib docs for how this is determined.
    $HOME is first priority on *ALL* platforms.

    Parameters
    ----------

    require_writable : bool [default: False]
        if True:
            guarantees the return value is a writable directory, otherwise
            raises HomeDirError
        if False:
            The path is resolved, but it is not guaranteed to exist or be writable.
    """

    homedir = os.path.expanduser('~')
    # Next line will make things work even when /home/ is a symlink to
    # /usr/home as it is on FreeBSD, for example
    homedir = os.path.realpath(homedir)

    if not _writable_dir(homedir) and os.name == 'nt':
        # expanduser failed, use the registry to get the 'My Documents' folder.
        try:
            try:
                import winreg as wreg  # Py 3
            except ImportError:
                import _winreg as wreg  # Py 2
            key = wreg.OpenKey(
                wreg.HKEY_CURRENT_USER,
                "Software\Microsoft\Windows\CurrentVersion\Explorer\Shell Folders"
            )
            homedir = wreg.QueryValueEx(key,'Personal')[0]
            key.Close()
        except:
            pass

    if (not require_writable) or _writable_dir(homedir):
        return py3compat.cast_unicode(homedir, fs_encoding)
    else:
        raise HomeDirError('%s is not a writable dir, '
                'set $HOME environment variable to override' % homedir)

def get_xdg_dir():
    """Return the XDG_CONFIG_HOME, if it is defined and exists, else None.

    This is only for non-OS X posix (Linux,Unix,etc.) systems.
    """

    env = os.environ

    if os.name == 'posix' and sys.platform != 'darwin':
        # Linux, Unix, AIX, etc.
        # use ~/.config if empty OR not set
        xdg = env.get("XDG_CONFIG_HOME", None) or os.path.join(get_home_dir(), '.config')
        if xdg and _writable_dir(xdg):
            return py3compat.cast_unicode(xdg, fs_encoding)

    return None


def get_xdg_cache_dir():
    """Return the XDG_CACHE_HOME, if it is defined and exists, else None.

    This is only for non-OS X posix (Linux,Unix,etc.) systems.
    """

    env = os.environ

    if os.name == 'posix' and sys.platform != 'darwin':
        # Linux, Unix, AIX, etc.
        # use ~/.cache if empty OR not set
        xdg = env.get("XDG_CACHE_HOME", None) or os.path.join(get_home_dir(), '.cache')
        if xdg and _writable_dir(xdg):
            return py3compat.cast_unicode(xdg, fs_encoding)

    return None


@undoc
def get_ipython_dir():
    warn("get_ipython_dir has moved to the IPython.paths module since IPython 4.0.", stacklevel=2)
    from IPython.paths import get_ipython_dir
    return get_ipython_dir()

@undoc
def get_ipython_cache_dir():
    warn("get_ipython_cache_dir has moved to the IPython.paths module since IPython 4.0.", stacklevel=2)
    from IPython.paths import get_ipython_cache_dir
    return get_ipython_cache_dir()

@undoc
def get_ipython_package_dir():
    warn("get_ipython_package_dir has moved to the IPython.paths module since IPython 4.0.", stacklevel=2)
    from IPython.paths import get_ipython_package_dir
    return get_ipython_package_dir()

@undoc
def get_ipython_module_path(module_str):
    warn("get_ipython_module_path has moved to the IPython.paths module since IPython 4.0.", stacklevel=2)
    from IPython.paths import get_ipython_module_path
    return get_ipython_module_path(module_str)

@undoc
def locate_profile(profile='default'):
    warn("locate_profile has moved to the IPython.paths module since IPython 4.0.", stacklevel=2)
    from IPython.paths import locate_profile
    return locate_profile(profile=profile)

def expand_path(s):
    """Expand $VARS and ~names in a string, like a shell

    :Examples:

       In [2]: os.environ['FOO']='test'

       In [3]: expand_path('variable FOO is $FOO')
       Out[3]: 'variable FOO is test'
    """
    # This is a pretty subtle hack. When expand user is given a UNC path
    # on Windows (\\server\share$\%username%), os.path.expandvars, removes
    # the $ to get (\\server\share\%username%). I think it considered $
    # alone an empty var. But, we need the $ to remains there (it indicates
    # a hidden share).
    if os.name=='nt':
        s = s.replace('$\\', 'IPYTHON_TEMP')
    s = os.path.expandvars(os.path.expanduser(s))
    if os.name=='nt':
        s = s.replace('IPYTHON_TEMP', '$\\')
    return s


def unescape_glob(string):
    """Unescape glob pattern in `string`."""
    def unescape(s):
        for pattern in '*[]!?':
            s = s.replace(r'\{0}'.format(pattern), pattern)
        return s
    return '\\'.join(map(unescape, string.split('\\\\')))


def shellglob(args):
    """
    Do glob expansion for each element in `args` and return a flattened list.

    Unmatched glob pattern will remain as-is in the returned list.

    """
    expanded = []
    # Do not unescape backslash in Windows as it is interpreted as
    # path separator:
    unescape = unescape_glob if sys.platform != 'win32' else lambda x: x
    for a in args:
        expanded.extend(glob.glob(a) or [unescape(a)])
    return expanded


def target_outdated(target,deps):
    """Determine whether a target is out of date.

    target_outdated(target,deps) -> 1/0

    deps: list of filenames which MUST exist.
    target: single filename which may or may not exist.

    If target doesn't exist or is older than any file listed in deps, return
    true, otherwise return false.
    """
    try:
        target_time = os.path.getmtime(target)
    except os.error:
        return 1
    for dep in deps:
        dep_time = os.path.getmtime(dep)
        if dep_time > target_time:
            #print "For target",target,"Dep failed:",dep # dbg
            #print "times (dep,tar):",dep_time,target_time # dbg
            return 1
    return 0


def target_update(target,deps,cmd):
    """Update a target with a given command given a list of dependencies.

    target_update(target,deps,cmd) -> runs cmd if target is outdated.

    This is just a wrapper around target_outdated() which calls the given
    command if target is outdated."""

    if target_outdated(target,deps):
        system(cmd)


ENOLINK = 1998

def link(src, dst):
    """Hard links ``src`` to ``dst``, returning 0 or errno.

    Note that the special errno ``ENOLINK`` will be returned if ``os.link`` isn't
    supported by the operating system.
    """

    if not hasattr(os, "link"):
        return ENOLINK
    link_errno = 0
    try:
        os.link(src, dst)
    except OSError as e:
        link_errno = e.errno
    return link_errno


def link_or_copy(src, dst):
    """Attempts to hardlink ``src`` to ``dst``, copying if the link fails.

    Attempts to maintain the semantics of ``shutil.copy``.

    Because ``os.link`` does not overwrite files, a unique temporary file
    will be used if the target already exists, then that file will be moved
    into place.
    """

    if os.path.isdir(dst):
        dst = os.path.join(dst, os.path.basename(src))

    link_errno = link(src, dst)
    if link_errno == errno.EEXIST:
        if os.stat(src).st_ino == os.stat(dst).st_ino:
            # dst is already a hard link to the correct file, so we don't need
            # to do anything else. If we try to link and rename the file
            # anyway, we get duplicate files - see http://bugs.python.org/issue21876
            return

        new_dst = dst + "-temp-%04X" %(random.randint(1, 16**4), )
        try:
            link_or_copy(src, new_dst)
        except:
            try:
                os.remove(new_dst)
            except OSError:
                pass
            raise
        os.rename(new_dst, dst)
    elif link_errno != 0:
        # Either link isn't supported, or the filesystem doesn't support
        # linking, or 'src' and 'dst' are on different filesystems.
        shutil.copy(src, dst)

def ensure_dir_exists(path, mode=0o755):
    """ensure that a directory exists

    If it doesn't exist, try to create it and protect against a race condition
    if another process is doing the same.

    The default permissions are 755, which differ from os.makedirs default of 777.
    """
    if not os.path.exists(path):
        try:
            os.makedirs(path, mode=mode)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise
    elif not os.path.isdir(path):
        raise IOError("%r exists but is not a directory" % path)
