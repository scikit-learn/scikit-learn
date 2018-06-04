"""TemporaryDirectory class, copied from Python 3

This is copied from the stdlib and will be standard in Python 3.2 and onwards.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

from __future__ import print_function

import os as _os
import warnings as _warnings
import sys as _sys

# This code should only be used in Python versions < 3.2, since after that we
# can rely on the stdlib itself.
try:
    from tempfile import TemporaryDirectory

except ImportError:
    from tempfile import mkdtemp, template

    class TemporaryDirectory(object):
        """Create and return a temporary directory.  This has the same
        behavior as mkdtemp but can be used as a context manager.  For
        example:

            with TemporaryDirectory() as tmpdir:
                ...

        Upon exiting the context, the directory and everthing contained
        in it are removed.
        """

        def __init__(self, suffix="", prefix=template, dir=None):
            self.name = mkdtemp(suffix, prefix, dir)
            self._closed = False

        def __enter__(self):
            return self.name

        def cleanup(self, _warn=False):
            if self.name and not self._closed:
                try:
                    self._rmtree(self.name)
                except (TypeError, AttributeError) as ex:
                    # Issue #10188: Emit a warning on stderr
                    # if the directory could not be cleaned
                    # up due to missing globals
                    if "None" not in str(ex):
                        raise
                    print("ERROR: {!r} while cleaning up {!r}".format(ex, self,),
                          file=_sys.stderr)
                    return
                self._closed = True
                if _warn:
                    self._warn("Implicitly cleaning up {!r}".format(self),
                               Warning)

        def __exit__(self, exc, value, tb):
            self.cleanup()

        def __del__(self):
            # Issue a ResourceWarning if implicit cleanup needed
            self.cleanup(_warn=True)


        # XXX (ncoghlan): The following code attempts to make
        # this class tolerant of the module nulling out process
        # that happens during CPython interpreter shutdown
        # Alas, it doesn't actually manage it. See issue #10188
        _listdir = staticmethod(_os.listdir)
        _path_join = staticmethod(_os.path.join)
        _isdir = staticmethod(_os.path.isdir)
        _remove = staticmethod(_os.remove)
        _rmdir = staticmethod(_os.rmdir)
        _os_error = _os.error
        _warn = _warnings.warn

        def _rmtree(self, path):
            # Essentially a stripped down version of shutil.rmtree.  We can't
            # use globals because they may be None'ed out at shutdown.
            for name in self._listdir(path):
                fullname = self._path_join(path, name)
                try:
                    isdir = self._isdir(fullname)
                except self._os_error:
                    isdir = False
                if isdir:
                    self._rmtree(fullname)
                else:
                    try:
                        self._remove(fullname)
                    except self._os_error:
                        pass
            try:
                self._rmdir(path)
            except self._os_error:
                pass

# extra temp-dir-related context managers

class NamedFileInTemporaryDirectory(object):

    def __init__(self, filename, mode='w+b', bufsize=-1, **kwds):
        """
        Open a file named `filename` in a temporary directory.

        This context manager is preferred over `NamedTemporaryFile` in
        stdlib `tempfile` when one needs to reopen the file.

        Arguments `mode` and `bufsize` are passed to `open`.
        Rest of the arguments are passed to `TemporaryDirectory`.

        """
        self._tmpdir = TemporaryDirectory(**kwds)
        path = _os.path.join(self._tmpdir.name, filename)
        self.file = open(path, mode, bufsize)

    def cleanup(self):
        self.file.close()
        self._tmpdir.cleanup()

    __del__ = cleanup

    def __enter__(self):
        return self.file

    def __exit__(self, type, value, traceback):
        self.cleanup()


class TemporaryWorkingDirectory(TemporaryDirectory):
    """
    Creates a temporary directory and sets the cwd to that directory.
    Automatically reverts to previous cwd upon cleanup.
    Usage example:

        with TemporaryWorkingDirectory() as tmpdir:
            ...
    """
    def __enter__(self):
        self.old_wd = _os.getcwd()
        _os.chdir(self.name)
        return super(TemporaryWorkingDirectory, self).__enter__()

    def __exit__(self, exc, value, tb):
        _os.chdir(self.old_wd)
        return super(TemporaryWorkingDirectory, self).__exit__(exc, value, tb)

