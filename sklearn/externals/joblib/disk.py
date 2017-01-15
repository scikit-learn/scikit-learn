"""
Disk management utilities.
"""

# Authors: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
#          Lars Buitinck
# Copyright (c) 2010 Gael Varoquaux
# License: BSD Style, 3 clauses.


import errno
import os
import shutil
import sys
import time


def disk_used(path):
    """ Return the disk usage in a directory."""
    size = 0
    for file in os.listdir(path) + ['.']:
        stat = os.stat(os.path.join(path, file))
        if hasattr(stat, 'st_blocks'):
            size += stat.st_blocks * 512
        else:
            # on some platform st_blocks is not available (e.g., Windows)
            # approximate by rounding to next multiple of 512
            size += (stat.st_size // 512 + 1) * 512
    # We need to convert to int to avoid having longs on some systems (we
    # don't want longs to avoid problems we SQLite)
    return int(size / 1024.)


def memstr_to_bytes(text):
    """ Convert a memory text to its value in bytes.
    """
    kilo = 1024
    units = dict(K=kilo, M=kilo ** 2, G=kilo ** 3)
    try:
        size = int(units[text[-1]] * float(text[:-1]))
    except (KeyError, ValueError):
        raise ValueError(
            "Invalid literal for size give: %s (type %s) should be "
            "alike '10G', '500M', '50K'." % (text, type(text)))
    return size


def mkdirp(d):
    """Ensure directory d exists (like mkdir -p on Unix)
    No guarantee that the directory is writable.
    """
    try:
        os.makedirs(d)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


# if a rmtree operation fails in rm_subdirs, wait for this much time (in secs),
# then retry once. if it still fails, raise the exception
RM_SUBDIRS_RETRY_TIME = 0.1


def rm_subdirs(path, onerror=None):
    """Remove all subdirectories in this path.

    The directory indicated by `path` is left in place, and its subdirectories
    are erased.

    If onerror is set, it is called to handle the error with arguments (func,
    path, exc_info) where func is os.listdir, os.remove, or os.rmdir;
    path is the argument to that function that caused it to fail; and
    exc_info is a tuple returned by sys.exc_info().  If onerror is None,
    an exception is raised.
    """

    # NOTE this code is adapted from the one in shutil.rmtree, and is
    # just as fast

    names = []
    try:
        names = os.listdir(path)
    except os.error as err:
        if onerror is not None:
            onerror(os.listdir, path, sys.exc_info())
        else:
            raise

    for name in names:
        fullname = os.path.join(path, name)
        if os.path.isdir(fullname):
            if onerror is not None:
                shutil.rmtree(fullname, False, onerror)
            else:
                # allow the rmtree to fail once, wait and re-try.
                # if the error is raised again, fail
                err_count = 0
                while True:
                    try:
                        shutil.rmtree(fullname, False, None)
                        break
                    except os.error:
                        if err_count > 0:
                            raise
                        err_count += 1
                        time.sleep(RM_SUBDIRS_RETRY_TIME)
