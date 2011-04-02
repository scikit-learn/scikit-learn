"""
Disk management utilities.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org> 
# Copyright (c) 2010 Gael Varoquaux
# License: BSD Style, 3 clauses.


import platform
import os
import shutil

def disk_used(path):
    """ Return the disk usage in a directory. 
    """
    size = 0
    for file in os.listdir(path) + ['.']:
        stat =  os.stat(os.path.join(path, file))
        if hasattr(stat, 'st_blocks'):
            size += stat.st_blocks * 512
        else:
            # on some platform st_blocks is not available (e.g., Windows)
            # approximate by rounding to next multiple of 512
            size += (stat.st_size // 512 + 1) * 512;
    # We need to convert to int to avoid having longs on some systems (we
    # don't want longs to avoid problems we SQLite)
    return int(size/1024.)


def memstr_to_kbytes(text):
    """ Convert a memory text to it's value in kilobytes.
    """
    kilo = 1024
    units = dict(K=1, M=kilo, G=kilo**2)
    try:
        size = int(units[text[-1]]*float(text[:-1]))
    except (KeyError, ValueError):
        raise ValueError(
                "Invalid literal for size give: %s (type %s) should be "
                "alike '10G', '500M', '50K'." % (text, type(text))
                )
    return size

def rm_subdirs(path, onerror=None):
    """Remove all subdirectories in this path.

    If onerror is set, it is called to handle the error with arguments (func,
    path, exc_info) where func is os.listdir, os.remove, or os.rmdir;
    path is the argument to that function that caused it to fail; and
    exc_info is a tuple returned by sys.exc_info().  If ignore_errors
    is false and onerror is None, an exception is raised.
    """

    # NOTE this code is adapted from the one in shutil.rmtree, and is
    # just as fast

    names = []
    try:
        names = os.listdir(path)
    except os.error, err:
        onerror(os.listdir, path, sys.exc_info())
        
    for name in names:
        fullname = os.path.join(path, name)
        if os.path.isdir(fullname):
            shutil.rmtree(fullname, False, onerror)
