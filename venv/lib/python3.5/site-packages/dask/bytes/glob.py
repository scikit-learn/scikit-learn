# The glob functionality here copies (and heavily modifies) code from the
# `glob` module to allow for use with arrow's hdfs. These functions are subject
# to the license found at https://docs.python.org/3/license.html, which is also
# included below:
#
#                PSF LICENSE AGREEMENT FOR PYTHON 3.6.4
#                ======================================
#
# 1. This LICENSE AGREEMENT is between the Python Software Foundation ("PSF"),
#    and the Individual or Organization ("Licensee") accessing and otherwise
#    using Python 3.6.4 software in source or binary form and its associated
#    documentation.
#
# 2. Subject to the terms and conditions of this License Agreement, PSF hereby
#    grants Licensee a nonexclusive, royalty-free, world-wide license to
#    reproduce, analyze, test, perform and/or display publicly, prepare
#    derivative works, distribute, and otherwise use Python 3.6.4 alone or in
#    any derivative version, provided, however, that PSF's License Agreement
#    and PSF's notice of copyright, i.e., "Copyright c 2001-2016 Python
#    Software Foundation; All Rights Reserved" are retained in Python 3.6.4
#    alone or in any derivative version prepared by Licensee.
#
# 3. In the event Licensee prepares a derivative work that is based on or
#    incorporates Python 3.6.4 or any part thereof, and wants to make the
#    derivative work available to others as provided herein, then Licensee
#    hereby agrees to include in any such work a brief summary of the changes
#    made to Python 3.6.4.
#
# 4. PSF is making Python 3.6.4 available to Licensee on an "AS IS" basis.  PSF
#    MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.  BY WAY OF
#    EXAMPLE, BUT NOT LIMITATION, PSF MAKES NO AND DISCLAIMS ANY REPRESENTATION
#    OR WARRANTY OF MERCHANTABILITY OR FITNESS FOR ANY PARTICULAR PURPOSE OR
#    THAT THE USE OF PYTHON 3.6.4 WILL NOT INFRINGE ANY THIRD PARTY RIGHTS.
#
# 5. PSF SHALL NOT BE LIABLE TO LICENSEE OR ANY OTHER USERS OF PYTHON 3.6.4 FOR
#    ANY INCIDENTAL, SPECIAL, OR CONSEQUENTIAL DAMAGES OR LOSS AS A RESULT OF
#    MODIFYING, DISTRIBUTING, OR OTHERWISE USING PYTHON 3.6.4, OR ANY
#    DERIVATIVE THEREOF, EVEN IF ADVISED OF THE POSSIBILITY THEREOF.
#
# 6. This License Agreement will automatically terminate upon a material breach
#    of its terms and conditions.
#
# 7. Nothing in this License Agreement shall be deemed to create any
#    relationship of agency, partnership, or joint venture between PSF and
#    Licensee.  This License Agreement does not grant permission to use PSF
#    trademarks or trade name in a trademark sense to endorse or promote
#    products or services of Licensee, or any third party.
#
# 8. By copying, installing or otherwise using Python 3.6.4, Licensee agrees to
#    be bound by the terms and conditions of this License Agreement.
#
# These functions are under copyright by the Python Software Foundation
#
#    Copyright 2001-2018 Python Software Foundation; All Rights Reserved

import fnmatch
import re


def generic_glob(fs, path_impl, pathname):
    """A filesystem agnostic glob implemention.

    Parameters
    ----------
    fs : filesystem
        The filesystem to search.
    path_impl : os.path like
        The path module implementation to use. Designed to pass in
        ``posixpath`` or ``ntpath`` modules directly.
    pathname : str
        The path or pattern to glob

    Returns
    -------
    paths : list
        A list of paths matching the given path or pattern.
    """
    dirname, basename = path_impl.split(pathname)
    if not dirname:
        raise ValueError("glob pattern must be an absolute path")
    if not _has_magic(pathname):
        if (not basename and _safe_isdir(fs, dirname) or
                basename and fs.exists(pathname)):
            return [pathname]
        return []
    if basename and _has_magic(dirname):
        # Directory is a pattern, collect all matching directories
        dirs = [d for d in generic_glob(fs, path_impl, dirname)
                if _safe_isdir(fs, d)]
    else:
        # No basename (pattern ends in `/`, must match directories only)
        # or no magic in dirname (use dirname directly)
        dirs = [dirname] if _safe_isdir(fs, dirname) else []
    glob_in_dir = _glob_pattern if _has_magic(basename) else _glob_path
    return [path_impl.join(dirname2, name)
            for dirname2 in dirs
            for name in glob_in_dir(fs, path_impl, dirname2, basename)]


def _safe_isdir(fs, dirname):
    try:
        return fs.isdir(dirname)
    except OSError:
        # pyarrow isdir raises if the directory doesn't exist
        return False


def _glob_pattern(fs, path_impl, dirname, pattern):
    names = [path_impl.split(f)[1] for f in fs.ls(dirname)]
    if not _ishidden(pattern):
        names = [x for x in names if not _ishidden(x)]
    return fnmatch.filter(names, pattern)


def _glob_path(fs, path_impl, dirname, basename):
    if (not basename and _safe_isdir(fs, dirname) or
            basename and fs.exists(path_impl.join(dirname, basename))):
        return [basename]
    return []


_magic_check = re.compile('([*?[])')


def _has_magic(s):
    return _magic_check.search(s) is not None


def _ishidden(path):
    return path[0] == '.'
