from __future__ import print_function
"""Utilities to find the site and prefix information of a Python executable, which may be Python 2.

This file MUST remain compatible with Python 2. Since we cannot make any assumptions about the
Python being executed, this module should not use *any* dependencies outside of the standard
library found in Python 2. This file is run each mypy run, so it should be kept as fast as
possible.
"""
import site
import sys

if __name__ == '__main__':
    sys.path = sys.path[1:]  # we don't want to pick up mypy.types

MYPY = False
if MYPY:
    from typing import List, Tuple


def getprefixes():
    # type: () -> Tuple[str, str]
    return sys.base_prefix, sys.prefix


def getsitepackages():
    # type: () -> List[str]
    res = []
    if hasattr(site, 'getsitepackages'):
        res.extend(site.getsitepackages())

        if hasattr(site, 'getusersitepackages') and site.ENABLE_USER_SITE:
            res.insert(0, site.getusersitepackages())
    else:
        from distutils.sysconfig import get_python_lib
        res = [get_python_lib()]
    return res


if __name__ == '__main__':
    if sys.argv[-1] == 'getsitepackages':
        print(repr(getsitepackages()))
    elif sys.argv[-1] == 'getprefixes':
        print(repr(getprefixes()))
    else:
        print("ERROR: incorrect argument to pyinfo.py.", file=sys.stderr)
        sys.exit(1)
