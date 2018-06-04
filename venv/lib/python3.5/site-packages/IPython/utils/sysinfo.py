# encoding: utf-8
"""
Utilities for getting information about IPython and the system it's running in.
"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os
import platform
import pprint
import sys
import subprocess

from IPython.core import release
from IPython.utils import py3compat, _sysinfo, encoding

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def pkg_commit_hash(pkg_path):
    """Get short form of commit hash given directory `pkg_path`

    We get the commit hash from (in order of preference):

    * IPython.utils._sysinfo.commit
    * git output, if we are in a git repository

    If these fail, we return a not-found placeholder tuple

    Parameters
    ----------
    pkg_path : str
       directory containing package
       only used for getting commit from active repo

    Returns
    -------
    hash_from : str
       Where we got the hash from - description
    hash_str : str
       short form of hash
    """
    # Try and get commit from written commit text file
    if _sysinfo.commit:
        return "installation", _sysinfo.commit

    # maybe we are in a repository
    proc = subprocess.Popen('git rev-parse --short HEAD',
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            cwd=pkg_path, shell=True)
    repo_commit, _ = proc.communicate()
    if repo_commit:
        return 'repository', repo_commit.strip().decode('ascii')
    return '(none found)', '<not found>'


def pkg_info(pkg_path):
    """Return dict describing the context of this package

    Parameters
    ----------
    pkg_path : str
       path containing __init__.py for package

    Returns
    -------
    context : dict
       with named parameters of interest
    """
    src, hsh = pkg_commit_hash(pkg_path)
    return dict(
        ipython_version=release.version,
        ipython_path=pkg_path,
        commit_source=src,
        commit_hash=hsh,
        sys_version=sys.version,
        sys_executable=sys.executable,
        sys_platform=sys.platform,
        platform=platform.platform(),
        os_name=os.name,
        default_encoding=encoding.DEFAULT_ENCODING,
        )

def get_sys_info():
    """Return useful information about IPython and the system, as a dict."""
    p = os.path
    path = p.realpath(p.dirname(p.abspath(p.join(__file__, '..'))))
    return pkg_info(path)

@py3compat.doctest_refactor_print
def sys_info():
    """Return useful information about IPython and the system, as a string.

    Examples
    --------
    ::
    
        In [2]: print sys_info()
        {'commit_hash': '144fdae',      # random
         'commit_source': 'repository',
         'ipython_path': '/home/fperez/usr/lib/python2.6/site-packages/IPython',
         'ipython_version': '0.11.dev',
         'os_name': 'posix',
         'platform': 'Linux-2.6.35-22-generic-i686-with-Ubuntu-10.10-maverick',
         'sys_executable': '/usr/bin/python',
         'sys_platform': 'linux2',
         'sys_version': '2.6.6 (r266:84292, Sep 15 2010, 15:52:39) \\n[GCC 4.4.5]'}
    """
    return pprint.pformat(get_sys_info())

def _num_cpus_unix():
    """Return the number of active CPUs on a Unix system."""
    return os.sysconf("SC_NPROCESSORS_ONLN")


def _num_cpus_darwin():
    """Return the number of active CPUs on a Darwin system."""
    p = subprocess.Popen(['sysctl','-n','hw.ncpu'],stdout=subprocess.PIPE)
    return p.stdout.read()


def _num_cpus_windows():
    """Return the number of active CPUs on a Windows system."""
    return os.environ.get("NUMBER_OF_PROCESSORS")


def num_cpus():
   """Return the effective number of CPUs in the system as an integer.

   This cross-platform function makes an attempt at finding the total number of
   available CPUs in the system, as returned by various underlying system and
   python calls.

   If it can't find a sensible answer, it returns 1 (though an error *may* make
   it return a large positive number that's actually incorrect).
   """

   # Many thanks to the Parallel Python project (http://www.parallelpython.com)
   # for the names of the keys we needed to look up for this function.  This
   # code was inspired by their equivalent function.

   ncpufuncs = {'Linux':_num_cpus_unix,
                'Darwin':_num_cpus_darwin,
                'Windows':_num_cpus_windows
                }

   ncpufunc = ncpufuncs.get(platform.system(),
                            # default to unix version (Solaris, AIX, etc)
                            _num_cpus_unix)

   try:
       ncpus = max(1,int(ncpufunc()))
   except:
       ncpus = 1
   return ncpus

