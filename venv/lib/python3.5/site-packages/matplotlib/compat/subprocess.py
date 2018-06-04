"""
A replacement wrapper around the subprocess module, with a number of
work-arounds:
- Provides a stub implementation of subprocess members on Google App Engine
  (which are missing in subprocess).
- Use subprocess32, backport from python 3.2 on Linux/Mac work-around for
  https://github.com/matplotlib/matplotlib/issues/5314

Instead of importing subprocess, other modules should use this as follows:

from matplotlib.compat import subprocess

This module is safe to import from anywhere within matplotlib.
"""

from __future__ import absolute_import    # Required to import subprocess
from __future__ import print_function
import os
import sys
if os.name == 'posix' and sys.version_info[0] < 3:
    # work around for https://github.com/matplotlib/matplotlib/issues/5314
    try:
        import subprocess32 as subprocess
    except ImportError:
        import subprocess
else:
    import subprocess

__all__ = ['Popen', 'PIPE', 'STDOUT', 'check_output', 'CalledProcessError']


if hasattr(subprocess, 'Popen'):
    Popen = subprocess.Popen
    # Assume that it also has the other constants.
    PIPE = subprocess.PIPE
    STDOUT = subprocess.STDOUT
    CalledProcessError = subprocess.CalledProcessError
    check_output = subprocess.check_output
else:
    # In restricted environments (such as Google App Engine), these are
    # non-existent. Replace them with dummy versions that always raise OSError.
    def Popen(*args, **kwargs):
        raise OSError("subprocess.Popen is not supported")

    def check_output(*args, **kwargs):
        raise OSError("subprocess.check_output is not supported")
    PIPE = -1
    STDOUT = -2
    # There is no need to catch CalledProcessError. These stubs cannot raise
    # it. None in an except clause will simply not match any exceptions.
    CalledProcessError = None
