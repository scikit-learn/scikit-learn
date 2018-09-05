# flake8: noqa
###############################################################################
# Compat file to import the correct modules for each platform and python
# version.
#
# author: Thomas Moreau and Olivier grisel
#
import sys

if sys.version_info[:2] >= (3, 3):
    import queue
else:
    import Queue as queue

from pickle import PicklingError

if sys.version_info >= (3, 4):
    from multiprocessing.process import BaseProcess
else:
    from multiprocessing.process import Process as BaseProcess

# Platform specific compat
if sys.platform == "win32":
    from .compat_win32 import *
else:
    from .compat_posix import *
