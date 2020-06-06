from __future__ import absolute_import

import sys
import pickle


from .cloudpickle import *
if sys.version_info[:2] >= (3, 8):
    from .cloudpickle_fast import CloudPickler, dumps, dump

__version__ = '1.4.1'
