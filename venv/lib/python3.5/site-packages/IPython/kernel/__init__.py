"""
Shim to maintain backwards compatibility with old IPython.kernel imports.
"""
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
from warnings import warn

from IPython.utils.shimmodule import ShimModule, ShimWarning

warn("The `IPython.kernel` package has been deprecated since IPython 4.0."
     "You should import from ipykernel or jupyter_client instead.", ShimWarning)


# zmq subdir is gone
sys.modules['IPython.kernel.zmq.session'] = ShimModule(
    src='IPython.kernel.zmq.session', mirror='jupyter_client.session')
sys.modules['IPython.kernel.zmq'] = ShimModule(
    src='IPython.kernel.zmq', mirror='ipykernel')

for pkg in ('comm', 'inprocess'):
    src = 'IPython.kernel.%s' % pkg
    sys.modules[src] = ShimModule(src=src, mirror='ipykernel.%s' % pkg)

for pkg in ('ioloop', 'blocking'):
    src = 'IPython.kernel.%s' % pkg
    sys.modules[src] = ShimModule(src=src, mirror='jupyter_client.%s' % pkg)

# required for `from IPython.kernel import PKG`
from ipykernel import comm, inprocess
from jupyter_client import ioloop, blocking
# public API
from ipykernel.connect import *
from jupyter_client import *
