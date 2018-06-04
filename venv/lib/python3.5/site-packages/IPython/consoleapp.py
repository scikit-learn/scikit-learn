"""
Shim to maintain backwards compatibility with old IPython.consoleapp imports.
"""
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

from warnings import warn

warn("The `IPython.consoleapp` package has been deprecated. "
     "You should import from jupyter_client.consoleapp instead.")

from jupyter_client.consoleapp import *
