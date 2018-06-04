"""
Shim to maintain backwards compatibility with old IPython.qt imports.
"""
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
from warnings import warn

from IPython.utils.shimmodule import ShimModule, ShimWarning

warn("The `IPython.qt` package has been deprecated since IPython 4.0. "
     "You should import from qtconsole instead.", ShimWarning)

# Unconditionally insert the shim into sys.modules so that further import calls
# trigger the custom attribute access above

_console = sys.modules['IPython.qt.console'] = ShimModule(
    src='IPython.qt.console', mirror='qtconsole')

_qt = ShimModule(src='IPython.qt', mirror='qtconsole')

_qt.console = _console
sys.modules['IPython.qt'] = _qt
