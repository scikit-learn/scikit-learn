"""
Shim to maintain backwards compatibility with old frontend imports.

We have moved all contents of the old `frontend` subpackage into top-level
subpackages (`html`, `qt` and `terminal`), and flattened the notebook into
just `IPython.html`, formerly `IPython.frontend.html.notebook`.

This will let code that was making `from IPython.frontend...` calls continue
working, though a warning will be printed.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
from warnings import warn

from IPython.utils.shimmodule import ShimModule, ShimWarning

warn("The top-level `frontend` package has been deprecated since IPython 1.0. "
     "All its subpackages have been moved to the top `IPython` level.", ShimWarning)

# Unconditionally insert the shim into sys.modules so that further import calls
# trigger the custom attribute access above

sys.modules['IPython.frontend.html.notebook'] = ShimModule(
    src='IPython.frontend.html.notebook', mirror='IPython.html')
sys.modules['IPython.frontend'] = ShimModule(
    src='IPython.frontend', mirror='IPython')
