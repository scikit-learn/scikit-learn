"""
Shim to maintain backwards compatibility with old IPython.html imports.
"""
# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.

import sys
from warnings import warn

from IPython.utils.shimmodule import ShimModule, ShimWarning

warn("The `IPython.html` package has been deprecated since IPython 4.0. "
     "You should import from `notebook` instead. "
     "`IPython.html.widgets` has moved to `ipywidgets`.", ShimWarning)

_widgets = sys.modules['IPython.html.widgets'] = ShimModule(
    src='IPython.html.widgets', mirror='ipywidgets')

_html = ShimModule(
    src='IPython.html', mirror='notebook')

# hook up widgets
_html.widgets = _widgets
sys.modules['IPython.html'] = _html

if __name__ == '__main__':
    from notebook import notebookapp as app
    app.launch_new_instance()
