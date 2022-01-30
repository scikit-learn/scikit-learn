"""
    sphinx.util.compat
    ~~~~~~~~~~~~~~~~~~

    modules for backward compatibility

    :copyright: Copyright 2007-2022 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""

import sys
from typing import TYPE_CHECKING, Any, Dict

if TYPE_CHECKING:
    from sphinx.application import Sphinx


def register_application_for_autosummary(app: "Sphinx") -> None:
    """Register application object to autosummary module.

    Since Sphinx-1.7, documenters and attrgetters are registered into
    application object.  As a result, the arguments of
    ``get_documenter()`` has been changed.  To keep compatibility,
    this handler registers application object to the module.
    """
    if 'sphinx.ext.autosummary' in sys.modules:
        from sphinx.ext import autosummary
        if hasattr(autosummary, '_objects'):
            autosummary._objects['_app'] = app  # type: ignore
        else:
            autosummary._app = app  # type: ignore


def setup(app: "Sphinx") -> Dict[str, Any]:
    app.connect('builder-inited', register_application_for_autosummary, priority=100)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
