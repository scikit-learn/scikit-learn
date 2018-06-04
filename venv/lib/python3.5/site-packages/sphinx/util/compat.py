# -*- coding: utf-8 -*-
"""
    sphinx.util.compat
    ~~~~~~~~~~~~~~~~~~

    modules for backward compatibility

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""
from __future__ import absolute_import

import sys


if False:
    # For type annotation
    from typing import Any, Dict  # NOQA
    from sphinx.application import Sphinx  # NOQA


def register_application_for_autosummary(app):
    # type: (Sphinx) -> None
    """Register application object to autosummary module.

    Since Sphinx-1.7, documenters and attrgetters are registered into
    applicaiton object.  As a result, the arguments of
    ``get_documenter()`` has been changed.  To keep compatibility,
    this handler registers application object to the module.
    """
    if 'sphinx.ext.autosummary' in sys.modules:
        from sphinx.ext import autosummary
        autosummary._app = app


def setup(app):
    # type: (Sphinx) -> Dict[unicode, Any]
    app.connect('builder-inited', register_application_for_autosummary)

    return {
        'version': 'builtin',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
