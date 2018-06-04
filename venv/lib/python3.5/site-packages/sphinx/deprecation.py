# -*- coding: utf-8 -*-
"""
    sphinx.deprecation
    ~~~~~~~~~~~~~~~~~~

    Sphinx deprecation classes and utilities.

    :copyright: Copyright 2007-2018 by the Sphinx team, see AUTHORS.
    :license: BSD, see LICENSE for details.
"""


class RemovedInSphinx18Warning(DeprecationWarning):
    pass


class RemovedInSphinx19Warning(PendingDeprecationWarning):
    pass


class RemovedInSphinx20Warning(PendingDeprecationWarning):
    pass


RemovedInNextVersionWarning = RemovedInSphinx18Warning
