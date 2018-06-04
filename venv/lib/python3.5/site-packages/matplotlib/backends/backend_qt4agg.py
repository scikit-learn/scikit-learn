"""
Render to qt from agg
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from .backend_qt5agg import (
    _BackendQT5Agg, FigureCanvasQTAgg, FigureManagerQT, NavigationToolbar2QT)


@_BackendQT5Agg.export
class _BackendQT4Agg(_BackendQT5Agg):
    pass
