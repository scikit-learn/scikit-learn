"""
GTK+ Matplotlib interface using cairo (not GDK) drawing operations.
Author: Steve Chaplin
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import gtk
if gtk.pygtk_version < (2, 7, 0):
    import cairo.gtk

from matplotlib import cbook
from matplotlib.backends import backend_cairo
from matplotlib.backends.backend_gtk import *
from matplotlib.backends.backend_gtk import _BackendGTK

backend_version = ('PyGTK(%d.%d.%d) ' % gtk.pygtk_version
                   + 'Pycairo(%s)' % backend_cairo.backend_version)


class RendererGTKCairo (backend_cairo.RendererCairo):
    if gtk.pygtk_version >= (2,7,0):
        def set_pixmap (self, pixmap):
            self.gc.ctx = pixmap.cairo_create()
    else:
        def set_pixmap (self, pixmap):
            self.gc.ctx = cairo.gtk.gdk_cairo_create (pixmap)


class FigureCanvasGTKCairo(backend_cairo.FigureCanvasCairo, FigureCanvasGTK):
    filetypes = FigureCanvasGTK.filetypes.copy()
    filetypes.update(backend_cairo.FigureCanvasCairo.filetypes)

    def __init__(self, *args, **kwargs):
        warn_deprecated('2.2',
                        message=('The GTKCairo backend is deprecated. It is '
                                 'untested and will be removed in Matplotlib '
                                 '3.0. Use the GTK3Cairo backend instead. See '
                                 'Matplotlib usage FAQ for more info on '
                                 'backends.'),
                        alternative='GTK3Cairo')
        super(FigureCanvasGTKCairo, self).__init__(*args, **kwargs)

    def _renderer_init(self):
        """Override to use cairo (rather than GDK) renderer"""
        self._renderer = RendererGTKCairo(self.figure.dpi)


# This class has been unused for a while at least.
@cbook.deprecated("2.1")
class FigureManagerGTKCairo(FigureManagerGTK):
    def _get_toolbar(self, canvas):
        # must be inited after the window, drawingArea and figure
        # attrs are set
        if matplotlib.rcParams['toolbar']=='toolbar2':
            toolbar = NavigationToolbar2GTKCairo (canvas, self.window)
        else:
            toolbar = None
        return toolbar


# This class has been unused for a while at least.
@cbook.deprecated("2.1")
class NavigationToolbar2Cairo(NavigationToolbar2GTK):
    def _get_canvas(self, fig):
        return FigureCanvasGTKCairo(fig)


@_BackendGTK.export
class _BackendGTKCairo(_BackendGTK):
    FigureCanvas = FigureCanvasGTKCairo
    FigureManager = FigureManagerGTK
