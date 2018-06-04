"""
Render to gtk from agg
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import matplotlib
from matplotlib.cbook import warn_deprecated
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.backends.backend_gtk import (
    gtk, _BackendGTK, FigureCanvasGTK, FigureManagerGTK, NavigationToolbar2GTK,
    backend_version, error_msg_gtk, PIXELS_PER_INCH)
from matplotlib.backends._gtkagg import agg_to_gtk_drawable


class NavigationToolbar2GTKAgg(NavigationToolbar2GTK):
    def _get_canvas(self, fig):
        return FigureCanvasGTKAgg(fig)


class FigureManagerGTKAgg(FigureManagerGTK):
    def _get_toolbar(self, canvas):
        # must be inited after the window, drawingArea and figure
        # attrs are set
        if matplotlib.rcParams['toolbar']=='toolbar2':
            toolbar = NavigationToolbar2GTKAgg (canvas, self.window)
        else:
            toolbar = None
        return toolbar


class FigureCanvasGTKAgg(FigureCanvasGTK, FigureCanvasAgg):
    filetypes = FigureCanvasGTK.filetypes.copy()
    filetypes.update(FigureCanvasAgg.filetypes)

    def __init__(self, *args, **kwargs):
        warn_deprecated('2.2',
                        message=('The GTKAgg backend is deprecated. It is '
                                 'untested and will be removed in Matplotlib '
                                 '3.0. Use the GTK3Agg backend instead. See '
                                 'Matplotlib usage FAQ for more info on '
                                 'backends.'),
                        alternative='GTK3Agg')
        super(FigureCanvasGTKAgg, self).__init__(*args, **kwargs)

    def configure_event(self, widget, event=None):

        if widget.window is None:
            return
        try:
            del self.renderer
        except AttributeError:
            pass
        w,h = widget.window.get_size()
        if w==1 or h==1: return # empty fig

        # compute desired figure size in inches
        dpival = self.figure.dpi
        winch = w/dpival
        hinch = h/dpival
        self.figure.set_size_inches(winch, hinch, forward=False)
        self._need_redraw = True
        self.resize_event()
        return True

    def _render_figure(self, pixmap, width, height):
        FigureCanvasAgg.draw(self)

        buf = self.buffer_rgba()
        ren = self.get_renderer()
        w = int(ren.width)
        h = int(ren.height)

        pixbuf = gtk.gdk.pixbuf_new_from_data(
            buf, gtk.gdk.COLORSPACE_RGB,  True, 8, w, h, w*4)
        pixmap.draw_pixbuf(pixmap.new_gc(), pixbuf, 0, 0, 0, 0, w, h,
                           gtk.gdk.RGB_DITHER_NONE, 0, 0)

    def blit(self, bbox=None):
        agg_to_gtk_drawable(self._pixmap, self.renderer._renderer, bbox)
        x, y, w, h = self.allocation
        self.window.draw_drawable(self.style.fg_gc[self.state], self._pixmap,
                                  0, 0, 0, 0, w, h)

    def print_png(self, filename, *args, **kwargs):
        # Do this so we can save the resolution of figure in the PNG file
        agg = self.switch_backends(FigureCanvasAgg)
        return agg.print_png(filename, *args, **kwargs)


@_BackendGTK.export
class _BackendGTKAgg(_BackendGTK):
    FigureCanvas = FigureCanvasGTKAgg
    FigureManager = FigureManagerGTKAgg
