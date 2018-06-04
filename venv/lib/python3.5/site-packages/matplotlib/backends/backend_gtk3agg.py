from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import numpy as np
import warnings

from . import backend_agg, backend_gtk3
from .backend_cairo import cairo, HAS_CAIRO_CFFI
from .backend_gtk3 import _BackendGTK3
from matplotlib import transforms

if six.PY3 and not HAS_CAIRO_CFFI:
    warnings.warn(
        "The Gtk3Agg backend is known to not work on Python 3.x with pycairo. "
        "Try installing cairocffi.")


class FigureCanvasGTK3Agg(backend_gtk3.FigureCanvasGTK3,
                          backend_agg.FigureCanvasAgg):
    def __init__(self, figure):
        backend_gtk3.FigureCanvasGTK3.__init__(self, figure)
        self._bbox_queue = []

    def _renderer_init(self):
        pass

    def _render_figure(self, width, height):
        backend_agg.FigureCanvasAgg.draw(self)

    def on_draw_event(self, widget, ctx):
        """ GtkDrawable draw event, like expose_event in GTK 2.X
        """
        allocation = self.get_allocation()
        w, h = allocation.width, allocation.height

        if not len(self._bbox_queue):
            self._render_figure(w, h)
            bbox_queue = [transforms.Bbox([[0, 0], [w, h]])]
        else:
            bbox_queue = self._bbox_queue

        if HAS_CAIRO_CFFI and not isinstance(ctx, cairo.Context):
            ctx = cairo.Context._from_pointer(
                cairo.ffi.cast('cairo_t **',
                               id(ctx) + object.__basicsize__)[0],
                incref=True)

        for bbox in bbox_queue:
            area = self.copy_from_bbox(bbox)
            buf = np.fromstring(area.to_string_argb(), dtype='uint8')

            x = int(bbox.x0)
            y = h - int(bbox.y1)
            width = int(bbox.x1) - int(bbox.x0)
            height = int(bbox.y1) - int(bbox.y0)

            if HAS_CAIRO_CFFI:
                image = cairo.ImageSurface.create_for_data(
                    buf.data, cairo.FORMAT_ARGB32, width, height)
            else:
                image = cairo.ImageSurface.create_for_data(
                    buf, cairo.FORMAT_ARGB32, width, height)
            ctx.set_source_surface(image, x, y)
            ctx.paint()

        if len(self._bbox_queue):
            self._bbox_queue = []

        return False

    def blit(self, bbox=None):
        # If bbox is None, blit the entire canvas to gtk. Otherwise
        # blit only the area defined by the bbox.
        if bbox is None:
            bbox = self.figure.bbox

        allocation = self.get_allocation()
        w, h = allocation.width, allocation.height
        x = int(bbox.x0)
        y = h - int(bbox.y1)
        width = int(bbox.x1) - int(bbox.x0)
        height = int(bbox.y1) - int(bbox.y0)

        self._bbox_queue.append(bbox)
        self.queue_draw_area(x, y, width, height)

    def print_png(self, filename, *args, **kwargs):
        # Do this so we can save the resolution of figure in the PNG file
        agg = self.switch_backends(backend_agg.FigureCanvasAgg)
        return agg.print_png(filename, *args, **kwargs)


class FigureManagerGTK3Agg(backend_gtk3.FigureManagerGTK3):
    pass


@_BackendGTK3.export
class _BackendGTK3Cairo(_BackendGTK3):
    FigureCanvas = FigureCanvasGTK3Agg
    FigureManager = FigureManagerGTK3Agg
