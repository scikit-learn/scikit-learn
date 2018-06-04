from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

from . import backend_cairo, backend_gtk3
from .backend_cairo import cairo, HAS_CAIRO_CFFI
from .backend_gtk3 import _BackendGTK3
from matplotlib.backend_bases import cursors


class RendererGTK3Cairo(backend_cairo.RendererCairo):
    def set_context(self, ctx):
        if HAS_CAIRO_CFFI and not isinstance(ctx, cairo.Context):
            ctx = cairo.Context._from_pointer(
                cairo.ffi.cast(
                    'cairo_t **',
                    id(ctx) + object.__basicsize__)[0],
                incref=True)

        self.gc.ctx = ctx


class FigureCanvasGTK3Cairo(backend_gtk3.FigureCanvasGTK3,
                            backend_cairo.FigureCanvasCairo):

    def _renderer_init(self):
        """Use cairo renderer."""
        self._renderer = RendererGTK3Cairo(self.figure.dpi)

    def _render_figure(self, width, height):
        self._renderer.set_width_height(width, height)
        self.figure.draw(self._renderer)

    def on_draw_event(self, widget, ctx):
        """GtkDrawable draw event."""
        toolbar = self.toolbar
        # if toolbar:
        #     toolbar.set_cursor(cursors.WAIT)
        self._renderer.set_context(ctx)
        allocation = self.get_allocation()
        self._render_figure(allocation.width, allocation.height)
        # if toolbar:
        #     toolbar.set_cursor(toolbar._lastCursor)
        return False  # finish event propagation?


class FigureManagerGTK3Cairo(backend_gtk3.FigureManagerGTK3):
    pass


@_BackendGTK3.export
class _BackendGTK3Cairo(_BackendGTK3):
    FigureCanvas = FigureCanvasGTK3Cairo
    FigureManager = FigureManagerGTK3Cairo
