from __future__ import absolute_import, division, print_function

import sys

import numpy as np

from . import tkagg  # Paint image to Tk photo blitter extension.
from .backend_cairo import cairo, FigureCanvasCairo, RendererCairo
from ._backend_tk import _BackendTk, FigureCanvasTk


class FigureCanvasTkCairo(FigureCanvasCairo, FigureCanvasTk):
    def __init__(self, *args, **kwargs):
        super(FigureCanvasTkCairo, self).__init__(*args, **kwargs)
        self._renderer = RendererCairo(self.figure.dpi)

    def draw(self):
        width = int(self.figure.bbox.width)
        height = int(self.figure.bbox.height)
        surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width, height)
        self._renderer.set_ctx_from_surface(surface)
        self._renderer.set_width_height(width, height)
        self.figure.draw(self._renderer)
        buf = np.reshape(surface.get_data(), (height, width, 4))
        # Convert from ARGB32 to RGBA8888.  Using .take() instead of directly
        # indexing ensures C-contiguity of the result, which is needed by
        # tkagg.
        buf = buf.take(
            [2, 1, 0, 3] if sys.byteorder == "little" else [1, 2, 3, 0],
            axis=2)
        tkagg.blit(self._tkphoto, buf, colormode=2)
        self._master.update_idletasks()


@_BackendTk.export
class _BackendTkCairo(_BackendTk):
    FigureCanvas = FigureCanvasTkCairo
