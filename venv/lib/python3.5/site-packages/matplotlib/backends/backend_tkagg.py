from __future__ import absolute_import, division, print_function

from .. import cbook
from . import tkagg  # Paint image to Tk photo blitter extension.
from .backend_agg import FigureCanvasAgg
from ._backend_tk import (
    _BackendTk, FigureCanvasTk, FigureManagerTk, NavigationToolbar2Tk)


class FigureCanvasTkAgg(FigureCanvasAgg, FigureCanvasTk):
    def draw(self):
        super(FigureCanvasTkAgg, self).draw()
        tkagg.blit(self._tkphoto, self.renderer._renderer, colormode=2)
        self._master.update_idletasks()

    def blit(self, bbox=None):
        tkagg.blit(
            self._tkphoto, self.renderer._renderer, bbox=bbox, colormode=2)
        self._master.update_idletasks()


@cbook.deprecated("2.2")
class FigureManagerTkAgg(FigureManagerTk):
    pass


@cbook.deprecated("2.2")
class NavigationToolbar2TkAgg(NavigationToolbar2Tk):
    pass


@_BackendTk.export
class _BackendTkAgg(_BackendTk):
    FigureCanvas = FigureCanvasTkAgg
