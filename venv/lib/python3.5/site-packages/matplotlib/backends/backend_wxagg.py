from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import wx

import matplotlib
from matplotlib import cbook
from . import wx_compat as wxc
from .backend_agg import FigureCanvasAgg
from .backend_wx import (
    _BackendWx, _FigureCanvasWxBase, FigureFrameWx,
    NavigationToolbar2Wx as NavigationToolbar2WxAgg)


class FigureFrameWxAgg(FigureFrameWx):
    def get_canvas(self, fig):
        return FigureCanvasWxAgg(self, -1, fig)


class FigureCanvasWxAgg(FigureCanvasAgg, _FigureCanvasWxBase):
    """
    The FigureCanvas contains the figure and does event handling.

    In the wxPython backend, it is derived from wxPanel, and (usually)
    lives inside a frame instantiated by a FigureManagerWx. The parent
    window probably implements a wxSizer to control the displayed
    control size - but we give a hint as to our preferred minimum
    size.
    """

    def draw(self, drawDC=None):
        """
        Render the figure using agg.
        """
        FigureCanvasAgg.draw(self)

        self.bitmap = _convert_agg_to_wx_bitmap(self.get_renderer(), None)
        self._isDrawn = True
        self.gui_repaint(drawDC=drawDC, origin='WXAgg')

    def blit(self, bbox=None):
        """
        Transfer the region of the agg buffer defined by bbox to the display.
        If bbox is None, the entire buffer is transferred.
        """
        if bbox is None:
            self.bitmap = _convert_agg_to_wx_bitmap(self.get_renderer(), None)
            self.gui_repaint()
            return

        l, b, w, h = bbox.bounds
        r = l + w
        t = b + h
        x = int(l)
        y = int(self.bitmap.GetHeight() - t)

        srcBmp = _convert_agg_to_wx_bitmap(self.get_renderer(), None)
        srcDC = wx.MemoryDC()
        srcDC.SelectObject(srcBmp)

        destDC = wx.MemoryDC()
        destDC.SelectObject(self.bitmap)

        destDC.Blit(x, y, int(w), int(h), srcDC, x, y)

        destDC.SelectObject(wx.NullBitmap)
        srcDC.SelectObject(wx.NullBitmap)
        self.gui_repaint()

    filetypes = FigureCanvasAgg.filetypes


@cbook.deprecated("2.2", alternative="NavigationToolbar2WxAgg")
class Toolbar(NavigationToolbar2WxAgg):
    pass


# agg/wxPython image conversion functions (wxPython >= 2.8)

def _convert_agg_to_wx_image(agg, bbox):
    """
    Convert the region of the agg buffer bounded by bbox to a wx.Image.  If
    bbox is None, the entire buffer is converted.

    Note: agg must be a backend_agg.RendererAgg instance.
    """
    if bbox is None:
        # agg => rgb -> image
        image = wxc.EmptyImage(int(agg.width), int(agg.height))
        image.SetData(agg.tostring_rgb())
        return image
    else:
        # agg => rgba buffer -> bitmap => clipped bitmap => image
        return wx.ImageFromBitmap(_WX28_clipped_agg_as_bitmap(agg, bbox))


def _convert_agg_to_wx_bitmap(agg, bbox):
    """
    Convert the region of the agg buffer bounded by bbox to a wx.Bitmap.  If
    bbox is None, the entire buffer is converted.

    Note: agg must be a backend_agg.RendererAgg instance.
    """
    if bbox is None:
        # agg => rgba buffer -> bitmap
        return wxc.BitmapFromBuffer(int(agg.width), int(agg.height),
                                    agg.buffer_rgba())
    else:
        # agg => rgba buffer -> bitmap => clipped bitmap
        return _WX28_clipped_agg_as_bitmap(agg, bbox)


def _WX28_clipped_agg_as_bitmap(agg, bbox):
    """
    Convert the region of a the agg buffer bounded by bbox to a wx.Bitmap.

    Note: agg must be a backend_agg.RendererAgg instance.
    """
    l, b, width, height = bbox.bounds
    r = l + width
    t = b + height

    srcBmp = wxc.BitmapFromBuffer(int(agg.width), int(agg.height),
                                  agg.buffer_rgba())
    srcDC = wx.MemoryDC()
    srcDC.SelectObject(srcBmp)

    destBmp = wxc.EmptyBitmap(int(width), int(height))
    destDC = wx.MemoryDC()
    destDC.SelectObject(destBmp)

    x = int(l)
    y = int(int(agg.height) - t)
    destDC.Blit(0, 0, int(width), int(height), srcDC, x, y)

    srcDC.SelectObject(wx.NullBitmap)
    destDC.SelectObject(wx.NullBitmap)

    return destBmp


@_BackendWx.export
class _BackendWxAgg(_BackendWx):
    FigureCanvas = FigureCanvasWxAgg
    _frame_class = FigureFrameWxAgg
