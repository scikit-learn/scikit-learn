#!/usr/bin/env python
"""
A wx API adapter to hide differences between wxPython classic and phoenix.

It is assumed that the user code is selecting what version it wants to use,
here we just ensure that it meets the minimum required by matplotlib.

For an example see embedding_in_wx2.py
"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from distutils.version import StrictVersion, LooseVersion

missingwx = "Matplotlib backend_wx and backend_wxagg require wxPython>=2.9"


try:
    import wx
    backend_version = wx.VERSION_STRING
    is_phoenix = 'phoenix' in wx.PlatformInfo
except ImportError:
    raise ImportError(missingwx)

try:
    wx_version = StrictVersion(wx.VERSION_STRING)
except ValueError:
    wx_version = LooseVersion(wx.VERSION_STRING)

# Ensure we have the correct version imported
if wx_version < str("2.9"):
    raise ImportError(missingwx)

if is_phoenix:
    # define all the wxPython phoenix stuff

    # font styles, families and weight
    fontweights = {
        100: wx.FONTWEIGHT_LIGHT,
        200: wx.FONTWEIGHT_LIGHT,
        300: wx.FONTWEIGHT_LIGHT,
        400: wx.FONTWEIGHT_NORMAL,
        500: wx.FONTWEIGHT_NORMAL,
        600: wx.FONTWEIGHT_NORMAL,
        700: wx.FONTWEIGHT_BOLD,
        800: wx.FONTWEIGHT_BOLD,
        900: wx.FONTWEIGHT_BOLD,
        'ultralight': wx.FONTWEIGHT_LIGHT,
        'light': wx.FONTWEIGHT_LIGHT,
        'normal': wx.FONTWEIGHT_NORMAL,
        'medium': wx.FONTWEIGHT_NORMAL,
        'semibold': wx.FONTWEIGHT_NORMAL,
        'bold': wx.FONTWEIGHT_BOLD,
        'heavy': wx.FONTWEIGHT_BOLD,
        'ultrabold': wx.FONTWEIGHT_BOLD,
        'black': wx.FONTWEIGHT_BOLD
    }
    fontangles = {
        'italic': wx.FONTSTYLE_ITALIC,
        'normal': wx.FONTSTYLE_NORMAL,
        'oblique': wx.FONTSTYLE_SLANT}

    # wxPython allows for portable font styles, choosing them appropriately
    # for the target platform. Map some standard font names to the portable
    # styles
    # QUESTION: Is it be wise to agree standard fontnames across all backends?
    fontnames = {'Sans': wx.FONTFAMILY_SWISS,
                 'Roman': wx.FONTFAMILY_ROMAN,
                 'Script': wx.FONTFAMILY_SCRIPT,
                 'Decorative': wx.FONTFAMILY_DECORATIVE,
                 'Modern': wx.FONTFAMILY_MODERN,
                 'Courier': wx.FONTFAMILY_MODERN,
                 'courier': wx.FONTFAMILY_MODERN}

    dashd_wx = {'solid': wx.PENSTYLE_SOLID,
                'dashed': wx.PENSTYLE_SHORT_DASH,
                'dashdot': wx.PENSTYLE_DOT_DASH,
                'dotted': wx.PENSTYLE_DOT}

    # functions changes
    BitmapFromBuffer = wx.Bitmap.FromBufferRGBA
    EmptyBitmap = wx.Bitmap
    EmptyImage = wx.Image
    Cursor = wx.Cursor
    EventLoop = wx.GUIEventLoop
    NamedColour = wx.Colour
    StockCursor = wx.Cursor

else:
    # define all the wxPython classic stuff

    # font styles, families and weight
    fontweights = {
        100: wx.LIGHT,
        200: wx.LIGHT,
        300: wx.LIGHT,
        400: wx.NORMAL,
        500: wx.NORMAL,
        600: wx.NORMAL,
        700: wx.BOLD,
        800: wx.BOLD,
        900: wx.BOLD,
        'ultralight': wx.LIGHT,
        'light': wx.LIGHT,
        'normal': wx.NORMAL,
        'medium': wx.NORMAL,
        'semibold': wx.NORMAL,
        'bold': wx.BOLD,
        'heavy': wx.BOLD,
        'ultrabold': wx.BOLD,
        'black': wx.BOLD
    }
    fontangles = {
        'italic': wx.ITALIC,
        'normal': wx.NORMAL,
        'oblique': wx.SLANT}

    # wxPython allows for portable font styles, choosing them appropriately
    # for the target platform. Map some standard font names to the portable
    # styles
    # QUESTION: Is it be wise to agree standard fontnames across all backends?
    fontnames = {'Sans': wx.SWISS,
                 'Roman': wx.ROMAN,
                 'Script': wx.SCRIPT,
                 'Decorative': wx.DECORATIVE,
                 'Modern': wx.MODERN,
                 'Courier': wx.MODERN,
                 'courier': wx.MODERN}

    dashd_wx = {'solid': wx.SOLID,
                'dashed': wx.SHORT_DASH,
                'dashdot': wx.DOT_DASH,
                'dotted': wx.DOT}

    # functions changes
    BitmapFromBuffer = wx.BitmapFromBufferRGBA
    EmptyBitmap = wx.EmptyBitmap
    EmptyImage = wx.EmptyImage
    Cursor = wx.StockCursor
    EventLoop = wx.EventLoop
    NamedColour = wx.NamedColour
    StockCursor = wx.StockCursor


# wxPython Classic's DoAddTool has become AddTool in Phoenix. Otherwise
# they are the same, except for early betas and prerelease builds of
# Phoenix. This function provides a shim that does the RightThing based on
# which wxPython is in use.
def _AddTool(parent, wx_ids, text, bmp, tooltip_text):
    if text in ['Pan', 'Zoom']:
        kind = wx.ITEM_CHECK
    else:
        kind = wx.ITEM_NORMAL
    if is_phoenix:
        add_tool = parent.AddTool
    else:
        add_tool = parent.DoAddTool

    if not is_phoenix or wx_version >= str("4.0.0b2"):
        # NOTE: when support for Phoenix prior to 4.0.0b2 is dropped then
        # all that is needed is this clause, and the if and else clause can
        # be removed.
        kwargs = dict(label=text,
                      bitmap=bmp,
                      bmpDisabled=wx.NullBitmap,
                      shortHelp=text,
                      longHelp=tooltip_text,
                      kind=kind)
    else:
        kwargs = dict(label=text,
                      bitmap=bmp,
                      bmpDisabled=wx.NullBitmap,
                      shortHelpString=text,
                      longHelpString=tooltip_text,
                      kind=kind)

    return add_tool(wx_ids[text], **kwargs)
