"""
 A wxPython backend for matplotlib, based (very heavily) on
 backend_template.py and backend_gtk.py

 Author: Jeremy O'Donoghue (jeremy@o-donoghue.com)

 Derived from original copyright work by John Hunter
 (jdhunter@ace.bsd.uchicago.edu)

 Copyright (C) Jeremy O'Donoghue & John Hunter, 2003-4

 License: This work is licensed under a PSF compatible license. A copy
 should be included with this source code.

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import xrange
import six

import sys
import os
import os.path
import math
import weakref
import warnings

import matplotlib
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, GraphicsContextBase,
    NavigationToolbar2, RendererBase, TimerBase, cursors)
from matplotlib.backend_bases import _has_pil

from matplotlib._pylab_helpers import Gcf
from matplotlib.cbook import is_writable_file_like, warn_deprecated
from matplotlib.figure import Figure
from matplotlib.path import Path
from matplotlib.transforms import Affine2D
from matplotlib.widgets import SubplotTool
from matplotlib import cbook, rcParams, backend_tools

from . import wx_compat as wxc
import wx

# Debugging settings here...
# Debug level set here. If the debug level is less than 5, information
# messages (progressively more info for lower value) are printed. In addition,
# traceback is performed, and pdb activated, for all uncaught exceptions in
# this case
_DEBUG = 5
if _DEBUG < 5:
    import traceback
    import pdb
_DEBUG_lvls = {1: 'Low ', 2: 'Med ', 3: 'High', 4: 'Error'}


def DEBUG_MSG(string, lvl=3, o=None):
    if lvl >= _DEBUG:
        cls = o.__class__
        # Jeremy, often times the commented line won't print but the
        # one below does.  I think WX is redefining stderr, damned
        # beast
        # print("%s- %s in %s" % (_DEBUG_lvls[lvl], string, cls),
        #       file=sys.stderr)
        print("%s- %s in %s" % (_DEBUG_lvls[lvl], string, cls))


def debug_on_error(type, value, tb):
    """Code due to Thomas Heller - published in Python Cookbook (O'Reilley)"""
    traceback.print_exception(type, value, tb)
    print()
    pdb.pm()  # jdh uncomment


class fake_stderr(object):
    """
    Wx does strange things with stderr, as it makes the assumption that
    there is probably no console. This redirects stderr to the console, since
    we know that there is one!
    """

    def write(self, msg):
        print("Stderr: %s\n\r" % msg)


# the True dots per inch on the screen; should be display dependent
# see
# http://groups.google.com/groups?q=screen+dpi+x11&hl=en&lr=&ie=UTF-8&oe=UTF-8&safe=off&selm=7077.26e81ad5%40swift.cs.tcd.ie&rnum=5
# for some info about screen dpi
PIXELS_PER_INCH = 75

# Delay time for idle checks
IDLE_DELAY = 5


def error_msg_wx(msg, parent=None):
    """
    Signal an error condition -- in a GUI, popup a error dialog
    """
    dialog = wx.MessageDialog(parent=parent,
                              message=msg,
                              caption='Matplotlib backend_wx error',
                              style=wx.OK | wx.CENTRE)
    dialog.ShowModal()
    dialog.Destroy()
    return None


def raise_msg_to_str(msg):
    """msg is a return arg from a raise.  Join with new lines."""
    if not isinstance(msg, six.string_types):
        msg = '\n'.join(map(str, msg))
    return msg


class TimerWx(TimerBase):
    '''
    Subclass of :class:`backend_bases.TimerBase` that uses WxTimer events.

    Attributes
    ----------
    interval : int
        The time between timer events in milliseconds. Default is 1000 ms.
    single_shot : bool
        Boolean flag indicating whether this timer should operate as single
        shot (run once and then stop). Defaults to False.
    callbacks : list
        Stores list of (func, args) tuples that will be called upon timer
        events. This list can be manipulated directly, or the functions
        `add_callback` and `remove_callback` can be used.

    '''

    def __init__(self, parent, *args, **kwargs):
        TimerBase.__init__(self, *args, **kwargs)

        # Create a new timer and connect the timer event to our handler.
        # For WX, the events have to use a widget for binding.
        self.parent = parent
        self._timer = wx.Timer(self.parent, wx.NewId())
        self.parent.Bind(wx.EVT_TIMER, self._on_timer, self._timer)

     # Unbinding causes Wx to stop for some reason. Disabling for now.
#    def __del__(self):
#        TimerBase.__del__(self)
#        self.parent.Bind(wx.EVT_TIMER, None, self._timer)

    def _timer_start(self):
        self._timer.Start(self._interval, self._single)

    def _timer_stop(self):
        self._timer.Stop()

    def _timer_set_interval(self):
        self._timer_start()

    def _timer_set_single_shot(self):
        self._timer.Start()

    def _on_timer(self, *args):
        TimerBase._on_timer(self)


class RendererWx(RendererBase):
    """
    The renderer handles all the drawing primitives using a graphics
    context instance that controls the colors/styles. It acts as the
    'renderer' instance used by many classes in the hierarchy.
    """
    # In wxPython, drawing is performed on a wxDC instance, which will
    # generally be mapped to the client aread of the window displaying
    # the plot. Under wxPython, the wxDC instance has a wx.Pen which
    # describes the colour and weight of any lines drawn, and a wxBrush
    # which describes the fill colour of any closed polygon.

    fontweights = wxc.fontweights
    fontangles = wxc.fontangles

    # wxPython allows for portable font styles, choosing them appropriately
    # for the target platform. Map some standard font names to the portable
    # styles
    # QUESTION: Is it be wise to agree standard fontnames across all backends?
    fontnames = wxc.fontnames

    def __init__(self, bitmap, dpi):
        """
        Initialise a wxWindows renderer instance.
        """
        warn_deprecated('2.0', message="The WX backend is "
                        "deprecated. It's untested "
                        "and will be removed in Matplotlib 3.0. "
                        "Use the WXAgg backend instead. "
                        "See Matplotlib usage FAQ for more info on backends.",
                        alternative='WXAgg')
        RendererBase.__init__(self)
        DEBUG_MSG("__init__()", 1, self)
        self.width = bitmap.GetWidth()
        self.height = bitmap.GetHeight()
        self.bitmap = bitmap
        self.fontd = {}
        self.dpi = dpi
        self.gc = None

    def flipy(self):
        return True

    def offset_text_height(self):
        return True

    def get_text_width_height_descent(self, s, prop, ismath):
        """
        get the width and height in display coords of the string s
        with FontPropertry prop
        """
        # return 1, 1
        if ismath:
            s = self.strip_math(s)

        if self.gc is None:
            gc = self.new_gc()
        else:
            gc = self.gc
        gfx_ctx = gc.gfx_ctx
        font = self.get_wx_font(s, prop)
        gfx_ctx.SetFont(font, wx.BLACK)
        w, h, descent, leading = gfx_ctx.GetFullTextExtent(s)

        return w, h, descent

    def get_canvas_width_height(self):
        'return the canvas width and height in display coords'
        return self.width, self.height

    def handle_clip_rectangle(self, gc):
        new_bounds = gc.get_clip_rectangle()
        if new_bounds is not None:
            new_bounds = new_bounds.bounds
        gfx_ctx = gc.gfx_ctx
        if gfx_ctx._lastcliprect != new_bounds:
            gfx_ctx._lastcliprect = new_bounds
            if new_bounds is None:
                gfx_ctx.ResetClip()
            else:
                gfx_ctx.Clip(new_bounds[0],
                             self.height - new_bounds[1] - new_bounds[3],
                             new_bounds[2], new_bounds[3])

    @staticmethod
    def convert_path(gfx_ctx, path, transform):
        wxpath = gfx_ctx.CreatePath()
        for points, code in path.iter_segments(transform):
            if code == Path.MOVETO:
                wxpath.MoveToPoint(*points)
            elif code == Path.LINETO:
                wxpath.AddLineToPoint(*points)
            elif code == Path.CURVE3:
                wxpath.AddQuadCurveToPoint(*points)
            elif code == Path.CURVE4:
                wxpath.AddCurveToPoint(*points)
            elif code == Path.CLOSEPOLY:
                wxpath.CloseSubpath()
        return wxpath

    def draw_path(self, gc, path, transform, rgbFace=None):
        gc.select()
        self.handle_clip_rectangle(gc)
        gfx_ctx = gc.gfx_ctx
        transform = transform + \
            Affine2D().scale(1.0, -1.0).translate(0.0, self.height)
        wxpath = self.convert_path(gfx_ctx, path, transform)
        if rgbFace is not None:
            gfx_ctx.SetBrush(wx.Brush(gc.get_wxcolour(rgbFace)))
            gfx_ctx.DrawPath(wxpath)
        else:
            gfx_ctx.StrokePath(wxpath)
        gc.unselect()

    def draw_image(self, gc, x, y, im):
        bbox = gc.get_clip_rectangle()
        if bbox is not None:
            l, b, w, h = bbox.bounds
        else:
            l = 0
            b = 0
            w = self.width
            h = self.height
        rows, cols = im.shape[:2]
        bitmap = wxc.BitmapFromBuffer(cols, rows, im.tostring())
        gc = self.get_gc()
        gc.select()
        gc.gfx_ctx.DrawBitmap(bitmap, int(l), int(self.height - b),
                              int(w), int(-h))
        gc.unselect()

    def draw_text(self, gc, x, y, s, prop, angle, ismath=False, mtext=None):
        if ismath:
            s = self.strip_math(s)
        DEBUG_MSG("draw_text()", 1, self)
        gc.select()
        self.handle_clip_rectangle(gc)
        gfx_ctx = gc.gfx_ctx

        font = self.get_wx_font(s, prop)
        color = gc.get_wxcolour(gc.get_rgb())
        gfx_ctx.SetFont(font, color)

        w, h, d = self.get_text_width_height_descent(s, prop, ismath)
        x = int(x)
        y = int(y - h)

        if angle == 0.0:
            gfx_ctx.DrawText(s, x, y)
        else:
            rads = math.radians(angle)
            xo = h * math.sin(rads)
            yo = h * math.cos(rads)
            gfx_ctx.DrawRotatedText(s, x - xo, y - yo, rads)

        gc.unselect()

    def new_gc(self):
        """
        Return an instance of a GraphicsContextWx, and sets the current gc copy
        """
        DEBUG_MSG('new_gc()', 2, self)
        self.gc = GraphicsContextWx(self.bitmap, self)
        self.gc.select()
        self.gc.unselect()
        return self.gc

    def get_gc(self):
        """
        Fetch the locally cached gc.
        """
        # This is a dirty hack to allow anything with access to a renderer to
        # access the current graphics context
        assert self.gc is not None, "gc must be defined"
        return self.gc

    def get_wx_font(self, s, prop):
        """
        Return a wx font.  Cache instances in a font dictionary for
        efficiency
        """
        DEBUG_MSG("get_wx_font()", 1, self)

        key = hash(prop)
        fontprop = prop
        fontname = fontprop.get_name()

        font = self.fontd.get(key)
        if font is not None:
            return font

        # Allow use of platform independent and dependent font names
        wxFontname = self.fontnames.get(fontname, wx.ROMAN)
        wxFacename = ''  # Empty => wxPython chooses based on wx_fontname

        # Font colour is determined by the active wx.Pen
        # TODO: It may be wise to cache font information
        size = self.points_to_pixels(fontprop.get_size_in_points())

        font = wx.Font(int(size + 0.5),             # Size
                       wxFontname,                # 'Generic' name
                       self.fontangles[fontprop.get_style()],   # Angle
                       self.fontweights[fontprop.get_weight()],  # Weight
                       False,                     # Underline
                       wxFacename)                # Platform font name

        # cache the font and gc and return it
        self.fontd[key] = font

        return font

    def points_to_pixels(self, points):
        """
        convert point measures to pixes using dpi and the pixels per
        inch of the display
        """
        return points * (PIXELS_PER_INCH / 72.0 * self.dpi / 72.0)


class GraphicsContextWx(GraphicsContextBase):
    """
    The graphics context provides the color, line styles, etc...

    This class stores a reference to a wxMemoryDC, and a
    wxGraphicsContext that draws to it.  Creating a wxGraphicsContext
    seems to be fairly heavy, so these objects are cached based on the
    bitmap object that is passed in.

    The base GraphicsContext stores colors as a RGB tuple on the unit
    interval, e.g., (0.5, 0.0, 1.0).  wxPython uses an int interval, but
    since wxPython colour management is rather simple, I have not chosen
    to implement a separate colour manager class.
    """
    _capd = {'butt': wx.CAP_BUTT,
             'projecting': wx.CAP_PROJECTING,
             'round': wx.CAP_ROUND}

    _joind = {'bevel': wx.JOIN_BEVEL,
              'miter': wx.JOIN_MITER,
              'round': wx.JOIN_ROUND}

    _cache = weakref.WeakKeyDictionary()

    def __init__(self, bitmap, renderer):
        GraphicsContextBase.__init__(self)
        # assert self.Ok(), "wxMemoryDC not OK to use"
        DEBUG_MSG("__init__()", 1, self)
        DEBUG_MSG("__init__() 2: %s" % bitmap, 1, self)

        dc, gfx_ctx = self._cache.get(bitmap, (None, None))
        if dc is None:
            dc = wx.MemoryDC()
            dc.SelectObject(bitmap)
            gfx_ctx = wx.GraphicsContext.Create(dc)
            gfx_ctx._lastcliprect = None
            self._cache[bitmap] = dc, gfx_ctx

        self.bitmap = bitmap
        self.dc = dc
        self.gfx_ctx = gfx_ctx
        self._pen = wx.Pen('BLACK', 1, wx.SOLID)
        gfx_ctx.SetPen(self._pen)
        self._style = wx.SOLID
        self.renderer = renderer

    def select(self):
        """
        Select the current bitmap into this wxDC instance
        """

        if sys.platform == 'win32':
            self.dc.SelectObject(self.bitmap)
            self.IsSelected = True

    def unselect(self):
        """
        Select a Null bitmasp into this wxDC instance
        """
        if sys.platform == 'win32':
            self.dc.SelectObject(wx.NullBitmap)
            self.IsSelected = False

    def set_foreground(self, fg, isRGBA=None):
        """
        Set the foreground color.  fg can be a matlab format string, a
        html hex color string, an rgb unit tuple, or a float between 0
        and 1.  In the latter case, grayscale is used.
        """
        # Implementation note: wxPython has a separate concept of pen and
        # brush - the brush fills any outline trace left by the pen.
        # Here we set both to the same colour - if a figure is not to be
        # filled, the renderer will set the brush to be transparent
        # Same goes for text foreground...
        DEBUG_MSG("set_foreground()", 1, self)
        self.select()
        GraphicsContextBase.set_foreground(self, fg, isRGBA)

        self._pen.SetColour(self.get_wxcolour(self.get_rgb()))
        self.gfx_ctx.SetPen(self._pen)
        self.unselect()

    def set_linewidth(self, w):
        """
        Set the line width.
        """
        w = float(w)
        DEBUG_MSG("set_linewidth()", 1, self)
        self.select()
        if w > 0 and w < 1:
            w = 1
        GraphicsContextBase.set_linewidth(self, w)
        lw = int(self.renderer.points_to_pixels(self._linewidth))
        if lw == 0:
            lw = 1
        self._pen.SetWidth(lw)
        self.gfx_ctx.SetPen(self._pen)
        self.unselect()

    def set_capstyle(self, cs):
        """
        Set the capstyle as a string in ('butt', 'round', 'projecting')
        """
        DEBUG_MSG("set_capstyle()", 1, self)
        self.select()
        GraphicsContextBase.set_capstyle(self, cs)
        self._pen.SetCap(GraphicsContextWx._capd[self._capstyle])
        self.gfx_ctx.SetPen(self._pen)
        self.unselect()

    def set_joinstyle(self, js):
        """
        Set the join style to be one of ('miter', 'round', 'bevel')
        """
        DEBUG_MSG("set_joinstyle()", 1, self)
        self.select()
        GraphicsContextBase.set_joinstyle(self, js)
        self._pen.SetJoin(GraphicsContextWx._joind[self._joinstyle])
        self.gfx_ctx.SetPen(self._pen)
        self.unselect()

    @cbook.deprecated("2.1")
    def set_linestyle(self, ls):
        """
        Set the line style to be one of
        """
        DEBUG_MSG("set_linestyle()", 1, self)
        self.select()
        GraphicsContextBase.set_linestyle(self, ls)
        try:
            self._style = wxc.dashd_wx[ls]
        except KeyError:
            self._style = wx.LONG_DASH  # Style not used elsewhere...

        # On MS Windows platform, only line width of 1 allowed for dash lines
        if wx.Platform == '__WXMSW__':
            self.set_linewidth(1)

        self._pen.SetStyle(self._style)
        self.gfx_ctx.SetPen(self._pen)
        self.unselect()

    def get_wxcolour(self, color):
        """return a wx.Colour from RGB format"""
        DEBUG_MSG("get_wx_color()", 1, self)
        if len(color) == 3:
            r, g, b = color
            r *= 255
            g *= 255
            b *= 255
            return wx.Colour(red=int(r), green=int(g), blue=int(b))
        else:
            r, g, b, a = color
            r *= 255
            g *= 255
            b *= 255
            a *= 255
            return wx.Colour(
                red=int(r),
                green=int(g),
                blue=int(b),
                alpha=int(a))


class _FigureCanvasWxBase(FigureCanvasBase, wx.Panel):
    """
    The FigureCanvas contains the figure and does event handling.

    In the wxPython backend, it is derived from wxPanel, and (usually) lives
    inside a frame instantiated by a FigureManagerWx. The parent window
    probably implements a wx.Sizer to control the displayed control size - but
    we give a hint as to our preferred minimum size.
    """

    keyvald = {
        wx.WXK_CONTROL: 'control',
        wx.WXK_SHIFT: 'shift',
        wx.WXK_ALT: 'alt',
        wx.WXK_LEFT: 'left',
        wx.WXK_UP: 'up',
        wx.WXK_RIGHT: 'right',
        wx.WXK_DOWN: 'down',
        wx.WXK_ESCAPE: 'escape',
        wx.WXK_F1: 'f1',
        wx.WXK_F2: 'f2',
        wx.WXK_F3: 'f3',
        wx.WXK_F4: 'f4',
        wx.WXK_F5: 'f5',
        wx.WXK_F6: 'f6',
        wx.WXK_F7: 'f7',
        wx.WXK_F8: 'f8',
        wx.WXK_F9: 'f9',
        wx.WXK_F10: 'f10',
        wx.WXK_F11: 'f11',
        wx.WXK_F12: 'f12',
        wx.WXK_SCROLL: 'scroll_lock',
        wx.WXK_PAUSE: 'break',
        wx.WXK_BACK: 'backspace',
        wx.WXK_RETURN: 'enter',
        wx.WXK_INSERT: 'insert',
        wx.WXK_DELETE: 'delete',
        wx.WXK_HOME: 'home',
        wx.WXK_END: 'end',
        wx.WXK_PAGEUP: 'pageup',
        wx.WXK_PAGEDOWN: 'pagedown',
        wx.WXK_NUMPAD0: '0',
        wx.WXK_NUMPAD1: '1',
        wx.WXK_NUMPAD2: '2',
        wx.WXK_NUMPAD3: '3',
        wx.WXK_NUMPAD4: '4',
        wx.WXK_NUMPAD5: '5',
        wx.WXK_NUMPAD6: '6',
        wx.WXK_NUMPAD7: '7',
        wx.WXK_NUMPAD8: '8',
        wx.WXK_NUMPAD9: '9',
        wx.WXK_NUMPAD_ADD: '+',
        wx.WXK_NUMPAD_SUBTRACT: '-',
        wx.WXK_NUMPAD_MULTIPLY: '*',
        wx.WXK_NUMPAD_DIVIDE: '/',
        wx.WXK_NUMPAD_DECIMAL: 'dec',
        wx.WXK_NUMPAD_ENTER: 'enter',
        wx.WXK_NUMPAD_UP: 'up',
        wx.WXK_NUMPAD_RIGHT: 'right',
        wx.WXK_NUMPAD_DOWN: 'down',
        wx.WXK_NUMPAD_LEFT: 'left',
        wx.WXK_NUMPAD_PAGEUP: 'pageup',
        wx.WXK_NUMPAD_PAGEDOWN: 'pagedown',
        wx.WXK_NUMPAD_HOME: 'home',
        wx.WXK_NUMPAD_END: 'end',
        wx.WXK_NUMPAD_INSERT: 'insert',
        wx.WXK_NUMPAD_DELETE: 'delete',
    }

    def __init__(self, parent, id, figure):
        """
        Initialise a FigureWx instance.

        - Initialise the FigureCanvasBase and wxPanel parents.
        - Set event handlers for:
          EVT_SIZE  (Resize event)
          EVT_PAINT (Paint event)
        """

        FigureCanvasBase.__init__(self, figure)
        # Set preferred window size hint - helps the sizer (if one is
        # connected)
        l, b, w, h = figure.bbox.bounds
        w = int(math.ceil(w))
        h = int(math.ceil(h))

        wx.Panel.__init__(self, parent, id, size=wx.Size(w, h))

        def do_nothing(*args, **kwargs):
            warnings.warn(
                "could not find a setinitialsize function for backend_wx; "
                "please report your wxpython version=%s "
                "to the matplotlib developers list" %
                wxc.backend_version)
            pass

        # try to find the set size func across wx versions
        try:
            getattr(self, 'SetInitialSize')
        except AttributeError:
            self.SetInitialSize = getattr(self, 'SetBestFittingSize',
                                          do_nothing)

        if not hasattr(self, 'IsShownOnScreen'):
            self.IsShownOnScreen = getattr(self, 'IsVisible',
                                           lambda *args: True)

        # Create the drawing bitmap
        self.bitmap = wxc.EmptyBitmap(w, h)
        DEBUG_MSG("__init__() - bitmap w:%d h:%d" % (w, h), 2, self)
        # TODO: Add support for 'point' inspection and plot navigation.
        self._isDrawn = False

        self.Bind(wx.EVT_SIZE, self._onSize)
        self.Bind(wx.EVT_PAINT, self._onPaint)
        self.Bind(wx.EVT_KEY_DOWN, self._onKeyDown)
        self.Bind(wx.EVT_KEY_UP, self._onKeyUp)
        self.Bind(wx.EVT_RIGHT_DOWN, self._onRightButtonDown)
        self.Bind(wx.EVT_RIGHT_DCLICK, self._onRightButtonDClick)
        self.Bind(wx.EVT_RIGHT_UP, self._onRightButtonUp)
        self.Bind(wx.EVT_MOUSEWHEEL, self._onMouseWheel)
        self.Bind(wx.EVT_LEFT_DOWN, self._onLeftButtonDown)
        self.Bind(wx.EVT_LEFT_DCLICK, self._onLeftButtonDClick)
        self.Bind(wx.EVT_LEFT_UP, self._onLeftButtonUp)
        self.Bind(wx.EVT_MOTION, self._onMotion)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._onLeave)
        self.Bind(wx.EVT_ENTER_WINDOW, self._onEnter)
        # Add middle button events
        self.Bind(wx.EVT_MIDDLE_DOWN, self._onMiddleButtonDown)
        self.Bind(wx.EVT_MIDDLE_DCLICK, self._onMiddleButtonDClick)
        self.Bind(wx.EVT_MIDDLE_UP, self._onMiddleButtonUp)

        self.Bind(wx.EVT_MOUSE_CAPTURE_CHANGED, self._onCaptureLost)
        self.Bind(wx.EVT_MOUSE_CAPTURE_LOST, self._onCaptureLost)

        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)  # Reduce flicker.
        self.SetBackgroundColour(wx.WHITE)

        self.macros = {}  # dict from wx id to seq of macros

    def Destroy(self, *args, **kwargs):
        wx.Panel.Destroy(self, *args, **kwargs)

    def Copy_to_Clipboard(self, event=None):
        "copy bitmap of canvas to system clipboard"
        bmp_obj = wx.BitmapDataObject()
        bmp_obj.SetBitmap(self.bitmap)

        if not wx.TheClipboard.IsOpened():
            open_success = wx.TheClipboard.Open()
            if open_success:
                wx.TheClipboard.SetData(bmp_obj)
                wx.TheClipboard.Close()
                wx.TheClipboard.Flush()

    def draw_idle(self):
        """
        Delay rendering until the GUI is idle.
        """
        DEBUG_MSG("draw_idle()", 1, self)
        self._isDrawn = False  # Force redraw
        # Triggering a paint event is all that is needed to defer drawing
        # until later. The platform will send the event when it thinks it is
        # a good time (usually as soon as there are no other events pending).
        self.Refresh(eraseBackground=False)

    def new_timer(self, *args, **kwargs):
        """
        Creates a new backend-specific subclass of
        :class:`backend_bases.Timer`. This is useful for getting periodic
        events through the backend's native event loop. Implemented only
        for backends with GUIs.

        Other Parameters
        ----------------
        interval : scalar
            Timer interval in milliseconds
        callbacks : list
            Sequence of (func, args, kwargs) where ``func(*args, **kwargs)``
            will be executed by the timer every *interval*.

        """
        return TimerWx(self, *args, **kwargs)

    def flush_events(self):
        wx.Yield()

    def start_event_loop(self, timeout=0):
        """
        Start an event loop.  This is used to start a blocking event
        loop so that interactive functions, such as ginput and
        waitforbuttonpress, can wait for events.  This should not be
        confused with the main GUI event loop, which is always running
        and has nothing to do with this.

        This call blocks until a callback function triggers
        stop_event_loop() or *timeout* is reached.  If *timeout* is
        <=0, never timeout.

        Raises RuntimeError if event loop is already running.
        """
        if hasattr(self, '_event_loop'):
            raise RuntimeError("Event loop already running")
        id = wx.NewId()
        timer = wx.Timer(self, id=id)
        if timeout > 0:
            timer.Start(timeout * 1000, oneShot=True)
            self.Bind(wx.EVT_TIMER, self.stop_event_loop, id=id)

        # Event loop handler for start/stop event loop
        self._event_loop = wxc.EventLoop()
        self._event_loop.Run()
        timer.Stop()

    def stop_event_loop(self, event=None):
        """
        Stop an event loop.  This is used to stop a blocking event
        loop so that interactive functions, such as ginput and
        waitforbuttonpress, can wait for events.

        """
        if hasattr(self, '_event_loop'):
            if self._event_loop.IsRunning():
                self._event_loop.Exit()
            del self._event_loop

    def _get_imagesave_wildcards(self):
        'return the wildcard string for the filesave dialog'
        default_filetype = self.get_default_filetype()
        filetypes = self.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(filetypes.items())
        wildcards = []
        extensions = []
        filter_index = 0
        for i, (name, exts) in enumerate(sorted_filetypes):
            ext_list = ';'.join(['*.%s' % ext for ext in exts])
            extensions.append(exts[0])
            wildcard = '%s (%s)|%s' % (name, ext_list, ext_list)
            if default_filetype in exts:
                filter_index = i
            wildcards.append(wildcard)
        wildcards = '|'.join(wildcards)
        return wildcards, extensions, filter_index

    def gui_repaint(self, drawDC=None, origin='WX'):
        """
        Performs update of the displayed image on the GUI canvas, using the
        supplied wx.PaintDC device context.

        The 'WXAgg' backend sets origin accordingly.
        """
        DEBUG_MSG("gui_repaint()", 1, self)
        if self.IsShownOnScreen():
            if not drawDC:
                # not called from OnPaint use a ClientDC
                drawDC = wx.ClientDC(self)

            # following is for 'WX' backend on Windows
            # the bitmap can not be in use by another DC,
            # see GraphicsContextWx._cache
            if wx.Platform == '__WXMSW__' and origin == 'WX':
                img = self.bitmap.ConvertToImage()
                bmp = img.ConvertToBitmap()
                drawDC.DrawBitmap(bmp, 0, 0)
            else:
                drawDC.DrawBitmap(self.bitmap, 0, 0)

    filetypes = FigureCanvasBase.filetypes.copy()
    filetypes['bmp'] = 'Windows bitmap'
    filetypes['jpeg'] = 'JPEG'
    filetypes['jpg'] = 'JPEG'
    filetypes['pcx'] = 'PCX'
    filetypes['png'] = 'Portable Network Graphics'
    filetypes['tif'] = 'Tagged Image Format File'
    filetypes['tiff'] = 'Tagged Image Format File'
    filetypes['xpm'] = 'X pixmap'

    def print_figure(self, filename, *args, **kwargs):
        super(_FigureCanvasWxBase, self).print_figure(
            filename, *args, **kwargs)
        # Restore the current view; this is needed because the artist contains
        # methods rely on particular attributes of the rendered figure for
        # determining things like bounding boxes.
        if self._isDrawn:
            self.draw()

    def _onPaint(self, evt):
        """
        Called when wxPaintEvt is generated
        """

        DEBUG_MSG("_onPaint()", 1, self)
        drawDC = wx.PaintDC(self)
        if not self._isDrawn:
            self.draw(drawDC=drawDC)
        else:
            self.gui_repaint(drawDC=drawDC)
        drawDC.Destroy()

    def _onSize(self, evt):
        """
        Called when wxEventSize is generated.

        In this application we attempt to resize to fit the window, so it
        is better to take the performance hit and redraw the whole window.
        """

        DEBUG_MSG("_onSize()", 2, self)
        sz = self.GetParent().GetSizer()
        if sz:
            si = sz.GetItem(self)
        if sz and si and not si.Proportion and not si.Flag & wx.EXPAND:
            # managed by a sizer, but with a fixed size
            size = self.GetMinSize()
        else:
            # variable size
            size = self.GetClientSize()
        if getattr(self, "_width", None):
            if size == (self._width, self._height):
                # no change in size
                return
        self._width, self._height = size
        # Create a new, correctly sized bitmap
        self.bitmap = wxc.EmptyBitmap(self._width, self._height)

        self._isDrawn = False

        if self._width <= 1 or self._height <= 1:
            return  # Empty figure

        dpival = self.figure.dpi
        winch = self._width / dpival
        hinch = self._height / dpival
        self.figure.set_size_inches(winch, hinch, forward=False)

        # Rendering will happen on the associated paint event
        # so no need to do anything here except to make sure
        # the whole background is repainted.
        self.Refresh(eraseBackground=False)
        FigureCanvasBase.resize_event(self)

    def _get_key(self, evt):

        keyval = evt.KeyCode
        if keyval in self.keyvald:
            key = self.keyvald[keyval]
        elif keyval < 256:
            key = chr(keyval)
            # wx always returns an uppercase, so make it lowercase if the shift
            # key is not depressed (NOTE: this will not handle Caps Lock)
            if not evt.ShiftDown():
                key = key.lower()
        else:
            key = None

        for meth, prefix in (
                [evt.AltDown, 'alt'],
                [evt.ControlDown, 'ctrl'], ):
            if meth():
                key = '{0}+{1}'.format(prefix, key)

        return key

    def _onKeyDown(self, evt):
        """Capture key press."""
        key = self._get_key(evt)
        FigureCanvasBase.key_press_event(self, key, guiEvent=evt)
        if self:
            evt.Skip()

    def _onKeyUp(self, evt):
        """Release key."""
        key = self._get_key(evt)
        FigureCanvasBase.key_release_event(self, key, guiEvent=evt)
        if self:
            evt.Skip()

    def _set_capture(self, capture=True):
        """control wx mouse capture """
        if self.HasCapture():
            self.ReleaseMouse()
        if capture:
            self.CaptureMouse()

    def _onCaptureLost(self, evt):
        """Capture changed or lost"""
        self._set_capture(False)

    def _onRightButtonDown(self, evt):
        """Start measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(True)
        FigureCanvasBase.button_press_event(self, x, y, 3, guiEvent=evt)

    def _onRightButtonDClick(self, evt):
        """Start measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(True)
        FigureCanvasBase.button_press_event(self, x, y, 3,
                                            dblclick=True, guiEvent=evt)

    def _onRightButtonUp(self, evt):
        """End measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(False)
        FigureCanvasBase.button_release_event(self, x, y, 3, guiEvent=evt)

    def _onLeftButtonDown(self, evt):
        """Start measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(True)
        FigureCanvasBase.button_press_event(self, x, y, 1, guiEvent=evt)

    def _onLeftButtonDClick(self, evt):
        """Start measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(True)
        FigureCanvasBase.button_press_event(self, x, y, 1,
                                            dblclick=True, guiEvent=evt)

    def _onLeftButtonUp(self, evt):
        """End measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(False)
        FigureCanvasBase.button_release_event(self, x, y, 1, guiEvent=evt)

    # Add middle button events
    def _onMiddleButtonDown(self, evt):
        """Start measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(True)
        FigureCanvasBase.button_press_event(self, x, y, 2, guiEvent=evt)

    def _onMiddleButtonDClick(self, evt):
        """Start measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(True)
        FigureCanvasBase.button_press_event(self, x, y, 2,
                                            dblclick=True, guiEvent=evt)

    def _onMiddleButtonUp(self, evt):
        """End measuring on an axis."""
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        self._set_capture(False)
        FigureCanvasBase.button_release_event(self, x, y, 2, guiEvent=evt)

    def _onMouseWheel(self, evt):
        """Translate mouse wheel events into matplotlib events"""

        # Determine mouse location
        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()

        # Convert delta/rotation/rate into a floating point step size
        delta = evt.GetWheelDelta()
        rotation = evt.GetWheelRotation()
        rate = evt.GetLinesPerAction()
        step = rate * rotation / delta

        # Done handling event
        evt.Skip()

        # Mac is giving two events for every wheel event
        # Need to skip every second one
        if wx.Platform == '__WXMAC__':
            if not hasattr(self, '_skipwheelevent'):
                self._skipwheelevent = True
            elif self._skipwheelevent:
                self._skipwheelevent = False
                return  # Return without processing event
            else:
                self._skipwheelevent = True

        # Convert to mpl event
        FigureCanvasBase.scroll_event(self, x, y, step, guiEvent=evt)

    def _onMotion(self, evt):
        """Start measuring on an axis."""

        x = evt.GetX()
        y = self.figure.bbox.height - evt.GetY()
        evt.Skip()
        FigureCanvasBase.motion_notify_event(self, x, y, guiEvent=evt)

    def _onLeave(self, evt):
        """Mouse has left the window."""

        evt.Skip()
        FigureCanvasBase.leave_notify_event(self, guiEvent=evt)

    def _onEnter(self, evt):
        """Mouse has entered the window."""
        FigureCanvasBase.enter_notify_event(self, guiEvent=evt)


class FigureCanvasWx(_FigureCanvasWxBase):
    # Rendering to a Wx canvas using the deprecated Wx renderer.

    def draw(self, drawDC=None):
        """
        Render the figure using RendererWx instance renderer, or using a
        previously defined renderer if none is specified.
        """
        DEBUG_MSG("draw()", 1, self)
        self.renderer = RendererWx(self.bitmap, self.figure.dpi)
        self.figure.draw(self.renderer)
        self._isDrawn = True
        self.gui_repaint(drawDC=drawDC)

    def print_bmp(self, filename, *args, **kwargs):
        return self._print_image(filename, wx.BITMAP_TYPE_BMP, *args, **kwargs)

    if not _has_pil:
        def print_jpeg(self, filename, *args, **kwargs):
            return self._print_image(filename, wx.BITMAP_TYPE_JPEG,
                                     *args, **kwargs)
        print_jpg = print_jpeg

    def print_pcx(self, filename, *args, **kwargs):
        return self._print_image(filename, wx.BITMAP_TYPE_PCX, *args, **kwargs)

    def print_png(self, filename, *args, **kwargs):
        return self._print_image(filename, wx.BITMAP_TYPE_PNG, *args, **kwargs)

    if not _has_pil:
        def print_tiff(self, filename, *args, **kwargs):
            return self._print_image(filename, wx.BITMAP_TYPE_TIF,
                                     *args, **kwargs)
        print_tif = print_tiff

    def print_xpm(self, filename, *args, **kwargs):
        return self._print_image(filename, wx.BITMAP_TYPE_XPM, *args, **kwargs)

    def _print_image(self, filename, filetype, *args, **kwargs):
        origBitmap = self.bitmap

        l, b, width, height = self.figure.bbox.bounds
        width = int(math.ceil(width))
        height = int(math.ceil(height))

        self.bitmap = wxc.EmptyBitmap(width, height)

        renderer = RendererWx(self.bitmap, self.figure.dpi)

        gc = renderer.new_gc()

        self.figure.draw(renderer)

        # image is the object that we call SaveFile on.
        image = self.bitmap
        # set the JPEG quality appropriately.  Unfortunately, it is only
        # possible to set the quality on a wx.Image object.  So if we
        # are saving a JPEG, convert the wx.Bitmap to a wx.Image,
        # and set the quality.
        if filetype == wx.BITMAP_TYPE_JPEG:
            jpeg_quality = kwargs.get('quality',
                                      rcParams['savefig.jpeg_quality'])
            image = self.bitmap.ConvertToImage()
            image.SetOption(wx.IMAGE_OPTION_QUALITY, str(jpeg_quality))

        # Now that we have rendered into the bitmap, save it
        # to the appropriate file type and clean up
        if isinstance(filename, six.string_types):
            if not image.SaveFile(filename, filetype):
                DEBUG_MSG('print_figure() file save error', 4, self)
                raise RuntimeError(
                    'Could not save figure to %s\n' %
                    (filename))
        elif is_writable_file_like(filename):
            if not isinstance(image, wx.Image):
                image = image.ConvertToImage()
            if not image.SaveStream(filename, filetype):
                DEBUG_MSG('print_figure() file save error', 4, self)
                raise RuntimeError(
                    'Could not save figure to %s\n' %
                    (filename))

        # Restore everything to normal
        self.bitmap = origBitmap

        # Note: draw is required here since bits of state about the
        # last renderer are strewn about the artist draw methods.  Do
        # not remove the draw without first verifying that these have
        # been cleaned up.  The artist contains() methods will fail
        # otherwise.
        if self._isDrawn:
            self.draw()
        self.Refresh()


########################################################################
#
# The following functions and classes are for pylab compatibility
# mode (matplotlib.pylab) and implement figure managers, etc...
#
########################################################################


class FigureFrameWx(wx.Frame):
    def __init__(self, num, fig):
        # On non-Windows platform, explicitly set the position - fix
        # positioning bug on some Linux platforms
        if wx.Platform == '__WXMSW__':
            pos = wx.DefaultPosition
        else:
            pos = wx.Point(20, 20)
        l, b, w, h = fig.bbox.bounds
        wx.Frame.__init__(self, parent=None, id=-1, pos=pos,
                          title="Figure %d" % num)
        # Frame will be sized later by the Fit method
        DEBUG_MSG("__init__()", 1, self)
        self.num = num

        statbar = StatusBarWx(self)
        self.SetStatusBar(statbar)
        self.canvas = self.get_canvas(fig)
        self.canvas.SetInitialSize(wx.Size(fig.bbox.width, fig.bbox.height))
        self.canvas.SetFocus()
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)
        # By adding toolbar in sizer, we are able to put it at the bottom
        # of the frame - so appearance is closer to GTK version

        self.toolbar = self._get_toolbar(statbar)

        if self.toolbar is not None:
            self.toolbar.Realize()
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            if wxc.is_phoenix:
                tw, th = self.toolbar.GetSize()
                fw, fh = self.canvas.GetSize()
            else:
                tw, th = self.toolbar.GetSizeTuple()
                fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            self.toolbar.SetSize(wx.Size(fw, th))
            self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

        self.canvas.SetMinSize((2, 2))

        # give the window a matplotlib icon rather than the stock one.
        # This is not currently working on Linux and is untested elsewhere.
        # icon_path = os.path.join(matplotlib.rcParams['datapath'],
        #                         'images', 'matplotlib.png')
        # icon = wx.IconFromBitmap(wx.Bitmap(icon_path))
        #  for xpm type icons try:
        # icon = wx.Icon(icon_path, wx.BITMAP_TYPE_XPM)
        #  self.SetIcon(icon)

        self.figmgr = FigureManagerWx(self.canvas, num, self)

        self.Bind(wx.EVT_CLOSE, self._onClose)

    def _get_toolbar(self, statbar):
        if rcParams['toolbar'] == 'toolbar2':
            toolbar = NavigationToolbar2Wx(self.canvas)
            toolbar.set_status_bar(statbar)
        else:
            toolbar = None
        return toolbar

    def get_canvas(self, fig):
        return FigureCanvasWx(self, -1, fig)

    def get_figure_manager(self):
        DEBUG_MSG("get_figure_manager()", 1, self)
        return self.figmgr

    def _onClose(self, evt):
        DEBUG_MSG("onClose()", 1, self)
        self.canvas.close_event()
        self.canvas.stop_event_loop()
        Gcf.destroy(self.num)
        # self.Destroy()

    def GetToolBar(self):
        """Override wxFrame::GetToolBar as we don't have managed toolbar"""
        return self.toolbar

    def Destroy(self, *args, **kwargs):
        try:
            self.canvas.mpl_disconnect(self.toolbar._idDrag)
            # Rationale for line above: see issue 2941338.
        except AttributeError:
            pass  # classic toolbar lacks the attribute
        if not self.IsBeingDeleted():
            wx.Frame.Destroy(self, *args, **kwargs)
            if self.toolbar is not None:
                self.toolbar.Destroy()
            wxapp = wx.GetApp()
            if wxapp:
                wxapp.Yield()
        return True


class FigureManagerWx(FigureManagerBase):
    """
    This class contains the FigureCanvas and GUI frame

    It is instantiated by GcfWx whenever a new figure is created. GcfWx is
    responsible for managing multiple instances of FigureManagerWx.

    Attributes
    ----------
    canvas : `FigureCanvas`
        a FigureCanvasWx(wx.Panel) instance
    window : wxFrame
        a wxFrame instance - wxpython.org/Phoenix/docs/html/Frame.html

    """

    def __init__(self, canvas, num, frame):
        DEBUG_MSG("__init__()", 1, self)
        FigureManagerBase.__init__(self, canvas, num)
        self.frame = frame
        self.window = frame

        self.tb = frame.GetToolBar()
        self.toolbar = self.tb  # consistent with other backends

        def notify_axes_change(fig):
            'this will be called whenever the current axes is changed'
            if self.tb is not None:
                self.tb.update()
        self.canvas.figure.add_axobserver(notify_axes_change)

    def show(self):
        self.frame.Show()
        self.canvas.draw()

    def destroy(self, *args):
        DEBUG_MSG("destroy()", 1, self)
        self.frame.Destroy()
        wxapp = wx.GetApp()
        if wxapp:
            wxapp.Yield()

    def get_window_title(self):
        return self.window.GetTitle()

    def set_window_title(self, title):
        self.window.SetTitle(title)

    def resize(self, width, height):
        'Set the canvas size in pixels'
        self.canvas.SetInitialSize(wx.Size(width, height))
        self.window.GetSizer().Fit(self.window)

# Identifiers for toolbar controls - images_wx contains bitmaps for the images
# used in the controls. wxWindows does not provide any stock images, so I've
# 'stolen' those from GTK2, and transformed them into the appropriate format.
# import images_wx


_NTB_AXISMENU = wx.NewId()
_NTB_AXISMENU_BUTTON = wx.NewId()
_NTB_X_PAN_LEFT = wx.NewId()
_NTB_X_PAN_RIGHT = wx.NewId()
_NTB_X_ZOOMIN = wx.NewId()
_NTB_X_ZOOMOUT = wx.NewId()
_NTB_Y_PAN_UP = wx.NewId()
_NTB_Y_PAN_DOWN = wx.NewId()
_NTB_Y_ZOOMIN = wx.NewId()
_NTB_Y_ZOOMOUT = wx.NewId()
# _NTB_SUBPLOT            =wx.NewId()
_NTB_SAVE = wx.NewId()
_NTB_CLOSE = wx.NewId()


def _load_bitmap(filename):
    """
    Load a bitmap file from the backends/images subdirectory in which the
    matplotlib library is installed. The filename parameter should not
    contain any path information as this is determined automatically.

    Returns a wx.Bitmap object
    """

    basedir = os.path.join(rcParams['datapath'], 'images')

    bmpFilename = os.path.normpath(os.path.join(basedir, filename))
    if not os.path.exists(bmpFilename):
        raise IOError('Could not find bitmap file "%s"; dying' % bmpFilename)

    bmp = wx.Bitmap(bmpFilename)
    return bmp


class MenuButtonWx(wx.Button):
    """
    wxPython does not permit a menu to be incorporated directly into a toolbar.
    This class simulates the effect by associating a pop-up menu with a button
    in the toolbar, and managing this as though it were a menu.
    """

    def __init__(self, parent):

        wx.Button.__init__(self, parent, _NTB_AXISMENU_BUTTON, "Axes:        ",
                           style=wx.BU_EXACTFIT)
        self._toolbar = parent
        self._menu = wx.Menu()
        self._axisId = []
        # First two menu items never change...
        self._allId = wx.NewId()
        self._invertId = wx.NewId()
        self._menu.Append(self._allId, "All", "Select all axes", False)
        self._menu.Append(self._invertId, "Invert", "Invert axes selected",
                          False)
        self._menu.AppendSeparator()

        self.Bind(wx.EVT_BUTTON, self._onMenuButton, id=_NTB_AXISMENU_BUTTON)
        self.Bind(wx.EVT_MENU, self._handleSelectAllAxes, id=self._allId)
        self.Bind(wx.EVT_MENU, self._handleInvertAxesSelected,
                  id=self._invertId)

    def Destroy(self):
        self._menu.Destroy()
        self.Destroy()

    def _onMenuButton(self, evt):
        """Handle menu button pressed."""
        if wxc.is_phoenix:
            x, y = self.GetPosition()
            w, h = self.GetSize()
        else:
            x, y = self.GetPositionTuple()
            w, h = self.GetSizeTuple()
        self.PopupMenuXY(self._menu, x, y + h - 4)
        # When menu returned, indicate selection in button
        evt.Skip()

    def _handleSelectAllAxes(self, evt):
        """Called when the 'select all axes' menu item is selected."""
        if len(self._axisId) == 0:
            return
        for i in range(len(self._axisId)):
            self._menu.Check(self._axisId[i], True)
        self._toolbar.set_active(self.getActiveAxes())
        evt.Skip()

    def _handleInvertAxesSelected(self, evt):
        """Called when the invert all menu item is selected"""
        if len(self._axisId) == 0:
            return
        for i in range(len(self._axisId)):
            if self._menu.IsChecked(self._axisId[i]):
                self._menu.Check(self._axisId[i], False)
            else:
                self._menu.Check(self._axisId[i], True)
        self._toolbar.set_active(self.getActiveAxes())
        evt.Skip()

    def _onMenuItemSelected(self, evt):
        """Called whenever one of the specific axis menu items is selected"""
        current = self._menu.IsChecked(evt.GetId())
        if current:
            new = False
        else:
            new = True
        self._menu.Check(evt.GetId(), new)
        # Lines above would be deleted based on svn tracker ID 2841525;
        # not clear whether this matters or not.
        self._toolbar.set_active(self.getActiveAxes())
        evt.Skip()

    def updateAxes(self, maxAxis):
        """Ensures that there are entries for max_axis axes in the menu
        (selected by default)."""
        if maxAxis > len(self._axisId):
            for i in range(len(self._axisId) + 1, maxAxis + 1, 1):
                menuId = wx.NewId()
                self._axisId.append(menuId)
                self._menu.Append(menuId, "Axis %d" % i,
                                  "Select axis %d" % i,
                                  True)
                self._menu.Check(menuId, True)
                self.Bind(wx.EVT_MENU, self._onMenuItemSelected, id=menuId)
        elif maxAxis < len(self._axisId):
            for menuId in self._axisId[maxAxis:]:
                self._menu.Delete(menuId)
            self._axisId = self._axisId[:maxAxis]
        self._toolbar.set_active(list(xrange(maxAxis)))

    def getActiveAxes(self):
        """Return a list of the selected axes."""
        active = []
        for i in range(len(self._axisId)):
            if self._menu.IsChecked(self._axisId[i]):
                active.append(i)
        return active

    def updateButtonText(self, lst):
        """Update the list of selected axes in the menu button."""
        self.SetLabel(
            'Axes: ' + ','.join('%d' % (e + 1) for e in lst))


cursord = {
    cursors.MOVE: wx.CURSOR_HAND,
    cursors.HAND: wx.CURSOR_HAND,
    cursors.POINTER: wx.CURSOR_ARROW,
    cursors.SELECT_REGION: wx.CURSOR_CROSS,
    cursors.WAIT: wx.CURSOR_WAIT,
}


@cbook.deprecated("2.2")
class SubplotToolWX(wx.Frame):
    def __init__(self, targetfig):
        wx.Frame.__init__(self, None, -1, "Configure subplots")

        toolfig = Figure((6, 3))
        canvas = FigureCanvasWx(self, -1, toolfig)

        # Create a figure manager to manage things
        figmgr = FigureManager(canvas, 1, self)

        # Now put all into a sizer
        sizer = wx.BoxSizer(wx.VERTICAL)
        # This way of adding to sizer allows resizing
        sizer.Add(canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(sizer)
        self.Fit()
        tool = SubplotTool(targetfig, toolfig)


class NavigationToolbar2Wx(NavigationToolbar2, wx.ToolBar):
    def __init__(self, canvas):
        wx.ToolBar.__init__(self, canvas.GetParent(), -1)
        NavigationToolbar2.__init__(self, canvas)
        self.canvas = canvas
        self._idle = True
        self.statbar = None
        self.prevZoomRect = None
        # for now, use alternate zoom-rectangle drawing on all
        # Macs. N.B. In future versions of wx it may be possible to
        # detect Retina displays with window.GetContentScaleFactor()
        # and/or dc.GetContentScaleFactor()
        self.retinaFix = 'wxMac' in wx.PlatformInfo

    def get_canvas(self, frame, fig):
        return type(self.canvas)(frame, -1, fig)

    def _init_toolbar(self):
        DEBUG_MSG("_init_toolbar", 1, self)

        self._parent = self.canvas.GetParent()

        self.wx_ids = {}
        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.AddSeparator()
                continue
            self.wx_ids[text] = wx.NewId()
            wxc._AddTool(self, self.wx_ids, text,
                         _load_bitmap(image_file + '.png'),
                         tooltip_text)

            self.Bind(wx.EVT_TOOL, getattr(self, callback),
                      id=self.wx_ids[text])

        self.Realize()

    def zoom(self, *args):
        self.ToggleTool(self.wx_ids['Pan'], False)
        NavigationToolbar2.zoom(self, *args)

    def pan(self, *args):
        self.ToggleTool(self.wx_ids['Zoom'], False)
        NavigationToolbar2.pan(self, *args)

    def configure_subplots(self, evt):
        frame = wx.Frame(None, -1, "Configure subplots")

        toolfig = Figure((6, 3))
        canvas = self.get_canvas(frame, toolfig)

        # Create a figure manager to manage things
        figmgr = FigureManager(canvas, 1, frame)

        # Now put all into a sizer
        sizer = wx.BoxSizer(wx.VERTICAL)
        # This way of adding to sizer allows resizing
        sizer.Add(canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        frame.SetSizer(sizer)
        frame.Fit()
        tool = SubplotTool(self.canvas.figure, toolfig)
        frame.Show()

    def save_figure(self, *args):
        # Fetch the required filename and file type.
        filetypes, exts, filter_index = self.canvas._get_imagesave_wildcards()
        default_file = self.canvas.get_default_filename()
        dlg = wx.FileDialog(self._parent, "Save to file", "", default_file,
                            filetypes,
                            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        dlg.SetFilterIndex(filter_index)
        if dlg.ShowModal() == wx.ID_OK:
            dirname = dlg.GetDirectory()
            filename = dlg.GetFilename()
            DEBUG_MSG(
                'Save file dir:%s name:%s' %
                (dirname, filename), 3, self)
            format = exts[dlg.GetFilterIndex()]
            basename, ext = os.path.splitext(filename)
            if ext.startswith('.'):
                ext = ext[1:]
            if ext in ('svg', 'pdf', 'ps', 'eps', 'png') and format != ext:
                # looks like they forgot to set the image type drop
                # down, going with the extension.
                warnings.warn(
                    'extension %s did not match the selected '
                    'image type %s; going with %s' %
                    (ext, format, ext), stacklevel=0)
                format = ext
            try:
                self.canvas.figure.savefig(
                    os.path.join(dirname, filename), format=format)
            except Exception as e:
                error_msg_wx(str(e))

    def set_cursor(self, cursor):
        cursor = wxc.Cursor(cursord[cursor])
        self.canvas.SetCursor(cursor)
        self.canvas.Update()

    @cbook.deprecated("2.1", alternative="canvas.draw_idle")
    def dynamic_update(self):
        d = self._idle
        self._idle = False
        if d:
            self.canvas.draw()
            self._idle = True

    def press(self, event):
        if self._active == 'ZOOM':
            if not self.retinaFix:
                self.wxoverlay = wx.Overlay()
            else:
                if event.inaxes is not None:
                    self.savedRetinaImage = self.canvas.copy_from_bbox(
                        event.inaxes.bbox)
                    self.zoomStartX = event.xdata
                    self.zoomStartY = event.ydata
                    self.zoomAxes = event.inaxes

    def release(self, event):
        if self._active == 'ZOOM':
            # When the mouse is released we reset the overlay and it
            # restores the former content to the window.
            if not self.retinaFix:
                self.wxoverlay.Reset()
                del self.wxoverlay
            else:
                del self.savedRetinaImage
                if self.prevZoomRect:
                    self.prevZoomRect.pop(0).remove()
                    self.prevZoomRect = None
                if self.zoomAxes:
                    self.zoomAxes = None

    def draw_rubberband(self, event, x0, y0, x1, y1):
        if self.retinaFix:  # On Macs, use the following code
            # wx.DCOverlay does not work properly on Retina displays.
            rubberBandColor = '#C0C0FF'
            if self.prevZoomRect:
                self.prevZoomRect.pop(0).remove()
            self.canvas.restore_region(self.savedRetinaImage)
            X0, X1 = self.zoomStartX, event.xdata
            Y0, Y1 = self.zoomStartY, event.ydata
            lineX = (X0, X0, X1, X1, X0)
            lineY = (Y0, Y1, Y1, Y0, Y0)
            self.prevZoomRect = self.zoomAxes.plot(
                lineX, lineY, '-', color=rubberBandColor)
            self.zoomAxes.draw_artist(self.prevZoomRect[0])
            self.canvas.blit(self.zoomAxes.bbox)
            return

        # Use an Overlay to draw a rubberband-like bounding box.

        dc = wx.ClientDC(self.canvas)
        odc = wx.DCOverlay(self.wxoverlay, dc)
        odc.Clear()

        # Mac's DC is already the same as a GCDC, and it causes
        # problems with the overlay if we try to use an actual
        # wx.GCDC so don't try it.
        if 'wxMac' not in wx.PlatformInfo:
            dc = wx.GCDC(dc)

        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0

        if y1 < y0:
            y0, y1 = y1, y0
        if x1 < y0:
            x0, x1 = x1, x0

        w = x1 - x0
        h = y1 - y0
        rect = wx.Rect(x0, y0, w, h)

        rubberBandColor = '#C0C0FF'  # or load from config?

        # Set a pen for the border
        color = wxc.NamedColour(rubberBandColor)
        dc.SetPen(wx.Pen(color, 1))

        # use the same color, plus alpha for the brush
        r, g, b, a = color.Get(True)
        color.Set(r, g, b, 0x60)
        dc.SetBrush(wx.Brush(color))
        if wxc.is_phoenix:
            dc.DrawRectangle(rect)
        else:
            dc.DrawRectangleRect(rect)

    def set_status_bar(self, statbar):
        self.statbar = statbar

    def set_message(self, s):
        if self.statbar is not None:
            self.statbar.set_function(s)

    def set_history_buttons(self):
        can_backward = self._nav_stack._pos > 0
        can_forward = self._nav_stack._pos < len(self._nav_stack._elements) - 1
        self.EnableTool(self.wx_ids['Back'], can_backward)
        self.EnableTool(self.wx_ids['Forward'], can_forward)


@cbook.deprecated("2.2", alternative="NavigationToolbar2Wx")
class Toolbar(NavigationToolbar2Wx):
    pass


class StatusBarWx(wx.StatusBar):
    """
    A status bar is added to _FigureFrame to allow measurements and the
    previously selected scroll function to be displayed as a user
    convenience.
    """

    def __init__(self, parent):
        wx.StatusBar.__init__(self, parent, -1)
        self.SetFieldsCount(2)
        self.SetStatusText("None", 1)
        # self.SetStatusText("Measurement: None", 2)
        # self.Reposition()

    def set_function(self, string):
        self.SetStatusText("%s" % string, 1)

    # def set_measurement(self, string):
    #    self.SetStatusText("Measurement: %s" % string, 2)


# tools for matplotlib.backend_managers.ToolManager:
# for now only SaveFigure, SetCursor and Rubberband are implemented
# once a ToolbarWx is implemented, also FigureManagerWx needs to be
# modified, similar to pull request #9934

class SaveFigureWx(backend_tools.SaveFigureBase):
    def trigger(self, *args):
        # Fetch the required filename and file type.
        filetypes, exts, filter_index = self.canvas._get_imagesave_wildcards()
        default_dir = os.path.expanduser(
            matplotlib.rcParams['savefig.directory'])
        default_file = self.canvas.get_default_filename()
        dlg = wx.FileDialog(self.canvas.GetTopLevelParent(), "Save to file",
                            default_dir, default_file, filetypes,
                            wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        dlg.SetFilterIndex(filter_index)
        if dlg.ShowModal() != wx.ID_OK:
            return

        dirname = dlg.GetDirectory()
        filename = dlg.GetFilename()
        DEBUG_MSG('Save file dir:%s name:%s' % (dirname, filename), 3, self)
        format = exts[dlg.GetFilterIndex()]
        basename, ext = os.path.splitext(filename)
        if ext.startswith('.'):
            ext = ext[1:]
        if ext in ('svg', 'pdf', 'ps', 'eps', 'png') and format != ext:
            # looks like they forgot to set the image type drop
            # down, going with the extension.
            warnings.warn(
                'extension %s did not match the selected '
                'image type %s; going with %s' %
                (ext, format, ext), stacklevel=0)
            format = ext
        if default_dir != "":
            matplotlib.rcParams['savefig.directory'] = dirname
        try:
            self.canvas.figure.savefig(
                os.path.join(dirname, filename), format=format)
        except Exception as e:
            error_msg_wx(str(e))


class SetCursorWx(backend_tools.SetCursorBase):
    def set_cursor(self, cursor):
        cursor = wxc.Cursor(cursord[cursor])
        self.canvas.SetCursor(cursor)
        self.canvas.Update()


if 'wxMac' not in wx.PlatformInfo:
    # on most platforms, use overlay
    class RubberbandWx(backend_tools.RubberbandBase):
        def __init__(self, *args, **kwargs):
            backend_tools.RubberbandBase.__init__(self, *args, **kwargs)
            self.wxoverlay = None

        def draw_rubberband(self, x0, y0, x1, y1):
            # Use an Overlay to draw a rubberband-like bounding box.
            if self.wxoverlay is None:
                self.wxoverlay = wx.Overlay()
            dc = wx.ClientDC(self.canvas)
            odc = wx.DCOverlay(self.wxoverlay, dc)
            odc.Clear()

            dc = wx.GCDC(dc)

            height = self.canvas.figure.bbox.height
            y1 = height - y1
            y0 = height - y0

            if y1 < y0:
                y0, y1 = y1, y0
            if x1 < y0:
                x0, x1 = x1, x0

            w = x1 - x0
            h = y1 - y0
            rect = wx.Rect(x0, y0, w, h)

            rubberBandColor = '#C0C0FF'  # or load from config?

            # Set a pen for the border
            color = wxc.NamedColour(rubberBandColor)
            dc.SetPen(wx.Pen(color, 1))

            # use the same color, plus alpha for the brush
            r, g, b, a = color.Get(True)
            color.Set(r, g, b, 0x60)
            dc.SetBrush(wx.Brush(color))
            if wxc.is_phoenix:
                dc.DrawRectangle(rect)
            else:
                dc.DrawRectangleRect(rect)

        def remove_rubberband(self):
            if self.wxoverlay is None:
                return
            self.wxoverlay.Reset()
            self.wxoverlay = None

else:
    # on Mac OS retina displays DCOverlay does not work
    # and dc.SetLogicalFunction does not have an effect on any display
    # the workaround is to blit the full image for remove_rubberband
    class RubberbandWx(backend_tools.RubberbandBase):
        def __init__(self, *args, **kwargs):
            backend_tools.RubberbandBase.__init__(self, *args, **kwargs)
            self._rect = None

        def draw_rubberband(self, x0, y0, x1, y1):
            dc = wx.ClientDC(self.canvas)
            # this would be required if the Canvas is a ScrolledWindow,
            # which is not the case for now
            # self.PrepareDC(dc)

            # delete old rubberband
            if self._rect:
                self.remove_rubberband(dc)

            # draw new rubberband
            dc.SetPen(wx.Pen(wx.BLACK, 1, wx.SOLID))
            dc.SetBrush(wx.TRANSPARENT_BRUSH)
            self._rect = (x0, self.canvas._height-y0, x1-x0, -y1+y0)
            if wxc.is_phoenix:
                dc.DrawRectangle(self._rect)
            else:
                dc.DrawRectangleRect(self._rect)

        def remove_rubberband(self, dc=None):
            if not self._rect:
                return
            if self.canvas.bitmap:
                if dc is None:
                    dc = wx.ClientDC(self.canvas)
                dc.DrawBitmap(self.canvas.bitmap, 0, 0)
                #  for testing the method on Windows, use this code instead:
                # img = self.canvas.bitmap.ConvertToImage()
                # bmp = img.ConvertToBitmap()
                # dc.DrawBitmap(bmp, 0, 0)
            self._rect = None


backend_tools.ToolSaveFigure = SaveFigureWx
backend_tools.ToolSetCursor = SetCursorWx
backend_tools.ToolRubberband = RubberbandWx


# < Additions for printing support: Matt Newville

class PrintoutWx(wx.Printout):
    """
    Simple wrapper around wx Printout class -- all the real work
    here is scaling the matplotlib canvas bitmap to the current
    printer's definition.
    """

    def __init__(self, canvas, width=5.5, margin=0.5, title='matplotlib'):
        wx.Printout.__init__(self, title=title)
        self.canvas = canvas
        # width, in inches of output figure (approximate)
        self.width = width
        self.margin = margin

    def HasPage(self, page):
        # current only supports 1 page print
        return page == 1

    def GetPageInfo(self):
        return (1, 1, 1, 1)

    def OnPrintPage(self, page):
        self.canvas.draw()

        dc = self.GetDC()
        (ppw, pph) = self.GetPPIPrinter()      # printer's pixels per in
        (pgw, pgh) = self.GetPageSizePixels()  # page size in pixels
        (dcw, dch) = dc.GetSize()
        if wxc.is_phoenix:
            (grw, grh) = self.canvas.GetSize()
        else:
            (grw, grh) = self.canvas.GetSizeTuple()

        # save current figure dpi resolution and bg color,
        # so that we can temporarily set them to the dpi of
        # the printer, and the bg color to white
        bgcolor = self.canvas.figure.get_facecolor()
        fig_dpi = self.canvas.figure.dpi

        # draw the bitmap, scaled appropriately
        vscale = float(ppw) / fig_dpi

        # set figure resolution,bg color for printer
        self.canvas.figure.dpi = ppw
        self.canvas.figure.set_facecolor('#FFFFFF')

        renderer = RendererWx(self.canvas.bitmap, self.canvas.figure.dpi)
        self.canvas.figure.draw(renderer)
        self.canvas.bitmap.SetWidth(
            int(self.canvas.bitmap.GetWidth() * vscale))
        self.canvas.bitmap.SetHeight(
            int(self.canvas.bitmap.GetHeight() * vscale))
        self.canvas.draw()

        # page may need additional scaling on preview
        page_scale = 1.0
        if self.IsPreview():
            page_scale = float(dcw) / pgw

        # get margin in pixels = (margin in in) * (pixels/in)
        top_margin = int(self.margin * pph * page_scale)
        left_margin = int(self.margin * ppw * page_scale)

        # set scale so that width of output is self.width inches
        # (assuming grw is size of graph in inches....)
        user_scale = (self.width * fig_dpi * page_scale) / float(grw)

        dc.SetDeviceOrigin(left_margin, top_margin)
        dc.SetUserScale(user_scale, user_scale)

        # this cute little number avoid API inconsistencies in wx
        try:
            dc.DrawBitmap(self.canvas.bitmap, 0, 0)
        except Exception:
            try:
                dc.DrawBitmap(self.canvas.bitmap, (0, 0))
            except Exception:
                pass

        # restore original figure  resolution
        self.canvas.figure.set_facecolor(bgcolor)
        self.canvas.figure.dpi = fig_dpi
        self.canvas.draw()
        return True
# >


@_Backend.export
class _BackendWx(_Backend):
    FigureCanvas = FigureCanvasWx
    FigureManager = FigureManagerWx
    _frame_class = FigureFrameWx

    @staticmethod
    def trigger_manager_draw(manager):
        manager.canvas.draw_idle()

    @classmethod
    def new_figure_manager(cls, num, *args, **kwargs):
        # Create a wx.App instance if it has not been created sofar.
        wxapp = wx.GetApp()
        if wxapp is None:
            wxapp = wx.App(False)
            wxapp.SetExitOnFrameDelete(True)
            # Retain a reference to the app object so that it does not get
            # garbage collected.
            _BackendWx._theWxApp = wxapp
        return super(_BackendWx, cls).new_figure_manager(num, *args, **kwargs)

    @classmethod
    def new_figure_manager_given_figure(cls, num, figure):
        frame = cls._frame_class(num, figure)
        figmgr = frame.get_figure_manager()
        if matplotlib.is_interactive():
            figmgr.frame.Show()
            figure.canvas.draw_idle()
        return figmgr

    @staticmethod
    def mainloop():
        if not wx.App.IsMainLoopRunning():
            wxapp = wx.GetApp()
            if wxapp is not None:
                wxapp.MainLoop()
