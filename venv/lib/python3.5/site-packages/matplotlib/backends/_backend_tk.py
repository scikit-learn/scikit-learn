from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
from six.moves import tkinter as Tk

import logging
import os.path
import sys

# Paint image to Tk photo blitter extension
import matplotlib.backends.tkagg as tkagg

from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.backends.windowing as windowing

import matplotlib
from matplotlib import backend_tools, cbook, rcParams
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, NavigationToolbar2,
    StatusbarBase, TimerBase, ToolContainerBase, cursors)
from matplotlib.backend_managers import ToolManager
from matplotlib._pylab_helpers import Gcf
from matplotlib.figure import Figure
from matplotlib.widgets import SubplotTool


_log = logging.getLogger(__name__)

backend_version = Tk.TkVersion

# the true dots per inch on the screen; should be display dependent
# see http://groups.google.com/groups?q=screen+dpi+x11&hl=en&lr=&ie=UTF-8&oe=UTF-8&safe=off&selm=7077.26e81ad5%40swift.cs.tcd.ie&rnum=5 for some info about screen dpi
PIXELS_PER_INCH = 75

cursord = {
    cursors.MOVE: "fleur",
    cursors.HAND: "hand2",
    cursors.POINTER: "arrow",
    cursors.SELECT_REGION: "tcross",
    cursors.WAIT: "watch",
    }


def raise_msg_to_str(msg):
    """msg is a return arg from a raise.  Join with new lines"""
    if not isinstance(msg, six.string_types):
        msg = '\n'.join(map(str, msg))
    return msg

def error_msg_tkpaint(msg, parent=None):
    from six.moves import tkinter_messagebox as tkMessageBox
    tkMessageBox.showerror("matplotlib", msg)


class TimerTk(TimerBase):
    '''
    Subclass of :class:`backend_bases.TimerBase` that uses Tk's timer events.

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
        self.parent = parent
        self._timer = None

    def _timer_start(self):
        self._timer_stop()
        self._timer = self.parent.after(self._interval, self._on_timer)

    def _timer_stop(self):
        if self._timer is not None:
            self.parent.after_cancel(self._timer)
        self._timer = None

    def _on_timer(self):
        TimerBase._on_timer(self)

        # Tk after() is only a single shot, so we need to add code here to
        # reset the timer if we're not operating in single shot mode.  However,
        # if _timer is None, this means that _timer_stop has been called; so
        # don't recreate the timer in that case.
        if not self._single and self._timer:
            self._timer = self.parent.after(self._interval, self._on_timer)
        else:
            self._timer = None


class FigureCanvasTk(FigureCanvasBase):
    keyvald = {65507 : 'control',
               65505 : 'shift',
               65513 : 'alt',
               65515 : 'super',
               65508 : 'control',
               65506 : 'shift',
               65514 : 'alt',
               65361 : 'left',
               65362 : 'up',
               65363 : 'right',
               65364 : 'down',
               65307 : 'escape',
               65470 : 'f1',
               65471 : 'f2',
               65472 : 'f3',
               65473 : 'f4',
               65474 : 'f5',
               65475 : 'f6',
               65476 : 'f7',
               65477 : 'f8',
               65478 : 'f9',
               65479 : 'f10',
               65480 : 'f11',
               65481 : 'f12',
               65300 : 'scroll_lock',
               65299 : 'break',
               65288 : 'backspace',
               65293 : 'enter',
               65379 : 'insert',
               65535 : 'delete',
               65360 : 'home',
               65367 : 'end',
               65365 : 'pageup',
               65366 : 'pagedown',
               65438 : '0',
               65436 : '1',
               65433 : '2',
               65435 : '3',
               65430 : '4',
               65437 : '5',
               65432 : '6',
               65429 : '7',
               65431 : '8',
               65434 : '9',
               65451 : '+',
               65453 : '-',
               65450 : '*',
               65455 : '/',
               65439 : 'dec',
               65421 : 'enter',
               }

    _keycode_lookup = {
                       262145: 'control',
                       524320: 'alt',
                       524352: 'alt',
                       1048584: 'super',
                       1048592: 'super',
                       131074: 'shift',
                       131076: 'shift',
                       }
    """_keycode_lookup is used for badly mapped (i.e. no event.key_sym set)
       keys on apple keyboards."""

    def __init__(self, figure, master=None, resize_callback=None):
        super(FigureCanvasTk, self).__init__(figure)
        self._idle = True
        self._idle_callback = None
        t1,t2,w,h = self.figure.bbox.bounds
        w, h = int(w), int(h)
        self._tkcanvas = Tk.Canvas(
            master=master, background="white",
            width=w, height=h, borderwidth=0, highlightthickness=0)
        self._tkphoto = Tk.PhotoImage(
            master=self._tkcanvas, width=w, height=h)
        self._tkcanvas.create_image(w//2, h//2, image=self._tkphoto)
        self._resize_callback = resize_callback
        self._tkcanvas.bind("<Configure>", self.resize)
        self._tkcanvas.bind("<Key>", self.key_press)
        self._tkcanvas.bind("<Motion>", self.motion_notify_event)
        self._tkcanvas.bind("<KeyRelease>", self.key_release)
        for name in "<Button-1>", "<Button-2>", "<Button-3>":
            self._tkcanvas.bind(name, self.button_press_event)
        for name in "<Double-Button-1>", "<Double-Button-2>", "<Double-Button-3>":
            self._tkcanvas.bind(name, self.button_dblclick_event)
        for name in "<ButtonRelease-1>", "<ButtonRelease-2>", "<ButtonRelease-3>":
            self._tkcanvas.bind(name, self.button_release_event)

        # Mouse wheel on Linux generates button 4/5 events
        for name in "<Button-4>", "<Button-5>":
            self._tkcanvas.bind(name, self.scroll_event)
        # Mouse wheel for windows goes to the window with the focus.
        # Since the canvas won't usually have the focus, bind the
        # event to the window containing the canvas instead.
        # See http://wiki.tcl.tk/3893 (mousewheel) for details
        root = self._tkcanvas.winfo_toplevel()
        root.bind("<MouseWheel>", self.scroll_event_windows, "+")

        # Can't get destroy events by binding to _tkcanvas. Therefore, bind
        # to the window and filter.
        def filter_destroy(evt):
            if evt.widget is self._tkcanvas:
                self._master.update_idletasks()
                self.close_event()
        root.bind("<Destroy>", filter_destroy, "+")

        self._master = master
        self._tkcanvas.focus_set()

    def resize(self, event):
        width, height = event.width, event.height
        if self._resize_callback is not None:
            self._resize_callback(event)

        # compute desired figure size in inches
        dpival = self.figure.dpi
        winch = width/dpival
        hinch = height/dpival
        self.figure.set_size_inches(winch, hinch, forward=False)


        self._tkcanvas.delete(self._tkphoto)
        self._tkphoto = Tk.PhotoImage(
            master=self._tkcanvas, width=int(width), height=int(height))
        self._tkcanvas.create_image(int(width/2),int(height/2),image=self._tkphoto)
        self.resize_event()
        self.draw()

        # a resizing will in general move the pointer position
        # relative to the canvas, so process it as a motion notify
        # event.  An intended side effect of this call is to allow
        # window raises (which trigger a resize) to get the cursor
        # position to the mpl event framework so key presses which are
        # over the axes will work w/o clicks or explicit motion
        self._update_pointer_position(event)

    def _update_pointer_position(self, guiEvent=None):
        """
        Figure out if we are inside the canvas or not and update the
        canvas enter/leave events
        """
        # if the pointer if over the canvas, set the lastx and lasty
        # attrs of the canvas so it can process event w/o mouse click
        # or move

        # the window's upper, left coords in screen coords
        xw = self._tkcanvas.winfo_rootx()
        yw = self._tkcanvas.winfo_rooty()
        # the pointer's location in screen coords
        xp, yp = self._tkcanvas.winfo_pointerxy()

        # not figure out the canvas coordinates of the pointer
        xc = xp - xw
        yc = yp - yw

        # flip top/bottom
        yc = self.figure.bbox.height - yc

        # JDH: this method was written originally to get the pointer
        # location to the backend lastx and lasty attrs so that events
        # like KeyEvent can be handled without mouse events.  e.g., if
        # the cursor is already above the axes, then key presses like
        # 'g' should toggle the grid.  In order for this to work in
        # backend_bases, the canvas needs to know _lastx and _lasty.
        # There are three ways to get this info the canvas:
        #
        # 1) set it explicitly
        #
        # 2) call enter/leave events explicitly.  The downside of this
        #    in the impl below is that enter could be repeatedly
        #    triggered if the mouse is over the axes and one is
        #    resizing with the keyboard.  This is not entirely bad,
        #    because the mouse position relative to the canvas is
        #    changing, but it may be surprising to get repeated entries
        #    without leaves
        #
        # 3) process it as a motion notify event.  This also has pros
        #    and cons.  The mouse is moving relative to the window, but
        #    this may surpise an event handler writer who is getting
        #   motion_notify_events even if the mouse has not moved

        # here are the three scenarios
        if 1:
            # just manually set it
            self._lastx, self._lasty = xc, yc
        elif 0:
            # alternate implementation: process it as a motion
            FigureCanvasBase.motion_notify_event(self, xc, yc, guiEvent)
        elif 0:
            # alternate implementation -- process enter/leave events
            # instead of motion/notify
            if self.figure.bbox.contains(xc, yc):
                self.enter_notify_event(guiEvent, xy=(xc,yc))
            else:
                self.leave_notify_event(guiEvent)

    show = cbook.deprecated("2.2", name="FigureCanvasTk.show",
                            alternative="FigureCanvasTk.draw")(
                                lambda self: self.draw())

    def draw_idle(self):
        'update drawing area only if idle'
        if self._idle is False:
            return

        self._idle = False

        def idle_draw(*args):
            try:
                self.draw()
            finally:
                self._idle = True

        self._idle_callback = self._tkcanvas.after_idle(idle_draw)

    def get_tk_widget(self):
        """returns the Tk widget used to implement FigureCanvasTkAgg.
        Although the initial implementation uses a Tk canvas,  this routine
        is intended to hide that fact.
        """
        return self._tkcanvas

    def motion_notify_event(self, event):
        x = event.x
        # flipy so y=0 is bottom of canvas
        y = self.figure.bbox.height - event.y
        FigureCanvasBase.motion_notify_event(self, x, y, guiEvent=event)


    def button_press_event(self, event, dblclick=False):
        x = event.x
        # flipy so y=0 is bottom of canvas
        y = self.figure.bbox.height - event.y
        num = getattr(event, 'num', None)

        if sys.platform=='darwin':
            # 2 and 3 were reversed on the OSX platform I
            # tested under tkagg
            if   num==2: num=3
            elif num==3: num=2

        FigureCanvasBase.button_press_event(self, x, y, num, dblclick=dblclick, guiEvent=event)

    def button_dblclick_event(self,event):
        self.button_press_event(event,dblclick=True)

    def button_release_event(self, event):
        x = event.x
        # flipy so y=0 is bottom of canvas
        y = self.figure.bbox.height - event.y

        num = getattr(event, 'num', None)

        if sys.platform=='darwin':
            # 2 and 3 were reversed on the OSX platform I
            # tested under tkagg
            if   num==2: num=3
            elif num==3: num=2

        FigureCanvasBase.button_release_event(self, x, y, num, guiEvent=event)

    def scroll_event(self, event):
        x = event.x
        y = self.figure.bbox.height - event.y
        num = getattr(event, 'num', None)
        if   num==4: step = +1
        elif num==5: step = -1
        else:        step =  0

        FigureCanvasBase.scroll_event(self, x, y, step, guiEvent=event)

    def scroll_event_windows(self, event):
        """MouseWheel event processor"""
        # need to find the window that contains the mouse
        w = event.widget.winfo_containing(event.x_root, event.y_root)
        if w == self._tkcanvas:
            x = event.x_root - w.winfo_rootx()
            y = event.y_root - w.winfo_rooty()
            y = self.figure.bbox.height - y
            step = event.delta/120.
            FigureCanvasBase.scroll_event(self, x, y, step, guiEvent=event)

    def _get_key(self, event):
        val = event.keysym_num
        if val in self.keyvald:
            key = self.keyvald[val]
        elif val == 0 and sys.platform == 'darwin' and \
                                        event.keycode in self._keycode_lookup:
            key = self._keycode_lookup[event.keycode]
        elif val < 256:
            key = chr(val)
        else:
            key = None

        # add modifier keys to the key string. Bit details originate from
        # http://effbot.org/tkinterbook/tkinter-events-and-bindings.htm
        # BIT_SHIFT = 0x001; BIT_CAPSLOCK = 0x002; BIT_CONTROL = 0x004;
        # BIT_LEFT_ALT = 0x008; BIT_NUMLOCK = 0x010; BIT_RIGHT_ALT = 0x080;
        # BIT_MB_1 = 0x100; BIT_MB_2 = 0x200; BIT_MB_3 = 0x400;
        # In general, the modifier key is excluded from the modifier flag,
        # however this is not the case on "darwin", so double check that
        # we aren't adding repeat modifier flags to a modifier key.
        if sys.platform == 'win32':
            modifiers = [(17, 'alt', 'alt'),
                         (2, 'ctrl', 'control'),
                         ]
        elif sys.platform == 'darwin':
            modifiers = [(3, 'super', 'super'),
                         (4, 'alt', 'alt'),
                         (2, 'ctrl', 'control'),
                         ]
        else:
            modifiers = [(6, 'super', 'super'),
                         (3, 'alt', 'alt'),
                         (2, 'ctrl', 'control'),
                         ]

        if key is not None:
            # note, shift is not added to the keys as this is already accounted for
            for bitmask, prefix, key_name in modifiers:
                if event.state & (1 << bitmask) and key_name not in key:
                    key = '{0}+{1}'.format(prefix, key)

        return key

    def key_press(self, event):
        key = self._get_key(event)
        FigureCanvasBase.key_press_event(self, key, guiEvent=event)

    def key_release(self, event):
        key = self._get_key(event)
        FigureCanvasBase.key_release_event(self, key, guiEvent=event)

    def new_timer(self, *args, **kwargs):
        """
        Creates a new backend-specific subclass of :class:`backend_bases.Timer`.
        This is useful for getting periodic events through the backend's native
        event loop. Implemented only for backends with GUIs.

        Other Parameters
        ----------------
        interval : scalar
            Timer interval in milliseconds
        callbacks : list
            Sequence of (func, args, kwargs) where ``func(*args, **kwargs)``
            will be executed by the timer every *interval*.

        """
        return TimerTk(self._tkcanvas, *args, **kwargs)

    def flush_events(self):
        self._master.update()


class FigureManagerTk(FigureManagerBase):
    """
    Attributes
    ----------
    canvas : `FigureCanvas`
        The FigureCanvas instance
    num : int or str
        The Figure number
    toolbar : tk.Toolbar
        The tk.Toolbar
    window : tk.Window
        The tk.Window

    """
    def __init__(self, canvas, num, window):
        FigureManagerBase.__init__(self, canvas, num)
        self.window = window
        self.window.withdraw()
        self.set_window_title("Figure %d" % num)
        self.canvas = canvas
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self._num = num

        self.toolmanager = self._get_toolmanager()
        self.toolbar = self._get_toolbar()
        self.statusbar = None

        if self.toolmanager:
            backend_tools.add_tools_to_manager(self.toolmanager)
            if self.toolbar:
                backend_tools.add_tools_to_container(self.toolbar)
                self.statusbar = StatusbarTk(self.window, self.toolmanager)

        self._shown = False

        def notify_axes_change(fig):
            'this will be called whenever the current axes is changed'
            if self.toolmanager is not None:
                pass
            elif self.toolbar is not None:
                self.toolbar.update()
        self.canvas.figure.add_axobserver(notify_axes_change)

    def _get_toolbar(self):
        if matplotlib.rcParams['toolbar'] == 'toolbar2':
            toolbar = NavigationToolbar2Tk(self.canvas, self.window)
        elif matplotlib.rcParams['toolbar'] == 'toolmanager':
            toolbar = ToolbarTk(self.toolmanager, self.window)
        else:
            toolbar = None
        return toolbar

    def _get_toolmanager(self):
        if rcParams['toolbar'] == 'toolmanager':
            toolmanager = ToolManager(self.canvas.figure)
        else:
            toolmanager = None
        return toolmanager

    def resize(self, width, height=None):
        # before 09-12-22, the resize method takes a single *event*
        # parameter. On the other hand, the resize method of other
        # FigureManager class takes *width* and *height* parameter,
        # which is used to change the size of the window. For the
        # Figure.set_size_inches with forward=True work with Tk
        # backend, I changed the function signature but tried to keep
        # it backward compatible. -JJL

        # when a single parameter is given, consider it as a event
        if height is None:
            cbook.warn_deprecated("2.2", "FigureManagerTkAgg.resize now takes "
                                  "width and height as separate arguments")
            width = width.width
        else:
            self.canvas._tkcanvas.master.geometry("%dx%d" % (width, height))

        if self.toolbar is not None:
            self.toolbar.configure(width=width)

    def show(self):
        """
        this function doesn't segfault but causes the
        PyEval_RestoreThread: NULL state bug on win32
        """
        _focus = windowing.FocusManager()
        if not self._shown:
            def destroy(*args):
                self.window = None
                Gcf.destroy(self._num)
            self.canvas._tkcanvas.bind("<Destroy>", destroy)
            self.window.deiconify()
        else:
            self.canvas.draw_idle()
        # Raise the new window.
        self.canvas.manager.window.attributes('-topmost', 1)
        self.canvas.manager.window.attributes('-topmost', 0)
        self._shown = True

    def destroy(self, *args):
        if self.window is not None:
            #self.toolbar.destroy()
            if self.canvas._idle_callback:
                self.canvas._tkcanvas.after_cancel(self.canvas._idle_callback)
            self.window.destroy()
        if Gcf.get_num_fig_managers()==0:
            if self.window is not None:
                self.window.quit()
        self.window = None

    def get_window_title(self):
        return self.window.wm_title()

    def set_window_title(self, title):
        self.window.wm_title(title)

    def full_screen_toggle(self):
        is_fullscreen = bool(self.window.attributes('-fullscreen'))
        self.window.attributes('-fullscreen', not is_fullscreen)


@cbook.deprecated("2.2")
class AxisMenu(object):
    def __init__(self, master, naxes):
        self._master = master
        self._naxes = naxes
        self._mbar = Tk.Frame(master=master, relief=Tk.RAISED, borderwidth=2)
        self._mbar.pack(side=Tk.LEFT)
        self._mbutton = Tk.Menubutton(
            master=self._mbar, text="Axes", underline=0)
        self._mbutton.pack(side=Tk.LEFT, padx="2m")
        self._mbutton.menu = Tk.Menu(self._mbutton)
        self._mbutton.menu.add_command(
            label="Select All", command=self.select_all)
        self._mbutton.menu.add_command(
            label="Invert All", command=self.invert_all)
        self._axis_var = []
        self._checkbutton = []
        for i in range(naxes):
            self._axis_var.append(Tk.IntVar())
            self._axis_var[i].set(1)
            self._checkbutton.append(self._mbutton.menu.add_checkbutton(
                label = "Axis %d" % (i+1),
                variable=self._axis_var[i],
                command=self.set_active))
            self._mbutton.menu.invoke(self._mbutton.menu.index("Select All"))
        self._mbutton['menu'] = self._mbutton.menu
        self._mbar.tk_menuBar(self._mbutton)
        self.set_active()

    def adjust(self, naxes):
        if self._naxes < naxes:
            for i in range(self._naxes, naxes):
                self._axis_var.append(Tk.IntVar())
                self._axis_var[i].set(1)
                self._checkbutton.append( self._mbutton.menu.add_checkbutton(
                    label = "Axis %d" % (i+1),
                    variable=self._axis_var[i],
                    command=self.set_active))
        elif self._naxes > naxes:
            for i in range(self._naxes-1, naxes-1, -1):
                del self._axis_var[i]
                self._mbutton.menu.forget(self._checkbutton[i])
                del self._checkbutton[i]
        self._naxes = naxes
        self.set_active()

    def get_indices(self):
        a = [i for i in range(len(self._axis_var)) if self._axis_var[i].get()]
        return a

    def set_active(self):
        self._master.set_active(self.get_indices())

    def invert_all(self):
        for a in self._axis_var:
            a.set(not a.get())
        self.set_active()

    def select_all(self):
        for a in self._axis_var:
            a.set(1)
        self.set_active()


class NavigationToolbar2Tk(NavigationToolbar2, Tk.Frame):
    """
    Attributes
    ----------
    canvas : `FigureCanvas`
        the figure canvas on which to operate
    win : tk.Window
        the tk.Window which owns this toolbar

    """
    def __init__(self, canvas, window):
        self.canvas = canvas
        self.window = window
        NavigationToolbar2.__init__(self, canvas)

    def destroy(self, *args):
        del self.message
        Tk.Frame.destroy(self, *args)

    def set_message(self, s):
        self.message.set(s)

    def draw_rubberband(self, event, x0, y0, x1, y1):
        height = self.canvas.figure.bbox.height
        y0 = height - y0
        y1 = height - y1
        if hasattr(self, "lastrect"):
            self.canvas._tkcanvas.delete(self.lastrect)
        self.lastrect = self.canvas._tkcanvas.create_rectangle(x0, y0, x1, y1)

        #self.canvas.draw()

    def release(self, event):
        try: self.lastrect
        except AttributeError: pass
        else:
            self.canvas._tkcanvas.delete(self.lastrect)
            del self.lastrect

    def set_cursor(self, cursor):
        self.window.configure(cursor=cursord[cursor])
        self.window.update_idletasks()

    def _Button(self, text, file, command, extension='.gif'):
        img_file = os.path.join(
            rcParams['datapath'], 'images', file + extension)
        im = Tk.PhotoImage(master=self, file=img_file)
        b = Tk.Button(
            master=self, text=text, padx=2, pady=2, image=im, command=command)
        b._ntimage = im
        b.pack(side=Tk.LEFT)
        return b

    def _Spacer(self):
        # Buttons are 30px high, so make this 26px tall with padding to center it
        s = Tk.Frame(
            master=self, height=26, relief=Tk.RIDGE, pady=2, bg="DarkGray")
        s.pack(side=Tk.LEFT, padx=5)
        return s

    def _init_toolbar(self):
        xmin, xmax = self.canvas.figure.bbox.intervalx
        height, width = 50, xmax-xmin
        Tk.Frame.__init__(self, master=self.window,
                          width=int(width), height=int(height),
                          borderwidth=2)

        self.update()  # Make axes menu

        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                # Add a spacer; return value is unused.
                self._Spacer()
            else:
                button = self._Button(text=text, file=image_file,
                                      command=getattr(self, callback))
                if tooltip_text is not None:
                    ToolTip.createToolTip(button, tooltip_text)

        self.message = Tk.StringVar(master=self)
        self._message_label = Tk.Label(master=self, textvariable=self.message)
        self._message_label.pack(side=Tk.RIGHT)
        self.pack(side=Tk.BOTTOM, fill=Tk.X)

    def configure_subplots(self):
        toolfig = Figure(figsize=(6,3))
        window = Tk.Toplevel()
        canvas = type(self.canvas)(toolfig, master=window)
        toolfig.subplots_adjust(top=0.9)
        canvas.tool = SubplotTool(self.canvas.figure, toolfig)
        canvas.draw()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        window.grab_set()

    def save_figure(self, *args):
        from six.moves import tkinter_tkfiledialog, tkinter_messagebox
        filetypes = self.canvas.get_supported_filetypes().copy()
        default_filetype = self.canvas.get_default_filetype()

        # Tk doesn't provide a way to choose a default filetype,
        # so we just have to put it first
        default_filetype_name = filetypes.pop(default_filetype)
        sorted_filetypes = ([(default_filetype, default_filetype_name)]
                            + sorted(six.iteritems(filetypes)))
        tk_filetypes = [(name, '*.%s' % ext) for ext, name in sorted_filetypes]

        # adding a default extension seems to break the
        # asksaveasfilename dialog when you choose various save types
        # from the dropdown.  Passing in the empty string seems to
        # work - JDH!
        #defaultextension = self.canvas.get_default_filetype()
        defaultextension = ''
        initialdir = os.path.expanduser(rcParams['savefig.directory'])
        initialfile = self.canvas.get_default_filename()
        fname = tkinter_tkfiledialog.asksaveasfilename(
            master=self.window,
            title='Save the figure',
            filetypes=tk_filetypes,
            defaultextension=defaultextension,
            initialdir=initialdir,
            initialfile=initialfile,
            )

        if fname in ["", ()]:
            return
        # Save dir for next time, unless empty str (i.e., use cwd).
        if initialdir != "":
            rcParams['savefig.directory'] = (
                os.path.dirname(six.text_type(fname)))
        try:
            # This method will handle the delegation to the correct type
            self.canvas.figure.savefig(fname)
        except Exception as e:
            tkinter_messagebox.showerror("Error saving file", str(e))

    def set_active(self, ind):
        self._ind = ind
        self._active = [self._axes[i] for i in self._ind]

    def update(self):
        _focus = windowing.FocusManager()
        self._axes = self.canvas.figure.axes
        NavigationToolbar2.update(self)


class ToolTip(object):
    """
    Tooltip recipe from
    http://www.voidspace.org.uk/python/weblog/arch_d7_2006_07_01.shtml#e387
    """
    @staticmethod
    def createToolTip(widget, text):
        toolTip = ToolTip(widget)
        def enter(event):
            toolTip.showtip(text)
        def leave(event):
            toolTip.hidetip()
        widget.bind('<Enter>', enter)
        widget.bind('<Leave>', leave)

    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, _, _ = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 27
        y = y + self.widget.winfo_rooty()
        self.tipwindow = tw = Tk.Toplevel(self.widget)
        tw.wm_overrideredirect(1)
        tw.wm_geometry("+%d+%d" % (x, y))
        try:
            # For Mac OS
            tw.tk.call("::tk::unsupported::MacWindowStyle",
                       "style", tw._w,
                       "help", "noActivates")
        except Tk.TclError:
            pass
        label = Tk.Label(tw, text=self.text, justify=Tk.LEFT,
                         background="#ffffe0", relief=Tk.SOLID, borderwidth=1)
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()


class RubberbandTk(backend_tools.RubberbandBase):
    def __init__(self, *args, **kwargs):
        backend_tools.RubberbandBase.__init__(self, *args, **kwargs)

    def draw_rubberband(self, x0, y0, x1, y1):
        height = self.figure.canvas.figure.bbox.height
        y0 = height - y0
        y1 = height - y1
        if hasattr(self, "lastrect"):
            self.figure.canvas._tkcanvas.delete(self.lastrect)
        self.lastrect = self.figure.canvas._tkcanvas.create_rectangle(
            x0, y0, x1, y1)

    def remove_rubberband(self):
        if hasattr(self, "lastrect"):
            self.figure.canvas._tkcanvas.delete(self.lastrect)
            del self.lastrect


class SetCursorTk(backend_tools.SetCursorBase):
    def set_cursor(self, cursor):
        self.figure.canvas.manager.window.configure(cursor=cursord[cursor])


class ToolbarTk(ToolContainerBase, Tk.Frame):
    _icon_extension = '.gif'
    def __init__(self, toolmanager, window):
        ToolContainerBase.__init__(self, toolmanager)
        xmin, xmax = self.toolmanager.canvas.figure.bbox.intervalx
        height, width = 50, xmax - xmin
        Tk.Frame.__init__(self, master=window,
                          width=int(width), height=int(height),
                          borderwidth=2)
        self._toolitems = {}
        self.pack(side=Tk.TOP, fill=Tk.X)
        self._groups = {}

    def add_toolitem(
            self, name, group, position, image_file, description, toggle):
        frame = self._get_groupframe(group)
        button = self._Button(name, image_file, toggle, frame)
        if description is not None:
            ToolTip.createToolTip(button, description)
        self._toolitems.setdefault(name, [])
        self._toolitems[name].append(button)

    def _get_groupframe(self, group):
        if group not in self._groups:
            if self._groups:
                self._add_separator()
            frame = Tk.Frame(master=self, borderwidth=0)
            frame.pack(side=Tk.LEFT, fill=Tk.Y)
            self._groups[group] = frame
        return self._groups[group]

    def _add_separator(self):
        separator = Tk.Frame(master=self, bd=5, width=1, bg='black')
        separator.pack(side=Tk.LEFT, fill=Tk.Y, padx=2)

    def _Button(self, text, image_file, toggle, frame):
        if image_file is not None:
            im = Tk.PhotoImage(master=self, file=image_file)
        else:
            im = None

        if not toggle:
            b = Tk.Button(master=frame, text=text, padx=2, pady=2, image=im,
                          command=lambda: self._button_click(text))
        else:
            # There is a bug in tkinter included in some python 3.6 versions
            # that without this variable, produces a "visual" toggling of
            # other near checkbuttons
            # https://bugs.python.org/issue29402
            # https://bugs.python.org/issue25684
            var = Tk.IntVar()
            b = Tk.Checkbutton(master=frame, text=text, padx=2, pady=2,
                               image=im, indicatoron=False,
                               command=lambda: self._button_click(text),
                               variable=var)
        b._ntimage = im
        b.pack(side=Tk.LEFT)
        return b

    def _button_click(self, name):
        self.trigger_tool(name)

    def toggle_toolitem(self, name, toggled):
        if name not in self._toolitems:
            return
        for toolitem in self._toolitems[name]:
            if toggled:
                toolitem.select()
            else:
                toolitem.deselect()

    def remove_toolitem(self, name):
        for toolitem in self._toolitems[name]:
            toolitem.pack_forget()
        del self._toolitems[name]


class StatusbarTk(StatusbarBase, Tk.Frame):
    def __init__(self, window, *args, **kwargs):
        StatusbarBase.__init__(self, *args, **kwargs)
        xmin, xmax = self.toolmanager.canvas.figure.bbox.intervalx
        height, width = 50, xmax - xmin
        Tk.Frame.__init__(self, master=window,
                          width=int(width), height=int(height),
                          borderwidth=2)
        self._message = Tk.StringVar(master=self)
        self._message_label = Tk.Label(master=self, textvariable=self._message)
        self._message_label.pack(side=Tk.RIGHT)
        self.pack(side=Tk.TOP, fill=Tk.X)

    def set_message(self, s):
        self._message.set(s)


class SaveFigureTk(backend_tools.SaveFigureBase):
    def trigger(self, *args):
        from six.moves import tkinter_tkfiledialog, tkinter_messagebox
        filetypes = self.figure.canvas.get_supported_filetypes().copy()
        default_filetype = self.figure.canvas.get_default_filetype()

        # Tk doesn't provide a way to choose a default filetype,
        # so we just have to put it first
        default_filetype_name = filetypes.pop(default_filetype)
        sorted_filetypes = ([(default_filetype, default_filetype_name)]
                            + sorted(six.iteritems(filetypes)))
        tk_filetypes = [(name, '*.%s' % ext) for ext, name in sorted_filetypes]

        # adding a default extension seems to break the
        # asksaveasfilename dialog when you choose various save types
        # from the dropdown.  Passing in the empty string seems to
        # work - JDH!
        # defaultextension = self.figure.canvas.get_default_filetype()
        defaultextension = ''
        initialdir = os.path.expanduser(rcParams['savefig.directory'])
        initialfile = self.figure.canvas.get_default_filename()
        fname = tkinter_tkfiledialog.asksaveasfilename(
            master=self.figure.canvas.manager.window,
            title='Save the figure',
            filetypes=tk_filetypes,
            defaultextension=defaultextension,
            initialdir=initialdir,
            initialfile=initialfile,
            )

        if fname == "" or fname == ():
            return
        else:
            if initialdir == '':
                # explicitly missing key or empty str signals to use cwd
                rcParams['savefig.directory'] = initialdir
            else:
                # save dir for next time
                rcParams['savefig.directory'] = os.path.dirname(
                    six.text_type(fname))
            try:
                # This method will handle the delegation to the correct type
                self.figure.savefig(fname)
            except Exception as e:
                tkinter_messagebox.showerror("Error saving file", str(e))


class ConfigureSubplotsTk(backend_tools.ConfigureSubplotsBase):
    def __init__(self, *args, **kwargs):
        backend_tools.ConfigureSubplotsBase.__init__(self, *args, **kwargs)
        self.window = None

    def trigger(self, *args):
        self.init_window()
        self.window.lift()

    def init_window(self):
        if self.window:
            return

        toolfig = Figure(figsize=(6, 3))
        self.window = Tk.Tk()

        canvas = type(self.canvas)(toolfig, master=self.window)
        toolfig.subplots_adjust(top=0.9)
        _tool = SubplotTool(self.figure, toolfig)
        canvas.draw()
        canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
        self.window.protocol("WM_DELETE_WINDOW", self.destroy)

    def destroy(self, *args, **kwargs):
        self.window.destroy()
        self.window = None


backend_tools.ToolSaveFigure = SaveFigureTk
backend_tools.ToolConfigureSubplots = ConfigureSubplotsTk
backend_tools.ToolSetCursor = SetCursorTk
backend_tools.ToolRubberband = RubberbandTk
Toolbar = ToolbarTk


@_Backend.export
class _BackendTk(_Backend):
    FigureManager = FigureManagerTk

    @classmethod
    def new_figure_manager_given_figure(cls, num, figure):
        """
        Create a new figure manager instance for the given figure.
        """
        _focus = windowing.FocusManager()
        window = Tk.Tk(className="matplotlib")
        window.withdraw()

        # Put a mpl icon on the window rather than the default tk icon.
        # Tkinter doesn't allow colour icons on linux systems, but tk>=8.5 has
        # a iconphoto command which we call directly. Source:
        # http://mail.python.org/pipermail/tkinter-discuss/2006-November/000954.html
        icon_fname = os.path.join(
            rcParams['datapath'], 'images', 'matplotlib.ppm')
        icon_img = Tk.PhotoImage(file=icon_fname)
        try:
            window.tk.call('wm', 'iconphoto', window._w, icon_img)
        except Exception as exc:
            # log the failure (due e.g. to Tk version), but carry on
            _log.info('Could not load matplotlib icon: %s', exc)

        canvas = cls.FigureCanvas(figure, master=window)
        manager = cls.FigureManager(canvas, num, window)
        if matplotlib.is_interactive():
            manager.show()
            canvas.draw_idle()
        return manager

    @staticmethod
    def trigger_manager_draw(manager):
        manager.show()

    @staticmethod
    def mainloop():
        Tk.mainloop()
