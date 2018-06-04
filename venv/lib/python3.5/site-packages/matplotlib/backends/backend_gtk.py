from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six

import logging
import os
import sys
import warnings

if six.PY3:
    warnings.warn(
        "The gtk* backends have not been tested with Python 3.x",
        ImportWarning)

try:
    import gobject
    import gtk; gdk = gtk.gdk
    import pango
except ImportError:
    raise ImportError("Gtk* backend requires pygtk to be installed.")

pygtk_version_required = (2,4,0)
if gtk.pygtk_version < pygtk_version_required:
    raise ImportError ("PyGTK %d.%d.%d is installed\n"
                      "PyGTK %d.%d.%d or later is required"
                      % (gtk.pygtk_version + pygtk_version_required))
del pygtk_version_required

_new_tooltip_api =  (gtk.pygtk_version[1] >= 12)

import matplotlib
from matplotlib._pylab_helpers import Gcf
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, NavigationToolbar2,
    TimerBase, cursors)

from matplotlib.backends.backend_gdk import RendererGDK, FigureCanvasGDK
from matplotlib.cbook import is_writable_file_like, warn_deprecated
from matplotlib.figure import Figure
from matplotlib.widgets import SubplotTool

from matplotlib import (
    cbook, colors as mcolors, lines, markers, rcParams)

_log = logging.getLogger(__name__)

backend_version = "%d.%d.%d" % gtk.pygtk_version

# the true dots per inch on the screen; should be display dependent
# see http://groups.google.com/groups?q=screen+dpi+x11&hl=en&lr=&ie=UTF-8&oe=UTF-8&safe=off&selm=7077.26e81ad5%40swift.cs.tcd.ie&rnum=5 for some info about screen dpi
PIXELS_PER_INCH = 96

# Hide the benign warning that it can't stat a file that doesn't
warnings.filterwarnings('ignore', '.*Unable to retrieve the file info for.*', gtk.Warning)

cursord = {
    cursors.MOVE          : gdk.Cursor(gdk.FLEUR),
    cursors.HAND          : gdk.Cursor(gdk.HAND2),
    cursors.POINTER       : gdk.Cursor(gdk.LEFT_PTR),
    cursors.SELECT_REGION : gdk.Cursor(gdk.TCROSS),
    cursors.WAIT          : gdk.Cursor(gdk.WATCH),
    }

# ref gtk+/gtk/gtkwidget.h
def GTK_WIDGET_DRAWABLE(w):
    flags = w.flags();
    return flags & gtk.VISIBLE != 0 and flags & gtk.MAPPED != 0


class TimerGTK(TimerBase):
    '''
    Subclass of :class:`backend_bases.TimerBase` using GTK for timer events.

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
    def _timer_start(self):
        # Need to stop it, otherwise we potentially leak a timer id that will
        # never be stopped.
        self._timer_stop()
        self._timer = gobject.timeout_add(self._interval, self._on_timer)

    def _timer_stop(self):
        if self._timer is not None:
            gobject.source_remove(self._timer)
            self._timer = None

    def _timer_set_interval(self):
        # Only stop and restart it if the timer has already been started
        if self._timer is not None:
            self._timer_stop()
            self._timer_start()

    def _on_timer(self):
        TimerBase._on_timer(self)

        # Gtk timeout_add() requires that the callback returns True if it
        # is to be called again.
        if len(self.callbacks) > 0 and not self._single:
            return True
        else:
            self._timer = None
            return False


class FigureCanvasGTK (gtk.DrawingArea, FigureCanvasBase):
    keyvald = {65507 : 'control',
               65505 : 'shift',
               65513 : 'alt',
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
               65511 : 'super',
               65512 : 'super',
               65406 : 'alt',
               65289 : 'tab',
               }

    # Setting this as a static constant prevents
    # this resulting expression from leaking
    event_mask = (gdk.BUTTON_PRESS_MASK   |
                  gdk.BUTTON_RELEASE_MASK |
                  gdk.EXPOSURE_MASK       |
                  gdk.KEY_PRESS_MASK      |
                  gdk.KEY_RELEASE_MASK    |
                  gdk.ENTER_NOTIFY_MASK   |
                  gdk.LEAVE_NOTIFY_MASK   |
                  gdk.POINTER_MOTION_MASK |
                  gdk.POINTER_MOTION_HINT_MASK)

    def __init__(self, figure):
        if self.__class__ == matplotlib.backends.backend_gtk.FigureCanvasGTK:
            warn_deprecated('2.0', message="The GTK backend is "
                            "deprecated. It is untested, known to be "
                            "broken and will be removed in Matplotlib 3.0. "
                            "Use the GTKAgg backend instead. "
                            "See Matplotlib usage FAQ for"
                            " more info on backends.",
                            alternative="GTKAgg")
        FigureCanvasBase.__init__(self, figure)
        gtk.DrawingArea.__init__(self)

        self._idle_draw_id  = 0
        self._need_redraw   = True
        self._pixmap_width  = -1
        self._pixmap_height = -1
        self._lastCursor    = None

        self.connect('scroll_event',         self.scroll_event)
        self.connect('button_press_event',   self.button_press_event)
        self.connect('button_release_event', self.button_release_event)
        self.connect('configure_event',      self.configure_event)
        self.connect('expose_event',         self.expose_event)
        self.connect('key_press_event',      self.key_press_event)
        self.connect('key_release_event',    self.key_release_event)
        self.connect('motion_notify_event',  self.motion_notify_event)
        self.connect('leave_notify_event',   self.leave_notify_event)
        self.connect('enter_notify_event',   self.enter_notify_event)

        self.set_events(self.__class__.event_mask)

        self.set_double_buffered(False)
        self.set_flags(gtk.CAN_FOCUS)
        self._renderer_init()

        self.last_downclick = {}

    def destroy(self):
        #gtk.DrawingArea.destroy(self)
        self.close_event()
        if self._idle_draw_id != 0:
            gobject.source_remove(self._idle_draw_id)

    def scroll_event(self, widget, event):
        x = event.x
        # flipy so y=0 is bottom of canvas
        y = self.allocation.height - event.y
        if event.direction==gdk.SCROLL_UP:
            step = 1
        else:
            step = -1
        FigureCanvasBase.scroll_event(self, x, y, step, guiEvent=event)
        return False  # finish event propagation?

    def button_press_event(self, widget, event):
        x = event.x
        # flipy so y=0 is bottom of canvas
        y = self.allocation.height - event.y
        dblclick = (event.type == gdk._2BUTTON_PRESS)
        if not dblclick:
            # GTK is the only backend that generates a DOWN-UP-DOWN-DBLCLICK-UP  event
            # sequence for a double click.  All other backends have a DOWN-UP-DBLCLICK-UP
            # sequence.  In order to provide consistency to matplotlib users, we will
            # eat the extra DOWN event in the case that we detect it is part of a double
            # click.
            # first, get the double click time in milliseconds.
            current_time  = event.get_time()
            last_time     = self.last_downclick.get(event.button,0)
            dblclick_time = gtk.settings_get_for_screen(gdk.screen_get_default()).get_property('gtk-double-click-time')
            delta_time    = current_time-last_time
            if delta_time < dblclick_time:
                del self.last_downclick[event.button] # we do not want to eat more than one event.
                return False                          # eat.
            self.last_downclick[event.button] = current_time
        FigureCanvasBase.button_press_event(self, x, y, event.button, dblclick=dblclick, guiEvent=event)
        return False  # finish event propagation?

    def button_release_event(self, widget, event):
        x = event.x
        # flipy so y=0 is bottom of canvas
        y = self.allocation.height - event.y
        FigureCanvasBase.button_release_event(self, x, y, event.button, guiEvent=event)
        return False  # finish event propagation?

    def key_press_event(self, widget, event):
        key = self._get_key(event)
        FigureCanvasBase.key_press_event(self, key, guiEvent=event)
        return True  # stop event propagation

    def key_release_event(self, widget, event):
        key = self._get_key(event)
        FigureCanvasBase.key_release_event(self, key, guiEvent=event)
        return True  # stop event propagation

    def motion_notify_event(self, widget, event):
        if event.is_hint:
            x, y, state = event.window.get_pointer()
        else:
            x, y, state = event.x, event.y, event.state

        # flipy so y=0 is bottom of canvas
        y = self.allocation.height - y
        FigureCanvasBase.motion_notify_event(self, x, y, guiEvent=event)
        return False  # finish event propagation?

    def leave_notify_event(self, widget, event):
        FigureCanvasBase.leave_notify_event(self, event)

    def enter_notify_event(self, widget, event):
        x, y, state = event.window.get_pointer()
        FigureCanvasBase.enter_notify_event(self, event, xy=(x, y))

    def _get_key(self, event):
        if event.keyval in self.keyvald:
            key = self.keyvald[event.keyval]
        elif event.keyval < 256:
            key = chr(event.keyval)
        else:
            key = None

        for key_mask, prefix in (
                                 [gdk.MOD4_MASK, 'super'],
                                 [gdk.MOD1_MASK, 'alt'],
                                 [gdk.CONTROL_MASK, 'ctrl'], ):
            if event.state & key_mask:
                key = '{0}+{1}'.format(prefix, key)

        return key

    def configure_event(self, widget, event):
        if widget.window is None:
            return
        w, h = event.width, event.height
        if w < 3 or h < 3:
            return # empty fig

        # resize the figure (in inches)
        dpi = self.figure.dpi
        self.figure.set_size_inches(w/dpi, h/dpi, forward=False)
        self._need_redraw = True

        return False  # finish event propagation?

    def draw(self):
        # Note: FigureCanvasBase.draw() is inconveniently named as it clashes
        # with the deprecated gtk.Widget.draw()

        self._need_redraw = True
        if GTK_WIDGET_DRAWABLE(self):
            self.queue_draw()
            # do a synchronous draw (its less efficient than an async draw,
            # but is required if/when animation is used)
            self.window.process_updates (False)

    def draw_idle(self):
        if self._idle_draw_id != 0:
            return
        def idle_draw(*args):
            try:
                self.draw()
            finally:
                self._idle_draw_id = 0
            return False
        self._idle_draw_id = gobject.idle_add(idle_draw)


    def _renderer_init(self):
        """Override by GTK backends to select a different renderer
        Renderer should provide the methods:
            set_pixmap ()
            set_width_height ()
        that are used by
            _render_figure() / _pixmap_prepare()
        """
        self._renderer = RendererGDK (self, self.figure.dpi)


    def _pixmap_prepare(self, width, height):
        """
        Make sure _._pixmap is at least width, height,
        create new pixmap if necessary
        """
        create_pixmap = False
        if width > self._pixmap_width:
            # increase the pixmap in 10%+ (rather than 1 pixel) steps
            self._pixmap_width  = max (int (self._pixmap_width  * 1.1),
                                       width)
            create_pixmap = True

        if height > self._pixmap_height:
            self._pixmap_height = max (int (self._pixmap_height * 1.1),
                                           height)
            create_pixmap = True

        if create_pixmap:
            self._pixmap = gdk.Pixmap (self.window, self._pixmap_width,
                                       self._pixmap_height)
            self._renderer.set_pixmap (self._pixmap)


    def _render_figure(self, pixmap, width, height):
        """used by GTK and GTKcairo. GTKAgg overrides
        """
        self._renderer.set_width_height (width, height)
        self.figure.draw (self._renderer)


    def expose_event(self, widget, event):
        """Expose_event for all GTK backends. Should not be overridden.
        """
        toolbar = self.toolbar
        # if toolbar:
        #     toolbar.set_cursor(cursors.WAIT)
        if GTK_WIDGET_DRAWABLE(self):
            if self._need_redraw:
                x, y, w, h = self.allocation
                self._pixmap_prepare (w, h)
                self._render_figure(self._pixmap, w, h)
                self._need_redraw = False
            x, y, w, h = event.area
            self.window.draw_drawable (self.style.fg_gc[self.state],
                                       self._pixmap, x, y, x, y, w, h)
        # if toolbar:
        #     toolbar.set_cursor(toolbar._lastCursor)
        return False  # finish event propagation?

    filetypes = FigureCanvasBase.filetypes.copy()
    filetypes['jpg'] = 'JPEG'
    filetypes['jpeg'] = 'JPEG'
    filetypes['png'] = 'Portable Network Graphics'

    def print_jpeg(self, filename, *args, **kwargs):
        return self._print_image(filename, 'jpeg')
    print_jpg = print_jpeg

    def print_png(self, filename, *args, **kwargs):
        return self._print_image(filename, 'png')

    def _print_image(self, filename, format, *args, **kwargs):
        if self.flags() & gtk.REALIZED == 0:
            # for self.window(for pixmap) and has a side effect of altering
            # figure width,height (via configure-event?)
            gtk.DrawingArea.realize(self)

        width, height = self.get_width_height()
        pixmap = gdk.Pixmap (self.window, width, height)
        self._renderer.set_pixmap (pixmap)
        self._render_figure(pixmap, width, height)

        # jpg colors don't match the display very well, png colors match
        # better
        pixbuf = gdk.Pixbuf(gdk.COLORSPACE_RGB, 0, 8, width, height)
        pixbuf.get_from_drawable(pixmap, pixmap.get_colormap(),
                                     0, 0, 0, 0, width, height)

        # set the default quality, if we are writing a JPEG.
        # http://www.pygtk.org/docs/pygtk/class-gdkpixbuf.html#method-gdkpixbuf--save
        options = {k: kwargs[k] for k in ['quality'] if k in kwargs}
        if format in ['jpg', 'jpeg']:
            options.setdefault('quality', rcParams['savefig.jpeg_quality'])
            options['quality'] = str(options['quality'])

        if isinstance(filename, six.string_types):
            try:
                pixbuf.save(filename, format, options=options)
            except gobject.GError as exc:
                error_msg_gtk('Save figure failure:\n%s' % (exc,), parent=self)
        elif is_writable_file_like(filename):
            if hasattr(pixbuf, 'save_to_callback'):
                def save_callback(buf, data=None):
                    data.write(buf)
                try:
                    pixbuf.save_to_callback(save_callback, format, user_data=filename, options=options)
                except gobject.GError as exc:
                    error_msg_gtk('Save figure failure:\n%s' % (exc,), parent=self)
            else:
                raise ValueError("Saving to a Python file-like object is only supported by PyGTK >= 2.8")
        else:
            raise ValueError("filename must be a path or a file-like object")

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
        return TimerGTK(*args, **kwargs)

    def flush_events(self):
        gtk.gdk.threads_enter()
        while gtk.events_pending():
            gtk.main_iteration(True)
        gtk.gdk.flush()
        gtk.gdk.threads_leave()


class FigureManagerGTK(FigureManagerBase):
    """
    Attributes
    ----------
    canvas : `FigureCanvas`
        The FigureCanvas instance
    num : int or str
        The Figure number
    toolbar : gtk.Toolbar
        The gtk.Toolbar  (gtk only)
    vbox : gtk.VBox
        The gtk.VBox containing the canvas and toolbar (gtk only)
    window : gtk.Window
        The gtk.Window   (gtk only)

    """
    def __init__(self, canvas, num):
        FigureManagerBase.__init__(self, canvas, num)

        self.window = gtk.Window()
        self.window.set_wmclass("matplotlib", "Matplotlib")
        self.set_window_title("Figure %d" % num)
        if window_icon:
            try:
                self.window.set_icon_from_file(window_icon)
            except:
                # some versions of gtk throw a glib.GError but not
                # all, so I am not sure how to catch it.  I am unhappy
                # diong a blanket catch here, but an not sure what a
                # better way is - JDH
                _log.info('Could not load matplotlib '
                         'icon: %s', sys.exc_info()[1])

        self.vbox = gtk.VBox()
        self.window.add(self.vbox)
        self.vbox.show()

        self.canvas.show()

        self.vbox.pack_start(self.canvas, True, True)

        self.toolbar = self._get_toolbar(canvas)

        # calculate size for window
        w = int (self.canvas.figure.bbox.width)
        h = int (self.canvas.figure.bbox.height)

        if self.toolbar is not None:
            self.toolbar.show()
            self.vbox.pack_end(self.toolbar, False, False)

            tb_w, tb_h = self.toolbar.size_request()
            h += tb_h
        self.window.set_default_size (w, h)

        def destroy(*args):
            Gcf.destroy(num)
        self.window.connect("destroy", destroy)
        self.window.connect("delete_event", destroy)
        if matplotlib.is_interactive():
            self.window.show()
            self.canvas.draw_idle()

        def notify_axes_change(fig):
            'this will be called whenever the current axes is changed'
            if self.toolbar is not None: self.toolbar.update()
        self.canvas.figure.add_axobserver(notify_axes_change)

        self.canvas.grab_focus()

    def destroy(self, *args):
        if hasattr(self, 'toolbar') and self.toolbar is not None:
            self.toolbar.destroy()
        if hasattr(self, 'vbox'):
            self.vbox.destroy()
        if hasattr(self, 'window'):
            self.window.destroy()
        if hasattr(self, 'canvas'):
            self.canvas.destroy()
        self.__dict__.clear()   #Is this needed? Other backends don't have it.

        if Gcf.get_num_fig_managers()==0 and \
               not matplotlib.is_interactive() and \
               gtk.main_level() >= 1:
            gtk.main_quit()

    def show(self):
        # show the figure window
        self.window.show()
        # raise the window above others and release the "above lock"
        self.window.set_keep_above(True)
        self.window.set_keep_above(False)

    def full_screen_toggle(self):
        self._full_screen_flag = not self._full_screen_flag
        if self._full_screen_flag:
            self.window.fullscreen()
        else:
            self.window.unfullscreen()
    _full_screen_flag = False


    def _get_toolbar(self, canvas):
        # must be inited after the window, drawingArea and figure
        # attrs are set
        if rcParams['toolbar'] == 'toolbar2':
            toolbar = NavigationToolbar2GTK (canvas, self.window)
        else:
            toolbar = None
        return toolbar

    def get_window_title(self):
        return self.window.get_title()

    def set_window_title(self, title):
        self.window.set_title(title)

    def resize(self, width, height):
        'set the canvas size in pixels'
        #_, _, cw, ch = self.canvas.allocation
        #_, _, ww, wh = self.window.allocation
        #self.window.resize (width-cw+ww, height-ch+wh)
        self.window.resize(width, height)


class NavigationToolbar2GTK(NavigationToolbar2, gtk.Toolbar):
    def __init__(self, canvas, window):
        self.win = window
        gtk.Toolbar.__init__(self)
        NavigationToolbar2.__init__(self, canvas)

    def set_message(self, s):
        self.message.set_label(s)

    def set_cursor(self, cursor):
        self.canvas.window.set_cursor(cursord[cursor])
        gtk.main_iteration()

    def release(self, event):
        try: del self._pixmapBack
        except AttributeError: pass

    def draw_rubberband(self, event, x0, y0, x1, y1):
        'adapted from http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/189744'
        drawable = self.canvas.window
        if drawable is None:
            return

        gc = drawable.new_gc()

        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0

        w = abs(x1 - x0)
        h = abs(y1 - y0)

        rect = [int(val)for val in (min(x0,x1), min(y0, y1), w, h)]
        try:
            lastrect, pixmapBack = self._pixmapBack
        except AttributeError:
            #snap image back
            if event.inaxes is None:
                return

            ax = event.inaxes
            l,b,w,h = [int(val) for val in ax.bbox.bounds]
            b = int(height)-(b+h)
            axrect = l,b,w,h
            self._pixmapBack = axrect, gtk.gdk.Pixmap(drawable, w, h)
            self._pixmapBack[1].draw_drawable(gc, drawable, l, b, 0, 0, w, h)
        else:
            drawable.draw_drawable(gc, pixmapBack, 0, 0, *lastrect)
        drawable.draw_rectangle(gc, False, *rect)


    def _init_toolbar(self):
        self.set_style(gtk.TOOLBAR_ICONS)
        self._init_toolbar2_4()


    def _init_toolbar2_4(self):
        basedir = os.path.join(rcParams['datapath'],'images')
        if not _new_tooltip_api:
            self.tooltips = gtk.Tooltips()

        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.insert( gtk.SeparatorToolItem(), -1 )
                continue
            fname = os.path.join(basedir, image_file + '.png')
            image = gtk.Image()
            image.set_from_file(fname)
            tbutton = gtk.ToolButton(image, text)
            self.insert(tbutton, -1)
            tbutton.connect('clicked', getattr(self, callback))
            if _new_tooltip_api:
                tbutton.set_tooltip_text(tooltip_text)
            else:
                tbutton.set_tooltip(self.tooltips, tooltip_text, 'Private')

        toolitem = gtk.SeparatorToolItem()
        self.insert(toolitem, -1)
        # set_draw() not making separator invisible,
        # bug #143692 fixed Jun 06 2004, will be in GTK+ 2.6
        toolitem.set_draw(False)
        toolitem.set_expand(True)

        toolitem = gtk.ToolItem()
        self.insert(toolitem, -1)
        self.message = gtk.Label()
        toolitem.add(self.message)

        self.show_all()

    def get_filechooser(self):
        fc = FileChooserDialog(
            title='Save the figure',
            parent=self.win,
            path=os.path.expanduser(rcParams['savefig.directory']),
            filetypes=self.canvas.get_supported_filetypes(),
            default_filetype=self.canvas.get_default_filetype())
        fc.set_current_name(self.canvas.get_default_filename())
        return fc

    def save_figure(self, *args):
        chooser = self.get_filechooser()
        fname, format = chooser.get_filename_from_user()
        chooser.destroy()
        if fname:
            startpath = os.path.expanduser(rcParams['savefig.directory'])
            # Save dir for next time, unless empty str (i.e., use cwd).
            if startpath != "":
                rcParams['savefig.directory'] = (
                    os.path.dirname(six.text_type(fname)))
            try:
                self.canvas.figure.savefig(fname, format=format)
            except Exception as e:
                error_msg_gtk(str(e), parent=self)

    def configure_subplots(self, button):
        toolfig = Figure(figsize=(6,3))
        canvas = self._get_canvas(toolfig)
        toolfig.subplots_adjust(top=0.9)
        tool =  SubplotTool(self.canvas.figure, toolfig)

        w = int(toolfig.bbox.width)
        h = int(toolfig.bbox.height)

        window = gtk.Window()
        if window_icon:
            try:
                window.set_icon_from_file(window_icon)
            except:
                # we presumably already logged a message on the
                # failure of the main plot, don't keep reporting
                pass
        window.set_title("Subplot Configuration Tool")
        window.set_default_size(w, h)
        vbox = gtk.VBox()
        window.add(vbox)
        vbox.show()

        canvas.show()
        vbox.pack_start(canvas, True, True)
        window.show()

    def _get_canvas(self, fig):
        return FigureCanvasGTK(fig)


class FileChooserDialog(gtk.FileChooserDialog):
    """GTK+ 2.4 file selector which presents the user with a menu
    of supported image formats
    """
    def __init__ (self,
                  title   = 'Save file',
                  parent  = None,
                  action  = gtk.FILE_CHOOSER_ACTION_SAVE,
                  buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL,
                             gtk.STOCK_SAVE,   gtk.RESPONSE_OK),
                  path    = None,
                  filetypes = [],
                  default_filetype = None
                  ):
        super(FileChooserDialog, self).__init__(title, parent, action, buttons)
        super(FileChooserDialog, self).set_do_overwrite_confirmation(True)
        self.set_default_response(gtk.RESPONSE_OK)

        if not path:
            path = os.getcwd() + os.sep

        # create an extra widget to list supported image formats
        self.set_current_folder (path)
        self.set_current_name ('image.' + default_filetype)

        hbox = gtk.HBox(spacing=10)
        hbox.pack_start(gtk.Label ("File Format:"), expand=False)

        liststore = gtk.ListStore(gobject.TYPE_STRING)
        cbox = gtk.ComboBox(liststore)
        cell = gtk.CellRendererText()
        cbox.pack_start(cell, True)
        cbox.add_attribute(cell, 'text', 0)
        hbox.pack_start(cbox)

        self.filetypes = filetypes
        self.sorted_filetypes = sorted(six.iteritems(filetypes))
        default = 0
        for i, (ext, name) in enumerate(self.sorted_filetypes):
            cbox.append_text("%s (*.%s)" % (name, ext))
            if ext == default_filetype:
                default = i
        cbox.set_active(default)
        self.ext = default_filetype

        def cb_cbox_changed (cbox, data=None):
            """File extension changed"""
            head, filename = os.path.split(self.get_filename())
            root, ext = os.path.splitext(filename)
            ext = ext[1:]
            new_ext = self.sorted_filetypes[cbox.get_active()][0]
            self.ext = new_ext

            if ext in self.filetypes:
                filename = root + '.' + new_ext
            elif ext == '':
                filename = filename.rstrip('.') + '.' + new_ext

            self.set_current_name(filename)
        cbox.connect("changed", cb_cbox_changed)

        hbox.show_all()
        self.set_extra_widget(hbox)

    def get_filename_from_user (self):
        while True:
            filename = None
            if self.run() != int(gtk.RESPONSE_OK):
                break
            filename = self.get_filename()
            break

        return filename, self.ext


class DialogLineprops(object):
    """
    A GUI dialog for controlling lineprops
    """
    signals = (
        'on_combobox_lineprops_changed',
        'on_combobox_linestyle_changed',
        'on_combobox_marker_changed',
        'on_colorbutton_linestyle_color_set',
        'on_colorbutton_markerface_color_set',
        'on_dialog_lineprops_okbutton_clicked',
        'on_dialog_lineprops_cancelbutton_clicked',
        )

    linestyles = [ls for ls in lines.Line2D.lineStyles if ls.strip()]
    linestyled = {s: i for i, s in enumerate(linestyles)}

    markers =  [m for m in markers.MarkerStyle.markers
                if isinstance(m, six.string_types)]
    markerd = {s: i for i, s in enumerate(markers)}

    def __init__(self, lines):
        import gtk.glade

        datadir = matplotlib.get_data_path()
        gladefile = os.path.join(datadir, 'lineprops.glade')
        if not os.path.exists(gladefile):
            raise IOError(
                'Could not find gladefile lineprops.glade in %s' % datadir)

        self._inited = False
        self._updateson = True # suppress updates when setting widgets manually
        self.wtree = gtk.glade.XML(gladefile, 'dialog_lineprops')
        self.wtree.signal_autoconnect(
            {s: getattr(self, s) for s in self.signals})

        self.dlg = self.wtree.get_widget('dialog_lineprops')

        self.lines = lines

        cbox = self.wtree.get_widget('combobox_lineprops')
        cbox.set_active(0)
        self.cbox_lineprops = cbox

        cbox = self.wtree.get_widget('combobox_linestyles')
        for ls in self.linestyles:
            cbox.append_text(ls)
        cbox.set_active(0)
        self.cbox_linestyles = cbox

        cbox = self.wtree.get_widget('combobox_markers')
        for m in self.markers:
            cbox.append_text(m)
        cbox.set_active(0)
        self.cbox_markers = cbox
        self._lastcnt = 0
        self._inited = True

    def show(self):
        'populate the combo box'
        self._updateson = False
        # flush the old
        cbox = self.cbox_lineprops
        for i in range(self._lastcnt-1,-1,-1):
            cbox.remove_text(i)

        # add the new
        for line in self.lines:
            cbox.append_text(line.get_label())
        cbox.set_active(0)

        self._updateson = True
        self._lastcnt = len(self.lines)
        self.dlg.show()

    def get_active_line(self):
        'get the active line'
        ind = self.cbox_lineprops.get_active()
        line = self.lines[ind]
        return line

    def get_active_linestyle(self):
        'get the active lineinestyle'
        ind = self.cbox_linestyles.get_active()
        ls = self.linestyles[ind]
        return ls

    def get_active_marker(self):
        'get the active lineinestyle'
        ind = self.cbox_markers.get_active()
        m = self.markers[ind]
        return m

    def _update(self):
        'update the active line props from the widgets'
        if not self._inited or not self._updateson: return
        line = self.get_active_line()
        ls = self.get_active_linestyle()
        marker = self.get_active_marker()
        line.set_linestyle(ls)
        line.set_marker(marker)

        button = self.wtree.get_widget('colorbutton_linestyle')
        color = button.get_color()
        r, g, b = [val/65535. for val in (color.red, color.green, color.blue)]
        line.set_color((r,g,b))

        button = self.wtree.get_widget('colorbutton_markerface')
        color = button.get_color()
        r, g, b = [val/65535. for val in (color.red, color.green, color.blue)]
        line.set_markerfacecolor((r,g,b))

        line.figure.canvas.draw()

    def on_combobox_lineprops_changed(self, item):
        'update the widgets from the active line'
        if not self._inited: return
        self._updateson = False
        line = self.get_active_line()

        ls = line.get_linestyle()
        if ls is None: ls = 'None'
        self.cbox_linestyles.set_active(self.linestyled[ls])

        marker = line.get_marker()
        if marker is None: marker = 'None'
        self.cbox_markers.set_active(self.markerd[marker])

        rgba = mcolors.to_rgba(line.get_color())
        color = gtk.gdk.Color(*[int(val*65535) for val in rgba[:3]])
        button = self.wtree.get_widget('colorbutton_linestyle')
        button.set_color(color)

        rgba = mcolors.to_rgba(line.get_markerfacecolor())
        color = gtk.gdk.Color(*[int(val*65535) for val in rgba[:3]])
        button = self.wtree.get_widget('colorbutton_markerface')
        button.set_color(color)
        self._updateson = True

    def on_combobox_linestyle_changed(self, item):
        self._update()

    def on_combobox_marker_changed(self, item):
        self._update()

    def on_colorbutton_linestyle_color_set(self, button):
        self._update()

    def on_colorbutton_markerface_color_set(self, button):
        'called colorbutton marker clicked'
        self._update()

    def on_dialog_lineprops_okbutton_clicked(self, button):
        self._update()
        self.dlg.hide()

    def on_dialog_lineprops_cancelbutton_clicked(self, button):
        self.dlg.hide()

# set icon used when windows are minimized
# Unfortunately, the SVG renderer (rsvg) leaks memory under earlier
# versions of pygtk, so we have to use a PNG file instead.
try:
    if gtk.pygtk_version < (2, 8, 0) or sys.platform == 'win32':
        icon_filename = 'matplotlib.png'
    else:
        icon_filename = 'matplotlib.svg'
    window_icon = os.path.join(rcParams['datapath'], 'images', icon_filename)
except:
    window_icon = None
    _log.info('Could not load matplotlib icon: %s', sys.exc_info()[1])

def error_msg_gtk(msg, parent=None):
    if parent is not None: # find the toplevel gtk.Window
        parent = parent.get_toplevel()
        if parent.flags() & gtk.TOPLEVEL == 0:
            parent = None

    if not isinstance(msg, six.string_types):
        msg = ','.join(map(str, msg))

    dialog = gtk.MessageDialog(
        parent         = parent,
        type           = gtk.MESSAGE_ERROR,
        buttons        = gtk.BUTTONS_OK,
        message_format = msg)
    dialog.run()
    dialog.destroy()


@_Backend.export
class _BackendGTK(_Backend):
    FigureCanvas = FigureCanvasGTK
    FigureManager = FigureManagerGTK

    @staticmethod
    def trigger_manager_draw(manager):
        manager.canvas.draw_idle()

    @staticmethod
    def mainloop():
        if gtk.main_level() == 0:
            gtk.main()
