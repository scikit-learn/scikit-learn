from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import six

import functools
import os
import re
import signal
import sys
from six import unichr
import traceback

import matplotlib

from matplotlib._pylab_helpers import Gcf
from matplotlib.backend_bases import (
    _Backend, FigureCanvasBase, FigureManagerBase, NavigationToolbar2,
    TimerBase, cursors, ToolContainerBase, StatusbarBase)
import matplotlib.backends.qt_editor.figureoptions as figureoptions
from matplotlib.backends.qt_editor.formsubplottool import UiSubplotTool
from matplotlib.figure import Figure
from matplotlib.backend_managers import ToolManager
from matplotlib import backend_tools

from .qt_compat import (
    QtCore, QtGui, QtWidgets, _getSaveFileName, is_pyqt5, __version__, QT_API)

backend_version = __version__

# SPECIAL_KEYS are keys that do *not* return their unicode name
# instead they have manually specified names
SPECIAL_KEYS = {QtCore.Qt.Key_Control: 'control',
                QtCore.Qt.Key_Shift: 'shift',
                QtCore.Qt.Key_Alt: 'alt',
                QtCore.Qt.Key_Meta: 'super',
                QtCore.Qt.Key_Return: 'enter',
                QtCore.Qt.Key_Left: 'left',
                QtCore.Qt.Key_Up: 'up',
                QtCore.Qt.Key_Right: 'right',
                QtCore.Qt.Key_Down: 'down',
                QtCore.Qt.Key_Escape: 'escape',
                QtCore.Qt.Key_F1: 'f1',
                QtCore.Qt.Key_F2: 'f2',
                QtCore.Qt.Key_F3: 'f3',
                QtCore.Qt.Key_F4: 'f4',
                QtCore.Qt.Key_F5: 'f5',
                QtCore.Qt.Key_F6: 'f6',
                QtCore.Qt.Key_F7: 'f7',
                QtCore.Qt.Key_F8: 'f8',
                QtCore.Qt.Key_F9: 'f9',
                QtCore.Qt.Key_F10: 'f10',
                QtCore.Qt.Key_F11: 'f11',
                QtCore.Qt.Key_F12: 'f12',
                QtCore.Qt.Key_Home: 'home',
                QtCore.Qt.Key_End: 'end',
                QtCore.Qt.Key_PageUp: 'pageup',
                QtCore.Qt.Key_PageDown: 'pagedown',
                QtCore.Qt.Key_Tab: 'tab',
                QtCore.Qt.Key_Backspace: 'backspace',
                QtCore.Qt.Key_Enter: 'enter',
                QtCore.Qt.Key_Insert: 'insert',
                QtCore.Qt.Key_Delete: 'delete',
                QtCore.Qt.Key_Pause: 'pause',
                QtCore.Qt.Key_SysReq: 'sysreq',
                QtCore.Qt.Key_Clear: 'clear', }

# define which modifier keys are collected on keyboard events.
# elements are (mpl names, Modifier Flag, Qt Key) tuples
SUPER = 0
ALT = 1
CTRL = 2
SHIFT = 3
MODIFIER_KEYS = [('super', QtCore.Qt.MetaModifier, QtCore.Qt.Key_Meta),
                 ('alt', QtCore.Qt.AltModifier, QtCore.Qt.Key_Alt),
                 ('ctrl', QtCore.Qt.ControlModifier, QtCore.Qt.Key_Control),
                 ('shift', QtCore.Qt.ShiftModifier, QtCore.Qt.Key_Shift),
                 ]

if sys.platform == 'darwin':
    # in OSX, the control and super (aka cmd/apple) keys are switched, so
    # switch them back.
    SPECIAL_KEYS.update({QtCore.Qt.Key_Control: 'cmd',  # cmd/apple key
                         QtCore.Qt.Key_Meta: 'control',
                         })
    MODIFIER_KEYS[0] = ('cmd', QtCore.Qt.ControlModifier,
                        QtCore.Qt.Key_Control)
    MODIFIER_KEYS[2] = ('ctrl', QtCore.Qt.MetaModifier,
                        QtCore.Qt.Key_Meta)


cursord = {
    cursors.MOVE: QtCore.Qt.SizeAllCursor,
    cursors.HAND: QtCore.Qt.PointingHandCursor,
    cursors.POINTER: QtCore.Qt.ArrowCursor,
    cursors.SELECT_REGION: QtCore.Qt.CrossCursor,
    cursors.WAIT: QtCore.Qt.WaitCursor,
    }


# make place holder
qApp = None


def _create_qApp():
    """
    Only one qApp can exist at a time, so check before creating one.
    """
    global qApp

    if qApp is None:
        app = QtWidgets.QApplication.instance()
        if app is None:
            # check for DISPLAY env variable on X11 build of Qt
            if is_pyqt5():
                try:
                    from PyQt5 import QtX11Extras
                    is_x11_build = True
                except ImportError:
                    is_x11_build = False
            else:
                is_x11_build = hasattr(QtGui, "QX11Info")
            if is_x11_build:
                display = os.environ.get('DISPLAY')
                if display is None or not re.search(r':\d', display):
                    raise RuntimeError('Invalid DISPLAY variable')

            qApp = QtWidgets.QApplication([b"matplotlib"])
            qApp.lastWindowClosed.connect(qApp.quit)
        else:
            qApp = app

    if is_pyqt5():
        try:
            qApp.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps)
            qApp.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
        except AttributeError:
            pass


def _allow_super_init(__init__):
    """
    Decorator for ``__init__`` to allow ``super().__init__`` on PyQt4/PySide2.
    """

    if QT_API == "PyQt5":

        return __init__

    else:
        # To work around lack of cooperative inheritance in PyQt4, PySide,
        # and PySide2, when calling FigureCanvasQT.__init__, we temporarily
        # patch QWidget.__init__ by a cooperative version, that first calls
        # QWidget.__init__ with no additional arguments, and then finds the
        # next class in the MRO with an __init__ that does support cooperative
        # inheritance (i.e., not defined by the PyQt4, PySide, PySide2, sip
        # or Shiboken packages), and manually call its `__init__`, once again
        # passing the additional arguments.

        qwidget_init = QtWidgets.QWidget.__init__

        def cooperative_qwidget_init(self, *args, **kwargs):
            qwidget_init(self)
            mro = type(self).__mro__
            next_coop_init = next(
                cls for cls in mro[mro.index(QtWidgets.QWidget) + 1:]
                if cls.__module__.split(".")[0] not in [
                    "PyQt4", "sip", "PySide", "PySide2", "Shiboken"])
            next_coop_init.__init__(self, *args, **kwargs)

        @functools.wraps(__init__)
        def wrapper(self, **kwargs):
            try:
                QtWidgets.QWidget.__init__ = cooperative_qwidget_init
                __init__(self, **kwargs)
            finally:
                # Restore __init__
                QtWidgets.QWidget.__init__ = qwidget_init

        return wrapper


class TimerQT(TimerBase):
    '''
    Subclass of :class:`backend_bases.TimerBase` that uses Qt timer events.

    Attributes
    ----------
    interval : int
        The time between timer events in milliseconds. Default is 1000 ms.
    single_shot : bool
        Boolean flag indicating whether this timer should
        operate as single shot (run once and then stop). Defaults to False.
    callbacks : list
        Stores list of (func, args) tuples that will be called upon timer
        events. This list can be manipulated directly, or the functions
        `add_callback` and `remove_callback` can be used.

    '''

    def __init__(self, *args, **kwargs):
        TimerBase.__init__(self, *args, **kwargs)

        # Create a new timer and connect the timeout() signal to the
        # _on_timer method.
        self._timer = QtCore.QTimer()
        self._timer.timeout.connect(self._on_timer)
        self._timer_set_interval()

    def _timer_set_single_shot(self):
        self._timer.setSingleShot(self._single)

    def _timer_set_interval(self):
        self._timer.setInterval(self._interval)

    def _timer_start(self):
        self._timer.start()

    def _timer_stop(self):
        self._timer.stop()


class FigureCanvasQT(QtWidgets.QWidget, FigureCanvasBase):

    # map Qt button codes to MouseEvent's ones:
    buttond = {QtCore.Qt.LeftButton: 1,
               QtCore.Qt.MidButton: 2,
               QtCore.Qt.RightButton: 3,
               # QtCore.Qt.XButton1: None,
               # QtCore.Qt.XButton2: None,
               }

    @_allow_super_init
    def __init__(self, figure):
        _create_qApp()
        super(FigureCanvasQT, self).__init__(figure=figure)

        self.figure = figure
        # We don't want to scale up the figure DPI more than once.
        # Note, we don't handle a signal for changing DPI yet.
        figure._original_dpi = figure.dpi
        self._update_figure_dpi()
        # In cases with mixed resolution displays, we need to be careful if the
        # dpi_ratio changes - in this case we need to resize the canvas
        # accordingly. We could watch for screenChanged events from Qt, but
        # the issue is that we can't guarantee this will be emitted *before*
        # the first paintEvent for the canvas, so instead we keep track of the
        # dpi_ratio value here and in paintEvent we resize the canvas if
        # needed.
        self._dpi_ratio_prev = None

        self._draw_pending = False
        self._is_drawing = False
        self._draw_rect_callback = lambda painter: None

        self.setAttribute(QtCore.Qt.WA_OpaquePaintEvent)
        self.setMouseTracking(True)
        self.resize(*self.get_width_height())
        # Key auto-repeat enabled by default
        self._keyautorepeat = True

        palette = QtGui.QPalette(QtCore.Qt.white)
        self.setPalette(palette)

    def _update_figure_dpi(self):
        dpi = self._dpi_ratio * self.figure._original_dpi
        self.figure._set_dpi(dpi, forward=False)

    @property
    def _dpi_ratio(self):
        # Not available on Qt4 or some older Qt5.
        try:
            # self.devicePixelRatio() returns 0 in rare cases
            return self.devicePixelRatio() or 1
        except AttributeError:
            return 1

    def _update_dpi(self):
        # As described in __init__ above, we need to be careful in cases with
        # mixed resolution displays if dpi_ratio is changing between painting
        # events.
        # Return whether we triggered a resizeEvent (and thus a paintEvent)
        # from within this function.
        if self._dpi_ratio != self._dpi_ratio_prev:
            # We need to update the figure DPI.
            self._update_figure_dpi()
            self._dpi_ratio_prev = self._dpi_ratio
            # The easiest way to resize the canvas is to emit a resizeEvent
            # since we implement all the logic for resizing the canvas for
            # that event.
            event = QtGui.QResizeEvent(self.size(), self.size())
            self.resizeEvent(event)
            # resizeEvent triggers a paintEvent itself, so we exit this one
            # (after making sure that the event is immediately handled).
            return True
        return False

    def get_width_height(self):
        w, h = FigureCanvasBase.get_width_height(self)
        return int(w / self._dpi_ratio), int(h / self._dpi_ratio)

    def enterEvent(self, event):
        FigureCanvasBase.enter_notify_event(self, guiEvent=event)

    def leaveEvent(self, event):
        QtWidgets.QApplication.restoreOverrideCursor()
        FigureCanvasBase.leave_notify_event(self, guiEvent=event)

    def mouseEventCoords(self, pos):
        """Calculate mouse coordinates in physical pixels

        Qt5 use logical pixels, but the figure is scaled to physical
        pixels for rendering.   Transform to physical pixels so that
        all of the down-stream transforms work as expected.

        Also, the origin is different and needs to be corrected.

        """
        dpi_ratio = self._dpi_ratio
        x = pos.x()
        # flip y so y=0 is bottom of canvas
        y = self.figure.bbox.height / dpi_ratio - pos.y()
        return x * dpi_ratio, y * dpi_ratio

    def mousePressEvent(self, event):
        x, y = self.mouseEventCoords(event.pos())
        button = self.buttond.get(event.button())
        if button is not None:
            FigureCanvasBase.button_press_event(self, x, y, button,
                                                guiEvent=event)

    def mouseDoubleClickEvent(self, event):
        x, y = self.mouseEventCoords(event.pos())
        button = self.buttond.get(event.button())
        if button is not None:
            FigureCanvasBase.button_press_event(self, x, y,
                                                button, dblclick=True,
                                                guiEvent=event)

    def mouseMoveEvent(self, event):
        x, y = self.mouseEventCoords(event)
        FigureCanvasBase.motion_notify_event(self, x, y, guiEvent=event)

    def mouseReleaseEvent(self, event):
        x, y = self.mouseEventCoords(event)
        button = self.buttond.get(event.button())
        if button is not None:
            FigureCanvasBase.button_release_event(self, x, y, button,
                                                  guiEvent=event)

    if is_pyqt5():
        def wheelEvent(self, event):
            x, y = self.mouseEventCoords(event)
            # from QWheelEvent::delta doc
            if event.pixelDelta().x() == 0 and event.pixelDelta().y() == 0:
                steps = event.angleDelta().y() / 120
            else:
                steps = event.pixelDelta().y()
            if steps:
                FigureCanvasBase.scroll_event(
                    self, x, y, steps, guiEvent=event)
    else:
        def wheelEvent(self, event):
            x = event.x()
            # flipy so y=0 is bottom of canvas
            y = self.figure.bbox.height - event.y()
            # from QWheelEvent::delta doc
            steps = event.delta() / 120
            if event.orientation() == QtCore.Qt.Vertical:
                FigureCanvasBase.scroll_event(
                    self, x, y, steps, guiEvent=event)

    def keyPressEvent(self, event):
        key = self._get_key(event)
        if key is not None:
            FigureCanvasBase.key_press_event(self, key, guiEvent=event)

    def keyReleaseEvent(self, event):
        key = self._get_key(event)
        if key is not None:
            FigureCanvasBase.key_release_event(self, key, guiEvent=event)

    @property
    def keyAutoRepeat(self):
        """
        If True, enable auto-repeat for key events.
        """
        return self._keyautorepeat

    @keyAutoRepeat.setter
    def keyAutoRepeat(self, val):
        self._keyautorepeat = bool(val)

    def resizeEvent(self, event):
        # _dpi_ratio_prev will be set the first time the canvas is painted, and
        # the rendered buffer is useless before anyways.
        if self._dpi_ratio_prev is None:
            return
        w = event.size().width() * self._dpi_ratio
        h = event.size().height() * self._dpi_ratio
        dpival = self.figure.dpi
        winch = w / dpival
        hinch = h / dpival
        self.figure.set_size_inches(winch, hinch, forward=False)
        # pass back into Qt to let it finish
        QtWidgets.QWidget.resizeEvent(self, event)
        # emit our resize events
        FigureCanvasBase.resize_event(self)

    def sizeHint(self):
        w, h = self.get_width_height()
        return QtCore.QSize(w, h)

    def minumumSizeHint(self):
        return QtCore.QSize(10, 10)

    def _get_key(self, event):
        if not self._keyautorepeat and event.isAutoRepeat():
            return None

        event_key = event.key()
        event_mods = int(event.modifiers())  # actually a bitmask

        # get names of the pressed modifier keys
        # bit twiddling to pick out modifier keys from event_mods bitmask,
        # if event_key is a MODIFIER, it should not be duplicated in mods
        mods = [name for name, mod_key, qt_key in MODIFIER_KEYS
                if event_key != qt_key and (event_mods & mod_key) == mod_key]
        try:
            # for certain keys (enter, left, backspace, etc) use a word for the
            # key, rather than unicode
            key = SPECIAL_KEYS[event_key]
        except KeyError:
            # unicode defines code points up to 0x0010ffff
            # QT will use Key_Codes larger than that for keyboard keys that are
            # are not unicode characters (like multimedia keys)
            # skip these
            # if you really want them, you should add them to SPECIAL_KEYS
            MAX_UNICODE = 0x10ffff
            if event_key > MAX_UNICODE:
                return None

            key = unichr(event_key)
            # qt delivers capitalized letters.  fix capitalization
            # note that capslock is ignored
            if 'shift' in mods:
                mods.remove('shift')
            else:
                key = key.lower()

        mods.reverse()
        return '+'.join(mods + [key])

    def new_timer(self, *args, **kwargs):
        """
        Creates a new backend-specific subclass of
        :class:`backend_bases.Timer`.  This is useful for getting
        periodic events through the backend's native event
        loop. Implemented only for backends with GUIs.

        Other Parameters
        ----------------
        interval : scalar
            Timer interval in milliseconds

        callbacks : list
            Sequence of (func, args, kwargs) where ``func(*args, **kwargs)``
            will be executed by the timer every *interval*.

        """
        return TimerQT(*args, **kwargs)

    def flush_events(self):
        qApp.processEvents()

    def start_event_loop(self, timeout=0):
        if hasattr(self, "_event_loop") and self._event_loop.isRunning():
            raise RuntimeError("Event loop already running")
        self._event_loop = event_loop = QtCore.QEventLoop()
        if timeout:
            timer = QtCore.QTimer.singleShot(timeout * 1000, event_loop.quit)
        event_loop.exec_()

    def stop_event_loop(self, event=None):
        if hasattr(self, "_event_loop"):
            self._event_loop.quit()

    def draw(self):
        """Render the figure, and queue a request for a Qt draw.
        """
        # The renderer draw is done here; delaying causes problems with code
        # that uses the result of the draw() to update plot elements.
        if self._is_drawing:
            return
        self._is_drawing = True
        try:
            super(FigureCanvasQT, self).draw()
        finally:
            self._is_drawing = False
        self.update()

    def draw_idle(self):
        """Queue redraw of the Agg buffer and request Qt paintEvent.
        """
        # The Agg draw needs to be handled by the same thread matplotlib
        # modifies the scene graph from. Post Agg draw request to the
        # current event loop in order to ensure thread affinity and to
        # accumulate multiple draw requests from event handling.
        # TODO: queued signal connection might be safer than singleShot
        if not (self._draw_pending or self._is_drawing):
            self._draw_pending = True
            QtCore.QTimer.singleShot(0, self._draw_idle)

    def _draw_idle(self):
        if self.height() < 0 or self.width() < 0:
            self._draw_pending = False
        if not self._draw_pending:
            return
        try:
            self.draw()
        except Exception:
            # Uncaught exceptions are fatal for PyQt5, so catch them instead.
            traceback.print_exc()
        finally:
            self._draw_pending = False

    def drawRectangle(self, rect):
        # Draw the zoom rectangle to the QPainter.  _draw_rect_callback needs
        # to be called at the end of paintEvent.
        if rect is not None:
            def _draw_rect_callback(painter):
                pen = QtGui.QPen(QtCore.Qt.black, 1 / self._dpi_ratio,
                                 QtCore.Qt.DotLine)
                painter.setPen(pen)
                painter.drawRect(*(pt / self._dpi_ratio for pt in rect))
        else:
            def _draw_rect_callback(painter):
                return
        self._draw_rect_callback = _draw_rect_callback
        self.update()


class MainWindow(QtWidgets.QMainWindow):
    closing = QtCore.Signal()

    def closeEvent(self, event):
        self.closing.emit()
        QtWidgets.QMainWindow.closeEvent(self, event)


class FigureManagerQT(FigureManagerBase):
    """
    Attributes
    ----------
    canvas : `FigureCanvas`
        The FigureCanvas instance
    num : int or str
        The Figure number
    toolbar : qt.QToolBar
        The qt.QToolBar
    window : qt.QMainWindow
        The qt.QMainWindow

    """

    def __init__(self, canvas, num):
        FigureManagerBase.__init__(self, canvas, num)
        self.canvas = canvas
        self.window = MainWindow()
        self.window.closing.connect(canvas.close_event)
        self.window.closing.connect(self._widgetclosed)

        self.window.setWindowTitle("Figure %d" % num)
        image = os.path.join(matplotlib.rcParams['datapath'],
                             'images', 'matplotlib.svg')
        self.window.setWindowIcon(QtGui.QIcon(image))

        # Give the keyboard focus to the figure instead of the
        # manager; StrongFocus accepts both tab and click to focus and
        # will enable the canvas to process event w/o clicking.
        # ClickFocus only takes the focus is the window has been
        # clicked
        # on. http://qt-project.org/doc/qt-4.8/qt.html#FocusPolicy-enum or
        # http://doc.qt.digia.com/qt/qt.html#FocusPolicy-enum
        self.canvas.setFocusPolicy(QtCore.Qt.StrongFocus)
        self.canvas.setFocus()

        self.window._destroying = False

        self.toolmanager = self._get_toolmanager()
        self.toolbar = self._get_toolbar(self.canvas, self.window)
        self.statusbar = None

        if self.toolmanager:
            backend_tools.add_tools_to_manager(self.toolmanager)
            if self.toolbar:
                backend_tools.add_tools_to_container(self.toolbar)
                self.statusbar = StatusbarQt(self.window, self.toolmanager)

        if self.toolbar is not None:
            self.window.addToolBar(self.toolbar)
            if not self.toolmanager:
                # add text label to status bar
                statusbar_label = QtWidgets.QLabel()
                self.window.statusBar().addWidget(statusbar_label)
                self.toolbar.message.connect(statusbar_label.setText)
            tbs_height = self.toolbar.sizeHint().height()
        else:
            tbs_height = 0

        # resize the main window so it will display the canvas with the
        # requested size:
        cs = canvas.sizeHint()
        sbs = self.window.statusBar().sizeHint()
        self._status_and_tool_height = tbs_height + sbs.height()
        height = cs.height() + self._status_and_tool_height
        self.window.resize(cs.width(), height)

        self.window.setCentralWidget(self.canvas)

        if matplotlib.is_interactive():
            self.window.show()
            self.canvas.draw_idle()

        def notify_axes_change(fig):
            # This will be called whenever the current axes is changed
            if self.toolbar is not None:
                self.toolbar.update()
        self.canvas.figure.add_axobserver(notify_axes_change)
        self.window.raise_()

    def full_screen_toggle(self):
        if self.window.isFullScreen():
            self.window.showNormal()
        else:
            self.window.showFullScreen()

    def _widgetclosed(self):
        if self.window._destroying:
            return
        self.window._destroying = True
        try:
            Gcf.destroy(self.num)
        except AttributeError:
            pass
            # It seems that when the python session is killed,
            # Gcf can get destroyed before the Gcf.destroy
            # line is run, leading to a useless AttributeError.

    def _get_toolbar(self, canvas, parent):
        # must be inited after the window, drawingArea and figure
        # attrs are set
        if matplotlib.rcParams['toolbar'] == 'toolbar2':
            toolbar = NavigationToolbar2QT(canvas, parent, False)
        elif matplotlib.rcParams['toolbar'] == 'toolmanager':
            toolbar = ToolbarQt(self.toolmanager, self.window)
        else:
            toolbar = None
        return toolbar

    def _get_toolmanager(self):
        if matplotlib.rcParams['toolbar'] == 'toolmanager':
            toolmanager = ToolManager(self.canvas.figure)
        else:
            toolmanager = None
        return toolmanager

    def resize(self, width, height):
        'set the canvas size in pixels'
        self.window.resize(width, height + self._status_and_tool_height)

    def show(self):
        self.window.show()
        self.window.activateWindow()
        self.window.raise_()

    def destroy(self, *args):
        # check for qApp first, as PySide deletes it in its atexit handler
        if QtWidgets.QApplication.instance() is None:
            return
        if self.window._destroying:
            return
        self.window._destroying = True
        self.window.destroyed.connect(self._widgetclosed)
        if self.toolbar:
            self.toolbar.destroy()
        self.window.close()

    def get_window_title(self):
        return six.text_type(self.window.windowTitle())

    def set_window_title(self, title):
        self.window.setWindowTitle(title)


class NavigationToolbar2QT(NavigationToolbar2, QtWidgets.QToolBar):
    message = QtCore.Signal(str)

    def __init__(self, canvas, parent, coordinates=True):
        """ coordinates: should we show the coordinates on the right? """
        self.canvas = canvas
        self.parent = parent
        self.coordinates = coordinates
        self._actions = {}
        """A mapping of toolitem method names to their QActions"""

        QtWidgets.QToolBar.__init__(self, parent)
        NavigationToolbar2.__init__(self, canvas)

    def _icon(self, name):
        if is_pyqt5():
            name = name.replace('.png', '_large.png')
        pm = QtGui.QPixmap(os.path.join(self.basedir, name))
        if hasattr(pm, 'setDevicePixelRatio'):
            pm.setDevicePixelRatio(self.canvas._dpi_ratio)
        return QtGui.QIcon(pm)

    def _init_toolbar(self):
        self.basedir = os.path.join(matplotlib.rcParams['datapath'], 'images')

        for text, tooltip_text, image_file, callback in self.toolitems:
            if text is None:
                self.addSeparator()
            else:
                a = self.addAction(self._icon(image_file + '.png'),
                                   text, getattr(self, callback))
                self._actions[callback] = a
                if callback in ['zoom', 'pan']:
                    a.setCheckable(True)
                if tooltip_text is not None:
                    a.setToolTip(tooltip_text)
                if text == 'Subplots':
                    a = self.addAction(self._icon("qt4_editor_options.png"),
                                       'Customize', self.edit_parameters)
                    a.setToolTip('Edit axis, curve and image parameters')

        self.buttons = {}

        # Add the x,y location widget at the right side of the toolbar
        # The stretch factor is 1 which means any resizing of the toolbar
        # will resize this label instead of the buttons.
        if self.coordinates:
            self.locLabel = QtWidgets.QLabel("", self)
            self.locLabel.setAlignment(
                    QtCore.Qt.AlignRight | QtCore.Qt.AlignTop)
            self.locLabel.setSizePolicy(
                QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                      QtWidgets.QSizePolicy.Ignored))
            labelAction = self.addWidget(self.locLabel)
            labelAction.setVisible(True)

        # reference holder for subplots_adjust window
        self.adj_window = None

        # Esthetic adjustments - we need to set these explicitly in PyQt5
        # otherwise the layout looks different - but we don't want to set it if
        # not using HiDPI icons otherwise they look worse than before.
        if is_pyqt5():
            self.setIconSize(QtCore.QSize(24, 24))
            self.layout().setSpacing(12)

    if is_pyqt5():
        # For some reason, self.setMinimumHeight doesn't seem to carry over to
        # the actual sizeHint, so override it instead in order to make the
        # aesthetic adjustments noted above.
        def sizeHint(self):
            size = super(NavigationToolbar2QT, self).sizeHint()
            size.setHeight(max(48, size.height()))
            return size

    def edit_parameters(self):
        allaxes = self.canvas.figure.get_axes()
        if not allaxes:
            QtWidgets.QMessageBox.warning(
                self.parent, "Error", "There are no axes to edit.")
            return
        elif len(allaxes) == 1:
            axes, = allaxes
        else:
            titles = []
            for axes in allaxes:
                name = (axes.get_title() or
                        " - ".join(filter(None, [axes.get_xlabel(),
                                                 axes.get_ylabel()])) or
                        "<anonymous {} (id: {:#x})>".format(
                            type(axes).__name__, id(axes)))
                titles.append(name)
            item, ok = QtWidgets.QInputDialog.getItem(
                self.parent, 'Customize', 'Select axes:', titles, 0, False)
            if ok:
                axes = allaxes[titles.index(six.text_type(item))]
            else:
                return

        figureoptions.figure_edit(axes, self)

    def _update_buttons_checked(self):
        # sync button checkstates to match active mode
        self._actions['pan'].setChecked(self._active == 'PAN')
        self._actions['zoom'].setChecked(self._active == 'ZOOM')

    def pan(self, *args):
        super(NavigationToolbar2QT, self).pan(*args)
        self._update_buttons_checked()

    def zoom(self, *args):
        super(NavigationToolbar2QT, self).zoom(*args)
        self._update_buttons_checked()

    def set_message(self, s):
        self.message.emit(s)
        if self.coordinates:
            self.locLabel.setText(s)

    def set_cursor(self, cursor):
        self.canvas.setCursor(cursord[cursor])

    def draw_rubberband(self, event, x0, y0, x1, y1):
        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0
        rect = [int(val) for val in (x0, y0, x1 - x0, y1 - y0)]
        self.canvas.drawRectangle(rect)

    def remove_rubberband(self):
        self.canvas.drawRectangle(None)

    def configure_subplots(self):
        image = os.path.join(matplotlib.rcParams['datapath'],
                             'images', 'matplotlib.png')
        dia = SubplotToolQt(self.canvas.figure, self.parent)
        dia.setWindowIcon(QtGui.QIcon(image))
        dia.exec_()

    def save_figure(self, *args):
        filetypes = self.canvas.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(six.iteritems(filetypes))
        default_filetype = self.canvas.get_default_filetype()

        startpath = os.path.expanduser(
            matplotlib.rcParams['savefig.directory'])
        start = os.path.join(startpath, self.canvas.get_default_filename())
        filters = []
        selectedFilter = None
        for name, exts in sorted_filetypes:
            exts_list = " ".join(['*.%s' % ext for ext in exts])
            filter = '%s (%s)' % (name, exts_list)
            if default_filetype in exts:
                selectedFilter = filter
            filters.append(filter)
        filters = ';;'.join(filters)

        fname, filter = _getSaveFileName(self.parent,
                                         "Choose a filename to save to",
                                         start, filters, selectedFilter)
        if fname:
            # Save dir for next time, unless empty str (i.e., use cwd).
            if startpath != "":
                matplotlib.rcParams['savefig.directory'] = (
                    os.path.dirname(six.text_type(fname)))
            try:
                self.canvas.figure.savefig(six.text_type(fname))
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self, "Error saving file", six.text_type(e),
                    QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)


class SubplotToolQt(UiSubplotTool):
    def __init__(self, targetfig, parent):
        UiSubplotTool.__init__(self, None)

        self._figure = targetfig

        for lower, higher in [("bottom", "top"), ("left", "right")]:
            self._widgets[lower].valueChanged.connect(
                lambda val: self._widgets[higher].setMinimum(val + .001))
            self._widgets[higher].valueChanged.connect(
                lambda val: self._widgets[lower].setMaximum(val - .001))

        self._attrs = ["top", "bottom", "left", "right", "hspace", "wspace"]
        self._defaults = {attr: vars(self._figure.subplotpars)[attr]
                          for attr in self._attrs}

        # Set values after setting the range callbacks, but before setting up
        # the redraw callbacks.
        self._reset()

        for attr in self._attrs:
            self._widgets[attr].valueChanged.connect(self._on_value_changed)
        for action, method in [("Export values", self._export_values),
                               ("Tight layout", self._tight_layout),
                               ("Reset", self._reset),
                               ("Close", self.close)]:
            self._widgets[action].clicked.connect(method)

    def _export_values(self):
        # Explicitly round to 3 decimals (which is also the spinbox precision)
        # to avoid numbers of the form 0.100...001.
        dialog = QtWidgets.QDialog()
        layout = QtWidgets.QVBoxLayout()
        dialog.setLayout(layout)
        text = QtWidgets.QPlainTextEdit()
        text.setReadOnly(True)
        layout.addWidget(text)
        text.setPlainText(
            ",\n".join("{}={:.3}".format(attr, self._widgets[attr].value())
                       for attr in self._attrs))
        # Adjust the height of the text widget to fit the whole text, plus
        # some padding.
        size = text.maximumSize()
        size.setHeight(
            QtGui.QFontMetrics(text.document().defaultFont())
            .size(0, text.toPlainText()).height() + 20)
        text.setMaximumSize(size)
        dialog.exec_()

    def _on_value_changed(self):
        self._figure.subplots_adjust(**{attr: self._widgets[attr].value()
                                        for attr in self._attrs})
        self._figure.canvas.draw_idle()

    def _tight_layout(self):
        self._figure.tight_layout()
        for attr in self._attrs:
            widget = self._widgets[attr]
            widget.blockSignals(True)
            widget.setValue(vars(self._figure.subplotpars)[attr])
            widget.blockSignals(False)
        self._figure.canvas.draw_idle()

    def _reset(self):
        for attr, value in self._defaults.items():
            self._widgets[attr].setValue(value)


class ToolbarQt(ToolContainerBase, QtWidgets.QToolBar):
    def __init__(self, toolmanager, parent):
        ToolContainerBase.__init__(self, toolmanager)
        QtWidgets.QToolBar.__init__(self, parent)
        self._toolitems = {}
        self._groups = {}
        self._last = None

    @property
    def _icon_extension(self):
        if is_pyqt5():
            return '_large.png'
        return '.png'

    def add_toolitem(
            self, name, group, position, image_file, description, toggle):

        button = QtWidgets.QToolButton(self)
        button.setIcon(self._icon(image_file))
        button.setText(name)
        if description:
            button.setToolTip(description)

        def handler():
            self.trigger_tool(name)
        if toggle:
            button.setCheckable(True)
            button.toggled.connect(handler)
        else:
            button.clicked.connect(handler)

        self._last = button
        self._toolitems.setdefault(name, [])
        self._add_to_group(group, name, button, position)
        self._toolitems[name].append((button, handler))

    def _add_to_group(self, group, name, button, position):
        gr = self._groups.get(group, [])
        if not gr:
            sep = self.addSeparator()
            gr.append(sep)
        before = gr[position]
        widget = self.insertWidget(before, button)
        gr.insert(position, widget)
        self._groups[group] = gr

    def _icon(self, name):
        pm = QtGui.QPixmap(name)
        if hasattr(pm, 'setDevicePixelRatio'):
            pm.setDevicePixelRatio(self.toolmanager.canvas._dpi_ratio)
        return QtGui.QIcon(pm)

    def toggle_toolitem(self, name, toggled):
        if name not in self._toolitems:
            return
        for button, handler in self._toolitems[name]:
            button.toggled.disconnect(handler)
            button.setChecked(toggled)
            button.toggled.connect(handler)

    def remove_toolitem(self, name):
        for button, handler in self._toolitems[name]:
            button.setParent(None)
        del self._toolitems[name]


class StatusbarQt(StatusbarBase, QtWidgets.QLabel):
    def __init__(self, window, *args, **kwargs):
        StatusbarBase.__init__(self, *args, **kwargs)
        QtWidgets.QLabel.__init__(self)
        window.statusBar().addWidget(self)

    def set_message(self, s):
        self.setText(s)


class ConfigureSubplotsQt(backend_tools.ConfigureSubplotsBase):
    def trigger(self, *args):
        image = os.path.join(matplotlib.rcParams['datapath'],
                             'images', 'matplotlib.png')
        parent = self.canvas.manager.window
        dia = SubplotToolQt(self.figure, parent)
        dia.setWindowIcon(QtGui.QIcon(image))
        dia.exec_()


class SaveFigureQt(backend_tools.SaveFigureBase):
    def trigger(self, *args):
        filetypes = self.canvas.get_supported_filetypes_grouped()
        sorted_filetypes = sorted(six.iteritems(filetypes))
        default_filetype = self.canvas.get_default_filetype()

        startpath = os.path.expanduser(
            matplotlib.rcParams['savefig.directory'])
        start = os.path.join(startpath, self.canvas.get_default_filename())
        filters = []
        selectedFilter = None
        for name, exts in sorted_filetypes:
            exts_list = " ".join(['*.%s' % ext for ext in exts])
            filter = '%s (%s)' % (name, exts_list)
            if default_filetype in exts:
                selectedFilter = filter
            filters.append(filter)
        filters = ';;'.join(filters)

        parent = self.canvas.manager.window
        fname, filter = _getSaveFileName(parent,
                                         "Choose a filename to save to",
                                         start, filters, selectedFilter)
        if fname:
            # Save dir for next time, unless empty str (i.e., use cwd).
            if startpath != "":
                matplotlib.rcParams['savefig.directory'] = (
                    os.path.dirname(six.text_type(fname)))
            try:
                self.canvas.figure.savefig(six.text_type(fname))
            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    self, "Error saving file", six.text_type(e),
                    QtWidgets.QMessageBox.Ok, QtWidgets.QMessageBox.NoButton)


class SetCursorQt(backend_tools.SetCursorBase):
    def set_cursor(self, cursor):
        self.canvas.setCursor(cursord[cursor])


class RubberbandQt(backend_tools.RubberbandBase):
    def draw_rubberband(self, x0, y0, x1, y1):
        height = self.canvas.figure.bbox.height
        y1 = height - y1
        y0 = height - y0
        rect = [int(val) for val in (x0, y0, x1 - x0, y1 - y0)]
        self.canvas.drawRectangle(rect)

    def remove_rubberband(self):
        self.canvas.drawRectangle(None)


backend_tools.ToolSaveFigure = SaveFigureQt
backend_tools.ToolConfigureSubplots = ConfigureSubplotsQt
backend_tools.ToolSetCursor = SetCursorQt
backend_tools.ToolRubberband = RubberbandQt


def error_msg_qt(msg, parent=None):
    if not isinstance(msg, six.string_types):
        msg = ','.join(map(str, msg))

    QtWidgets.QMessageBox.warning(None, "Matplotlib",
                                  msg, QtGui.QMessageBox.Ok)


def exception_handler(type, value, tb):
    """Handle uncaught exceptions
    It does not catch SystemExit
    """
    msg = ''
    # get the filename attribute if available (for IOError)
    if hasattr(value, 'filename') and value.filename is not None:
        msg = value.filename + ': '
    if hasattr(value, 'strerror') and value.strerror is not None:
        msg += value.strerror
    else:
        msg += six.text_type(value)

    if len(msg):
        error_msg_qt(msg)


@_Backend.export
class _BackendQT5(_Backend):
    FigureCanvas = FigureCanvasQT
    FigureManager = FigureManagerQT

    @staticmethod
    def trigger_manager_draw(manager):
        manager.canvas.draw_idle()

    @staticmethod
    def mainloop():
        # allow KeyboardInterrupt exceptions to close the plot window.
        signal.signal(signal.SIGINT, signal.SIG_DFL)
        qApp.exec_()
