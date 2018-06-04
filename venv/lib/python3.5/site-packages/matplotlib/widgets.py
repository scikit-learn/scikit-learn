"""
GUI neutral widgets
===================

Widgets that are designed to work for any of the GUI backends.
All of these widgets require you to predefine a :class:`matplotlib.axes.Axes`
instance and pass that as the first arg.  matplotlib doesn't try to
be too smart with respect to layout -- you will have to figure out how
wide and tall you want your Axes to be to accommodate your widget.
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import copy
import six
from six.moves import zip

import numpy as np
from matplotlib import rcParams

from .patches import Circle, Rectangle, Ellipse
from .lines import Line2D
from .transforms import blended_transform_factory


class LockDraw(object):
    """
    Some widgets, like the cursor, draw onto the canvas, and this is not
    desirable under all circumstances, like when the toolbar is in
    zoom-to-rect mode and drawing a rectangle.  The module level "lock"
    allows someone to grab the lock and prevent other widgets from
    drawing.  Use ``matplotlib.widgets.lock(someobj)`` to prevent
    other widgets from drawing while you're interacting with the canvas.
    """

    def __init__(self):
        self._owner = None

    def __call__(self, o):
        """reserve the lock for *o*"""
        if not self.available(o):
            raise ValueError('already locked')
        self._owner = o

    def release(self, o):
        """release the lock"""
        if not self.available(o):
            raise ValueError('you do not own this lock')
        self._owner = None

    def available(self, o):
        """drawing is available to *o*"""
        return not self.locked() or self.isowner(o)

    def isowner(self, o):
        """Return True if *o* owns this lock"""
        return self._owner is o

    def locked(self):
        """Return True if the lock is currently held by an owner"""
        return self._owner is not None


class Widget(object):
    """
    Abstract base class for GUI neutral widgets
    """
    drawon = True
    eventson = True
    _active = True

    def set_active(self, active):
        """Set whether the widget is active.
        """
        self._active = active

    def get_active(self):
        """Get whether the widget is active.
        """
        return self._active

    # set_active is overridden by SelectorWidgets.
    active = property(get_active, lambda self, active: self.set_active(active),
                      doc="Is the widget active?")

    def ignore(self, event):
        """Return True if event should be ignored.

        This method (or a version of it) should be called at the beginning
        of any event callback.
        """
        return not self.active


class AxesWidget(Widget):
    """Widget that is connected to a single
    :class:`~matplotlib.axes.Axes`.

    To guarantee that the widget remains responsive and not garbage-collected,
    a reference to the object should be maintained by the user.

    This is necessary because the callback registry
    maintains only weak-refs to the functions, which are member
    functions of the widget.  If there are no references to the widget
    object it may be garbage collected which will disconnect the
    callbacks.

    Attributes:

    *ax* : :class:`~matplotlib.axes.Axes`
        The parent axes for the widget
    *canvas* : :class:`~matplotlib.backend_bases.FigureCanvasBase` subclass
        The parent figure canvas for the widget.
    *active* : bool
        If False, the widget does not respond to events.
    """
    def __init__(self, ax):
        self.ax = ax
        self.canvas = ax.figure.canvas
        self.cids = []

    def connect_event(self, event, callback):
        """Connect callback with an event.

        This should be used in lieu of `figure.canvas.mpl_connect` since this
        function stores callback ids for later clean up.
        """
        cid = self.canvas.mpl_connect(event, callback)
        self.cids.append(cid)

    def disconnect_events(self):
        """Disconnect all events created by this widget."""
        for c in self.cids:
            self.canvas.mpl_disconnect(c)


class Button(AxesWidget):
    """
    A GUI neutral button.

    For the button to remain responsive you must keep a reference to it.
    Call :meth:`on_clicked` to connect to the button.

    Attributes
    ----------
    ax :
        The :class:`matplotlib.axes.Axes` the button renders into.
    label :
        A :class:`matplotlib.text.Text` instance.
    color :
        The color of the button when not hovering.
    hovercolor :
        The color of the button when hovering.
    """

    def __init__(self, ax, label, image=None,
                 color='0.85', hovercolor='0.95'):
        """
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The :class:`matplotlib.axes.Axes` instance the button
            will be placed into.

        label : str
            The button text. Accepts string.

        image : array, mpl image, Pillow Image
            The image to place in the button, if not *None*.
            Can be any legal arg to imshow (numpy array,
            matplotlib Image instance, or Pillow Image).

        color : color
            The color of the button when not activated

        hovercolor : color
            The color of the button when the mouse is over it
        """
        AxesWidget.__init__(self, ax)

        if image is not None:
            ax.imshow(image)
        self.label = ax.text(0.5, 0.5, label,
                             verticalalignment='center',
                             horizontalalignment='center',
                             transform=ax.transAxes)

        self.cnt = 0
        self.observers = {}

        self.connect_event('button_press_event', self._click)
        self.connect_event('button_release_event', self._release)
        self.connect_event('motion_notify_event', self._motion)
        ax.set_navigate(False)
        ax.set_facecolor(color)
        ax.set_xticks([])
        ax.set_yticks([])
        self.color = color
        self.hovercolor = hovercolor

        self._lastcolor = color

    def _click(self, event):
        if self.ignore(event):
            return
        if event.inaxes != self.ax:
            return
        if not self.eventson:
            return
        if event.canvas.mouse_grabber != self.ax:
            event.canvas.grab_mouse(self.ax)

    def _release(self, event):
        if self.ignore(event):
            return
        if event.canvas.mouse_grabber != self.ax:
            return
        event.canvas.release_mouse(self.ax)
        if not self.eventson:
            return
        if event.inaxes != self.ax:
            return
        for cid, func in six.iteritems(self.observers):
            func(event)

    def _motion(self, event):
        if self.ignore(event):
            return
        if event.inaxes == self.ax:
            c = self.hovercolor
        else:
            c = self.color
        if c != self._lastcolor:
            self.ax.set_facecolor(c)
            self._lastcolor = c
            if self.drawon:
                self.ax.figure.canvas.draw()

    def on_clicked(self, func):
        """
        When the button is clicked, call this *func* with event.

        A connection id is returned. It can be used to disconnect
        the button from its callback.
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass


class Slider(AxesWidget):
    """
    A slider representing a floating point range.

    Create a slider from *valmin* to *valmax* in axes *ax*. For the slider to
    remain responsive you must maintain a reference to it. Call
    :meth:`on_changed` to connect to the slider event.

    Attributes
    ----------
    val : float
        Slider value.
    """
    def __init__(self, ax, label, valmin, valmax, valinit=0.5, valfmt='%1.2f',
                 closedmin=True, closedmax=True, slidermin=None,
                 slidermax=None, dragging=True, valstep=None, **kwargs):
        """
        Parameters
        ----------
        ax : Axes
            The Axes to put the slider in.

        label : str
            Slider label.

        valmin : float
            The minimum value of the slider.

        valmax : float
            The maximum value of the slider.

        valinit : float, optional, default: 0.5
            The slider initial position.

        valfmt : str, optional, default: "%1.2f"
            Used to format the slider value, fprint format string.

        closedmin : bool, optional, default: True
            Indicate whether the slider interval is closed on the bottom.

        closedmax : bool, optional, default: True
            Indicate whether the slider interval is closed on the top.

        slidermin : Slider, optional, default: None
            Do not allow the current slider to have a value less than
            the value of the Slider `slidermin`.

        slidermax : Slider, optional, default: None
            Do not allow the current slider to have a value greater than
            the value of the Slider `slidermax`.

        dragging : bool, optional, default: True
            If True the slider can be dragged by the mouse.

        valstep : float, optional, default: None
            If given, the slider will snap to multiples of `valstep`.

        Notes
        -----
        Additional kwargs are passed on to ``self.poly`` which is the
        :class:`~matplotlib.patches.Rectangle` that draws the slider
        knob.  See the :class:`~matplotlib.patches.Rectangle` documentation for
        valid property names (e.g., `facecolor`, `edgecolor`, `alpha`).
        """
        AxesWidget.__init__(self, ax)

        if slidermin is not None and not hasattr(slidermin, 'val'):
            raise ValueError("Argument slidermin ({}) has no 'val'"
                             .format(type(slidermin)))
        if slidermax is not None and not hasattr(slidermax, 'val'):
            raise ValueError("Argument slidermax ({}) has no 'val'"
                             .format(type(slidermax)))
        self.closedmin = closedmin
        self.closedmax = closedmax
        self.slidermin = slidermin
        self.slidermax = slidermax
        self.drag_active = False
        self.valmin = valmin
        self.valmax = valmax
        self.valstep = valstep
        valinit = self._value_in_bounds(valinit)
        if valinit is None:
            valinit = valmin
        self.val = valinit
        self.valinit = valinit
        self.poly = ax.axvspan(valmin, valinit, 0, 1, **kwargs)
        self.vline = ax.axvline(valinit, 0, 1, color='r', lw=1)

        self.valfmt = valfmt
        ax.set_yticks([])
        ax.set_xlim((valmin, valmax))
        ax.set_xticks([])
        ax.set_navigate(False)

        self.connect_event('button_press_event', self._update)
        self.connect_event('button_release_event', self._update)
        if dragging:
            self.connect_event('motion_notify_event', self._update)
        self.label = ax.text(-0.02, 0.5, label, transform=ax.transAxes,
                             verticalalignment='center',
                             horizontalalignment='right')

        self.valtext = ax.text(1.02, 0.5, valfmt % valinit,
                               transform=ax.transAxes,
                               verticalalignment='center',
                               horizontalalignment='left')

        self.cnt = 0
        self.observers = {}

        self.set_val(valinit)

    def _value_in_bounds(self, val):
        """ Makes sure self.val is with given bounds."""
        if self.valstep:
            val = np.round((val - self.valmin)/self.valstep)*self.valstep
            val += self.valmin

        if val <= self.valmin:
            if not self.closedmin:
                return
            val = self.valmin
        elif val >= self.valmax:
            if not self.closedmax:
                return
            val = self.valmax

        if self.slidermin is not None and val <= self.slidermin.val:
            if not self.closedmin:
                return
            val = self.slidermin.val

        if self.slidermax is not None and val >= self.slidermax.val:
            if not self.closedmax:
                return
            val = self.slidermax.val
        return val

    def _update(self, event):
        """update the slider position"""
        if self.ignore(event):
            return

        if event.button != 1:
            return

        if event.name == 'button_press_event' and event.inaxes == self.ax:
            self.drag_active = True
            event.canvas.grab_mouse(self.ax)

        if not self.drag_active:
            return

        elif ((event.name == 'button_release_event') or
              (event.name == 'button_press_event' and
               event.inaxes != self.ax)):
            self.drag_active = False
            event.canvas.release_mouse(self.ax)
            return
        val = self._value_in_bounds(event.xdata)
        if (val is not None) and (val != self.val):
            self.set_val(val)

    def set_val(self, val):
        """
        Set slider value to *val*

        Parameters
        ----------
        val : float
        """
        xy = self.poly.xy
        xy[2] = val, 1
        xy[3] = val, 0
        self.poly.xy = xy
        self.valtext.set_text(self.valfmt % val)
        if self.drawon:
            self.ax.figure.canvas.draw_idle()
        self.val = val
        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(val)

    def on_changed(self, func):
        """
        When the slider value is changed call *func* with the new
        slider value

        Parameters
        ----------
        func : callable
            Function to call when slider is changed.
            The function must accept a single float as its arguments.

        Returns
        -------
        cid : int
            Connection id (which can be used to disconnect *func*)
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """
        Remove the observer with connection id *cid*

        Parameters
        ----------
        cid : int
            Connection id of the observer to be removed
        """
        try:
            del self.observers[cid]
        except KeyError:
            pass

    def reset(self):
        """Reset the slider to the initial value"""
        if (self.val != self.valinit):
            self.set_val(self.valinit)


class CheckButtons(AxesWidget):
    """
    A GUI neutral set of check buttons.

    For the check buttons to remain responsive you must keep a
    reference to this object.

    The following attributes are exposed

     *ax*
        The :class:`matplotlib.axes.Axes` instance the buttons are
        located in

     *labels*
        List of :class:`matplotlib.text.Text` instances

     *lines*
        List of (line1, line2) tuples for the x's in the check boxes.
        These lines exist for each box, but have ``set_visible(False)``
        when its box is not checked.

     *rectangles*
        List of :class:`matplotlib.patches.Rectangle` instances

    Connect to the CheckButtons with the :meth:`on_clicked` method
    """
    def __init__(self, ax, labels, actives):
        """
        Add check buttons to :class:`matplotlib.axes.Axes` instance *ax*

        *labels*
            A len(buttons) list of labels as strings

        *actives*
            A len(buttons) list of booleans indicating whether
             the button is active
        """
        AxesWidget.__init__(self, ax)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_navigate(False)

        if len(labels) > 1:
            dy = 1. / (len(labels) + 1)
            ys = np.linspace(1 - dy, dy, len(labels))
        else:
            dy = 0.25
            ys = [0.5]

        cnt = 0
        axcolor = ax.get_facecolor()

        self.labels = []
        self.lines = []
        self.rectangles = []

        lineparams = {'color': 'k', 'linewidth': 1.25,
                      'transform': ax.transAxes, 'solid_capstyle': 'butt'}
        for y, label in zip(ys, labels):
            t = ax.text(0.25, y, label, transform=ax.transAxes,
                        horizontalalignment='left',
                        verticalalignment='center')

            w, h = dy / 2., dy / 2.
            x, y = 0.05, y - h / 2.

            p = Rectangle(xy=(x, y), width=w, height=h, edgecolor='black',
                          facecolor=axcolor, transform=ax.transAxes)

            l1 = Line2D([x, x + w], [y + h, y], **lineparams)
            l2 = Line2D([x, x + w], [y, y + h], **lineparams)

            l1.set_visible(actives[cnt])
            l2.set_visible(actives[cnt])
            self.labels.append(t)
            self.rectangles.append(p)
            self.lines.append((l1, l2))
            ax.add_patch(p)
            ax.add_line(l1)
            ax.add_line(l2)
            cnt += 1

        self.connect_event('button_press_event', self._clicked)

        self.cnt = 0
        self.observers = {}

    def _clicked(self, event):
        if self.ignore(event):
            return
        if event.button != 1:
            return
        if event.inaxes != self.ax:
            return

        for i, (p, t) in enumerate(zip(self.rectangles, self.labels)):
            if (t.get_window_extent().contains(event.x, event.y) or
                    p.get_window_extent().contains(event.x, event.y)):
                self.set_active(i)
                break
        else:
            return

    def set_active(self, index):
        """
        Directly (de)activate a check button by index.

        *index* is an index into the original label list
            that this object was constructed with.
            Raises ValueError if *index* is invalid.

        Callbacks will be triggered if :attr:`eventson` is True.

        """
        if 0 > index >= len(self.labels):
            raise ValueError("Invalid CheckButton index: %d" % index)

        l1, l2 = self.lines[index]
        l1.set_visible(not l1.get_visible())
        l2.set_visible(not l2.get_visible())

        if self.drawon:
            self.ax.figure.canvas.draw()

        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(self.labels[index].get_text())

    def get_status(self):
        """
        returns a tuple of the status (True/False) of all of the check buttons
        """
        return [l1.get_visible() for (l1, l2) in self.lines]

    def on_clicked(self, func):
        """
        When the button is clicked, call *func* with button label

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass


class TextBox(AxesWidget):
    """
    A GUI neutral text input box.

    For the text box to remain responsive you must keep a reference to it.

    The following attributes are accessible:

      *ax*
        The :class:`matplotlib.axes.Axes` the button renders into.

      *label*
        A :class:`matplotlib.text.Text` instance.

      *color*
        The color of the text box when not hovering.

      *hovercolor*
        The color of the text box when hovering.

    Call :meth:`on_text_change` to be updated whenever the text changes.

    Call :meth:`on_submit` to be updated whenever the user hits enter or
    leaves the text entry field.
    """

    def __init__(self, ax, label, initial='',
                 color='.95', hovercolor='1', label_pad=.01):
        """
        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The :class:`matplotlib.axes.Axes` instance the button
            will be placed into.

        label : str
            Label for this text box. Accepts string.

        initial : str
            Initial value in the text box

        color : color
            The color of the box

        hovercolor : color
            The color of the box when the mouse is over it

        label_pad : float
            the distance between the label and the right side of the textbox
        """
        AxesWidget.__init__(self, ax)

        self.DIST_FROM_LEFT = .05

        self.params_to_disable = [key for key in rcParams if u'keymap' in key]

        self.text = initial
        self.label = ax.text(-label_pad, 0.5, label,
                             verticalalignment='center',
                             horizontalalignment='right',
                             transform=ax.transAxes)
        self.text_disp = self._make_text_disp(self.text)

        self.cnt = 0
        self.change_observers = {}
        self.submit_observers = {}

        # If these lines are removed, the cursor won't appear the first
        # time the box is clicked:
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)

        self.cursor_index = 0

        # Because this is initialized, _render_cursor
        # can assume that cursor exists.
        self.cursor = self.ax.vlines(0, 0, 0)
        self.cursor.set_visible(False)

        self.connect_event('button_press_event', self._click)
        self.connect_event('button_release_event', self._release)
        self.connect_event('motion_notify_event', self._motion)
        self.connect_event('key_press_event', self._keypress)
        self.connect_event('resize_event', self._resize)
        ax.set_navigate(False)
        ax.set_facecolor(color)
        ax.set_xticks([])
        ax.set_yticks([])
        self.color = color
        self.hovercolor = hovercolor

        self._lastcolor = color

        self.capturekeystrokes = False

    def _make_text_disp(self, string):
        return self.ax.text(self.DIST_FROM_LEFT, 0.5, string,
                            verticalalignment='center',
                            horizontalalignment='left',
                            transform=self.ax.transAxes)

    def _rendercursor(self):
        # this is a hack to figure out where the cursor should go.
        # we draw the text up to where the cursor should go, measure
        # and save its dimensions, draw the real text, then put the cursor
        # at the saved dimensions

        widthtext = self.text[:self.cursor_index]
        no_text = False
        if(widthtext == "" or widthtext == " " or widthtext == "  "):
            no_text = widthtext == ""
            widthtext = ","

        wt_disp = self._make_text_disp(widthtext)

        self.ax.figure.canvas.draw()
        bb = wt_disp.get_window_extent()
        inv = self.ax.transData.inverted()
        bb = inv.transform(bb)
        wt_disp.set_visible(False)
        if no_text:
            bb[1, 0] = bb[0, 0]
        # hack done
        self.cursor.set_visible(False)

        self.cursor = self.ax.vlines(bb[1, 0], bb[0, 1], bb[1, 1])
        self.ax.figure.canvas.draw()

    def _notify_submit_observers(self):
        for cid, func in six.iteritems(self.submit_observers):
                func(self.text)

    def _release(self, event):
        if self.ignore(event):
            return
        if event.canvas.mouse_grabber != self.ax:
            return
        event.canvas.release_mouse(self.ax)

    def _keypress(self, event):
        if self.ignore(event):
            return
        if self.capturekeystrokes:
            key = event.key

            if(len(key) == 1):
                self.text = (self.text[:self.cursor_index] + key +
                             self.text[self.cursor_index:])
                self.cursor_index += 1
            elif key == "right":
                if self.cursor_index != len(self.text):
                    self.cursor_index += 1
            elif key == "left":
                if self.cursor_index != 0:
                    self.cursor_index -= 1
            elif key == "home":
                self.cursor_index = 0
            elif key == "end":
                self.cursor_index = len(self.text)
            elif(key == "backspace"):
                if self.cursor_index != 0:
                    self.text = (self.text[:self.cursor_index - 1] +
                                 self.text[self.cursor_index:])
                    self.cursor_index -= 1
            elif(key == "delete"):
                if self.cursor_index != len(self.text):
                    self.text = (self.text[:self.cursor_index] +
                                 self.text[self.cursor_index + 1:])

            self.text_disp.remove()
            self.text_disp = self._make_text_disp(self.text)
            self._rendercursor()
            self._notify_change_observers()
            if key == "enter":
                self._notify_submit_observers()

    def set_val(self, val):
        newval = str(val)
        if self.text == newval:
            return
        self.text = newval
        self.text_disp.remove()
        self.text_disp = self._make_text_disp(self.text)
        self._rendercursor()
        self._notify_change_observers()
        self._notify_submit_observers()

    def _notify_change_observers(self):
        for cid, func in six.iteritems(self.change_observers):
            func(self.text)

    def begin_typing(self, x):
        self.capturekeystrokes = True
        # disable command keys so that the user can type without
        # command keys causing figure to be saved, etc
        self.reset_params = {}
        for key in self.params_to_disable:
            self.reset_params[key] = rcParams[key]
            rcParams[key] = []

    def stop_typing(self):
        notifysubmit = False
        # because _notify_submit_users might throw an error in the
        # user's code, we only want to call it once we've already done
        # our cleanup.
        if self.capturekeystrokes:
            # since the user is no longer typing,
            # reactivate the standard command keys
            for key in self.params_to_disable:
                rcParams[key] = self.reset_params[key]
            notifysubmit = True
        self.capturekeystrokes = False
        self.cursor.set_visible(False)
        self.ax.figure.canvas.draw()
        if notifysubmit:
            self._notify_submit_observers()

    def position_cursor(self, x):
        # now, we have to figure out where the cursor goes.
        # approximate it based on assuming all characters the same length
        if len(self.text) == 0:
            self.cursor_index = 0
        else:
            bb = self.text_disp.get_window_extent()

            trans = self.ax.transData
            inv = self.ax.transData.inverted()
            bb = trans.transform(inv.transform(bb))

            text_start = bb[0, 0]
            text_end = bb[1, 0]

            ratio = (x - text_start) / (text_end - text_start)

            if ratio < 0:
                ratio = 0
            if ratio > 1:
                ratio = 1

            self.cursor_index = int(len(self.text) * ratio)

        self._rendercursor()

    def _click(self, event):
        if self.ignore(event):
            return
        if event.inaxes != self.ax:
            self.stop_typing()
            return
        if not self.eventson:
            return
        if event.canvas.mouse_grabber != self.ax:
            event.canvas.grab_mouse(self.ax)
        if not self.capturekeystrokes:
            self.begin_typing(event.x)
        self.position_cursor(event.x)

    def _resize(self, event):
        self.stop_typing()

    def _motion(self, event):
        if self.ignore(event):
            return
        if event.inaxes == self.ax:
            c = self.hovercolor
        else:
            c = self.color
        if c != self._lastcolor:
            self.ax.set_facecolor(c)
            self._lastcolor = c
            if self.drawon:
                self.ax.figure.canvas.draw()

    def on_text_change(self, func):
        """
        When the text changes, call this *func* with event.

        A connection id is returned which can be used to disconnect.
        """
        cid = self.cnt
        self.change_observers[cid] = func
        self.cnt += 1
        return cid

    def on_submit(self, func):
        """
        When the user hits enter or leaves the submision box, call this
        *func* with event.

        A connection id is returned which can be used to disconnect.
        """
        cid = self.cnt
        self.submit_observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        for reg in (self.change_observers, self.submit_observers):
            try:
                del reg[cid]
            except KeyError:
                pass


class RadioButtons(AxesWidget):
    """
    A GUI neutral radio button.

    For the buttons to remain responsive
    you must keep a reference to this object.

    The following attributes are exposed:

     *ax*
        The :class:`matplotlib.axes.Axes` instance the buttons are in

     *activecolor*
        The color of the button when clicked

     *labels*
        A list of :class:`matplotlib.text.Text` instances

     *circles*
        A list of :class:`matplotlib.patches.Circle` instances

     *value_selected*
        A string listing the current value selected

    Connect to the RadioButtons with the :meth:`on_clicked` method
    """
    def __init__(self, ax, labels, active=0, activecolor='blue'):
        """
        Add radio buttons to :class:`matplotlib.axes.Axes` instance *ax*

        *labels*
            A len(buttons) list of labels as strings

        *active*
            The index into labels for the button that is active

        *activecolor*
            The color of the button when clicked
        """
        AxesWidget.__init__(self, ax)
        self.activecolor = activecolor
        self.value_selected = None

        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_navigate(False)
        dy = 1. / (len(labels) + 1)
        ys = np.linspace(1 - dy, dy, len(labels))
        cnt = 0
        axcolor = ax.get_facecolor()

        self.labels = []
        self.circles = []
        for y, label in zip(ys, labels):
            t = ax.text(0.25, y, label, transform=ax.transAxes,
                        horizontalalignment='left',
                        verticalalignment='center')

            if cnt == active:
                self.value_selected = label
                facecolor = activecolor
            else:
                facecolor = axcolor

            p = Circle(xy=(0.15, y), radius=0.05, edgecolor='black',
                       facecolor=facecolor, transform=ax.transAxes)

            self.labels.append(t)
            self.circles.append(p)
            ax.add_patch(p)
            cnt += 1

        self.connect_event('button_press_event', self._clicked)

        self.cnt = 0
        self.observers = {}

    def _clicked(self, event):
        if self.ignore(event):
            return
        if event.button != 1:
            return
        if event.inaxes != self.ax:
            return
        xy = self.ax.transAxes.inverted().transform_point((event.x, event.y))
        pclicked = np.array([xy[0], xy[1]])

        def inside(p):
            pcirc = np.array([p.center[0], p.center[1]])
            d = pclicked - pcirc
            return np.sqrt(np.dot(d, d)) < p.radius

        for i, (p, t) in enumerate(zip(self.circles, self.labels)):
            if t.get_window_extent().contains(event.x, event.y) or inside(p):
                self.set_active(i)
                break
        else:
            return

    def set_active(self, index):
        """
        Trigger which radio button to make active.

        *index* is an index into the original label list
            that this object was constructed with.
            Raise ValueError if the index is invalid.

        Callbacks will be triggered if :attr:`eventson` is True.

        """
        if 0 > index >= len(self.labels):
            raise ValueError("Invalid RadioButton index: %d" % index)

        self.value_selected = self.labels[index].get_text()

        for i, p in enumerate(self.circles):
            if i == index:
                color = self.activecolor
            else:
                color = self.ax.get_facecolor()
            p.set_facecolor(color)

        if self.drawon:
            self.ax.figure.canvas.draw()

        if not self.eventson:
            return
        for cid, func in six.iteritems(self.observers):
            func(self.labels[index].get_text())

    def on_clicked(self, func):
        """
        When the button is clicked, call *func* with button label

        A connection id is returned which can be used to disconnect
        """
        cid = self.cnt
        self.observers[cid] = func
        self.cnt += 1
        return cid

    def disconnect(self, cid):
        """remove the observer with connection id *cid*"""
        try:
            del self.observers[cid]
        except KeyError:
            pass


class SubplotTool(Widget):
    """
    A tool to adjust the subplot params of a :class:`matplotlib.figure.Figure`.
    """
    def __init__(self, targetfig, toolfig):
        """
        *targetfig*
            The figure instance to adjust.

        *toolfig*
            The figure instance to embed the subplot tool into. If
            *None*, a default figure will be created. If you are using
            this from the GUI
        """
        # FIXME: The docstring seems to just abruptly end without...

        self.targetfig = targetfig
        toolfig.subplots_adjust(left=0.2, right=0.9)

        class toolbarfmt:
            def __init__(self, slider):
                self.slider = slider

            def __call__(self, x, y):
                fmt = '%s=%s' % (self.slider.label.get_text(),
                                 self.slider.valfmt)
                return fmt % x

        self.axleft = toolfig.add_subplot(711)
        self.axleft.set_title('Click on slider to adjust subplot param')
        self.axleft.set_navigate(False)

        self.sliderleft = Slider(self.axleft, 'left',
                                 0, 1, targetfig.subplotpars.left,
                                 closedmax=False)
        self.sliderleft.on_changed(self.funcleft)

        self.axbottom = toolfig.add_subplot(712)
        self.axbottom.set_navigate(False)
        self.sliderbottom = Slider(self.axbottom,
                                   'bottom', 0, 1,
                                   targetfig.subplotpars.bottom,
                                   closedmax=False)
        self.sliderbottom.on_changed(self.funcbottom)

        self.axright = toolfig.add_subplot(713)
        self.axright.set_navigate(False)
        self.sliderright = Slider(self.axright, 'right', 0, 1,
                                  targetfig.subplotpars.right,
                                  closedmin=False)
        self.sliderright.on_changed(self.funcright)

        self.axtop = toolfig.add_subplot(714)
        self.axtop.set_navigate(False)
        self.slidertop = Slider(self.axtop, 'top', 0, 1,
                                targetfig.subplotpars.top,
                                closedmin=False)
        self.slidertop.on_changed(self.functop)

        self.axwspace = toolfig.add_subplot(715)
        self.axwspace.set_navigate(False)
        self.sliderwspace = Slider(self.axwspace, 'wspace',
                                   0, 1, targetfig.subplotpars.wspace,
                                   closedmax=False)
        self.sliderwspace.on_changed(self.funcwspace)

        self.axhspace = toolfig.add_subplot(716)
        self.axhspace.set_navigate(False)
        self.sliderhspace = Slider(self.axhspace, 'hspace',
                                   0, 1, targetfig.subplotpars.hspace,
                                   closedmax=False)
        self.sliderhspace.on_changed(self.funchspace)

        # constraints
        self.sliderleft.slidermax = self.sliderright
        self.sliderright.slidermin = self.sliderleft
        self.sliderbottom.slidermax = self.slidertop
        self.slidertop.slidermin = self.sliderbottom

        bax = toolfig.add_axes([0.8, 0.05, 0.15, 0.075])
        self.buttonreset = Button(bax, 'Reset')

        sliders = (self.sliderleft, self.sliderbottom, self.sliderright,
                   self.slidertop, self.sliderwspace, self.sliderhspace,)

        def func(event):
            thisdrawon = self.drawon

            self.drawon = False

            # store the drawon state of each slider
            bs = []
            for slider in sliders:
                bs.append(slider.drawon)
                slider.drawon = False

            # reset the slider to the initial position
            for slider in sliders:
                slider.reset()

            # reset drawon
            for slider, b in zip(sliders, bs):
                slider.drawon = b

            # draw the canvas
            self.drawon = thisdrawon
            if self.drawon:
                toolfig.canvas.draw()
                self.targetfig.canvas.draw()

        # during reset there can be a temporary invalid state
        # depending on the order of the reset so we turn off
        # validation for the resetting
        validate = toolfig.subplotpars.validate
        toolfig.subplotpars.validate = False
        self.buttonreset.on_clicked(func)
        toolfig.subplotpars.validate = validate

    def funcleft(self, val):
        self.targetfig.subplots_adjust(left=val)
        if self.drawon:
            self.targetfig.canvas.draw()

    def funcright(self, val):
        self.targetfig.subplots_adjust(right=val)
        if self.drawon:
            self.targetfig.canvas.draw()

    def funcbottom(self, val):
        self.targetfig.subplots_adjust(bottom=val)
        if self.drawon:
            self.targetfig.canvas.draw()

    def functop(self, val):
        self.targetfig.subplots_adjust(top=val)
        if self.drawon:
            self.targetfig.canvas.draw()

    def funcwspace(self, val):
        self.targetfig.subplots_adjust(wspace=val)
        if self.drawon:
            self.targetfig.canvas.draw()

    def funchspace(self, val):
        self.targetfig.subplots_adjust(hspace=val)
        if self.drawon:
            self.targetfig.canvas.draw()


class Cursor(AxesWidget):
    """
    A horizontal and vertical line that spans the axes and moves with
    the pointer.  You can turn off the hline or vline respectively with
    the following attributes:

      *horizOn*
        Controls the visibility of the horizontal line

      *vertOn*
        Controls the visibility of the horizontal line

    and the visibility of the cursor itself with the *visible* attribute.

    For the cursor to remain responsive you must keep a reference to
    it.
    """
    def __init__(self, ax, horizOn=True, vertOn=True, useblit=False,
                 **lineprops):
        """
        Add a cursor to *ax*.  If ``useblit=True``, use the backend-dependent
        blitting features for faster updates.  *lineprops* is a dictionary of
        line properties.
        """
        AxesWidget.__init__(self, ax)

        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('draw_event', self.clear)

        self.visible = True
        self.horizOn = horizOn
        self.vertOn = vertOn
        self.useblit = useblit and self.canvas.supports_blit

        if self.useblit:
            lineprops['animated'] = True
        self.lineh = ax.axhline(ax.get_ybound()[0], visible=False, **lineprops)
        self.linev = ax.axvline(ax.get_xbound()[0], visible=False, **lineprops)

        self.background = None
        self.needclear = False

    def clear(self, event):
        """clear the cursor"""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.linev.set_visible(False)
        self.lineh.set_visible(False)

    def onmove(self, event):
        """on mouse motion draw the cursor if visible"""
        if self.ignore(event):
            return
        if not self.canvas.widgetlock.available(self):
            return
        if event.inaxes != self.ax:
            self.linev.set_visible(False)
            self.lineh.set_visible(False)

            if self.needclear:
                self.canvas.draw()
                self.needclear = False
            return
        self.needclear = True
        if not self.visible:
            return
        self.linev.set_xdata((event.xdata, event.xdata))

        self.lineh.set_ydata((event.ydata, event.ydata))
        self.linev.set_visible(self.visible and self.vertOn)
        self.lineh.set_visible(self.visible and self.horizOn)

        self._update()

    def _update(self):

        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.linev)
            self.ax.draw_artist(self.lineh)
            self.canvas.blit(self.ax.bbox)
        else:

            self.canvas.draw_idle()

        return False


class MultiCursor(Widget):
    """
    Provide a vertical (default) and/or horizontal line cursor shared between
    multiple axes.

    For the cursor to remain responsive you must keep a reference to
    it.

    Example usage::

        from matplotlib.widgets import MultiCursor
        from pylab import figure, show, np

        t = np.arange(0.0, 2.0, 0.01)
        s1 = np.sin(2*np.pi*t)
        s2 = np.sin(4*np.pi*t)
        fig = figure()
        ax1 = fig.add_subplot(211)
        ax1.plot(t, s1)


        ax2 = fig.add_subplot(212, sharex=ax1)
        ax2.plot(t, s2)

        multi = MultiCursor(fig.canvas, (ax1, ax2), color='r', lw=1,
                            horizOn=False, vertOn=True)
        show()

    """
    def __init__(self, canvas, axes, useblit=True, horizOn=False, vertOn=True,
                 **lineprops):

        self.canvas = canvas
        self.axes = axes
        self.horizOn = horizOn
        self.vertOn = vertOn

        xmin, xmax = axes[-1].get_xlim()
        ymin, ymax = axes[-1].get_ylim()
        xmid = 0.5 * (xmin + xmax)
        ymid = 0.5 * (ymin + ymax)

        self.visible = True
        self.useblit = useblit and self.canvas.supports_blit
        self.background = None
        self.needclear = False

        if self.useblit:
            lineprops['animated'] = True

        if vertOn:
            self.vlines = [ax.axvline(xmid, visible=False, **lineprops)
                           for ax in axes]
        else:
            self.vlines = []

        if horizOn:
            self.hlines = [ax.axhline(ymid, visible=False, **lineprops)
                           for ax in axes]
        else:
            self.hlines = []

        self.connect()

    def connect(self):
        """connect events"""
        self._cidmotion = self.canvas.mpl_connect('motion_notify_event',
                                                  self.onmove)
        self._ciddraw = self.canvas.mpl_connect('draw_event', self.clear)

    def disconnect(self):
        """disconnect events"""
        self.canvas.mpl_disconnect(self._cidmotion)
        self.canvas.mpl_disconnect(self._ciddraw)

    def clear(self, event):
        """clear the cursor"""
        if self.ignore(event):
            return
        if self.useblit:
            self.background = (
                self.canvas.copy_from_bbox(self.canvas.figure.bbox))
        for line in self.vlines + self.hlines:
            line.set_visible(False)

    def onmove(self, event):
        if self.ignore(event):
            return
        if event.inaxes is None:
            return
        if not self.canvas.widgetlock.available(self):
            return
        self.needclear = True
        if not self.visible:
            return
        if self.vertOn:
            for line in self.vlines:
                line.set_xdata((event.xdata, event.xdata))
                line.set_visible(self.visible)
        if self.horizOn:
            for line in self.hlines:
                line.set_ydata((event.ydata, event.ydata))
                line.set_visible(self.visible)
        self._update()

    def _update(self):
        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            if self.vertOn:
                for ax, line in zip(self.axes, self.vlines):
                    ax.draw_artist(line)
            if self.horizOn:
                for ax, line in zip(self.axes, self.hlines):
                    ax.draw_artist(line)
            self.canvas.blit(self.canvas.figure.bbox)
        else:
            self.canvas.draw_idle()


class _SelectorWidget(AxesWidget):

    def __init__(self, ax, onselect, useblit=False, button=None,
                 state_modifier_keys=None):
        AxesWidget.__init__(self, ax)

        self.visible = True
        self.onselect = onselect
        self.useblit = useblit and self.canvas.supports_blit
        self.connect_default_events()

        self.state_modifier_keys = dict(move=' ', clear='escape',
                                        square='shift', center='control')
        self.state_modifier_keys.update(state_modifier_keys or {})

        self.background = None
        self.artists = []

        if isinstance(button, int):
            self.validButtons = [button]
        else:
            self.validButtons = button

        # will save the data (position at mouseclick)
        self.eventpress = None
        # will save the data (pos. at mouserelease)
        self.eventrelease = None
        self._prev_event = None
        self.state = set()

    def set_active(self, active):
        AxesWidget.set_active(self, active)
        if active:
            self.update_background(None)

    def update_background(self, event):
        """force an update of the background"""
        # If you add a call to `ignore` here, you'll want to check edge case:
        # `release` can call a draw event even when `ignore` is True.
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)

    def connect_default_events(self):
        """Connect the major canvas events to methods."""
        self.connect_event('motion_notify_event', self.onmove)
        self.connect_event('button_press_event', self.press)
        self.connect_event('button_release_event', self.release)
        self.connect_event('draw_event', self.update_background)
        self.connect_event('key_press_event', self.on_key_press)
        self.connect_event('key_release_event', self.on_key_release)
        self.connect_event('scroll_event', self.on_scroll)

    def ignore(self, event):
        """return *True* if *event* should be ignored"""
        if not self.active or not self.ax.get_visible():
            return True

        # If canvas was locked
        if not self.canvas.widgetlock.available(self):
            return True

        if not hasattr(event, 'button'):
            event.button = None

        # Only do rectangle selection if event was triggered
        # with a desired button
        if self.validButtons is not None:
            if event.button not in self.validButtons:
                return True

        # If no button was pressed yet ignore the event if it was out
        # of the axes
        if self.eventpress is None:
            return event.inaxes != self.ax

        # If a button was pressed, check if the release-button is the
        # same.
        if event.button == self.eventpress.button:
            return False

        # If a button was pressed, check if the release-button is the
        # same.
        return (event.inaxes != self.ax or
                event.button != self.eventpress.button)

    def update(self):
        """draw using newfangled blit or oldfangled draw depending on
        useblit

        """
        if not self.ax.get_visible():
            return False

        if self.useblit:
            if self.background is not None:
                self.canvas.restore_region(self.background)
            for artist in self.artists:
                self.ax.draw_artist(artist)

            self.canvas.blit(self.ax.bbox)

        else:
            self.canvas.draw_idle()
        return False

    def _get_data(self, event):
        """Get the xdata and ydata for event, with limits"""
        if event.xdata is None:
            return None, None
        x0, x1 = self.ax.get_xbound()
        y0, y1 = self.ax.get_ybound()
        xdata = max(x0, event.xdata)
        xdata = min(x1, xdata)
        ydata = max(y0, event.ydata)
        ydata = min(y1, ydata)
        return xdata, ydata

    def _clean_event(self, event):
        """Clean up an event

        Use prev event if there is no xdata
        Limit the xdata and ydata to the axes limits
        Set the prev event
        """
        if event.xdata is None:
            event = self._prev_event
        else:
            event = copy.copy(event)
        event.xdata, event.ydata = self._get_data(event)

        self._prev_event = event
        return event

    def press(self, event):
        """Button press handler and validator"""
        if not self.ignore(event):
            event = self._clean_event(event)
            self.eventpress = event
            self._prev_event = event
            key = event.key or ''
            key = key.replace('ctrl', 'control')
            # move state is locked in on a button press
            if key == self.state_modifier_keys['move']:
                self.state.add('move')
            self._press(event)
            return True
        return False

    def _press(self, event):
        """Button press handler"""
        pass

    def release(self, event):
        """Button release event handler and validator"""
        if not self.ignore(event) and self.eventpress:
            event = self._clean_event(event)
            self.eventrelease = event
            self._release(event)
            self.eventpress = None
            self.eventrelease = None
            self.state.discard('move')
            return True
        return False

    def _release(self, event):
        """Button release event handler"""
        pass

    def onmove(self, event):
        """Cursor move event handler and validator"""
        if not self.ignore(event) and self.eventpress:
            event = self._clean_event(event)
            self._onmove(event)
            return True
        return False

    def _onmove(self, event):
        """Cursor move event handler"""
        pass

    def on_scroll(self, event):
        """Mouse scroll event handler and validator"""
        if not self.ignore(event):
            self._on_scroll(event)

    def _on_scroll(self, event):
        """Mouse scroll event handler"""
        pass

    def on_key_press(self, event):
        """Key press event handler and validator for all selection widgets"""
        if self.active:
            key = event.key or ''
            key = key.replace('ctrl', 'control')
            if key == self.state_modifier_keys['clear']:
                for artist in self.artists:
                    artist.set_visible(False)
                self.update()
                return
            for (state, modifier) in self.state_modifier_keys.items():
                if modifier in key:
                    self.state.add(state)
            self._on_key_press(event)

    def _on_key_press(self, event):
        """Key press event handler - use for widget-specific key press actions.
        """
        pass

    def on_key_release(self, event):
        """Key release event handler and validator"""
        if self.active:
            key = event.key or ''
            for (state, modifier) in self.state_modifier_keys.items():
                if modifier in key:
                    self.state.discard(state)
            self._on_key_release(event)

    def _on_key_release(self, event):
        """Key release event handler"""
        pass

    def set_visible(self, visible):
        """ Set the visibility of our artists """
        self.visible = visible
        for artist in self.artists:
            artist.set_visible(visible)


class SpanSelector(_SelectorWidget):
    """
    Visually select a min/max range on a single axis and call a function with
    those values.

    To guarantee that the selector remains responsive, keep a reference to it.

    In order to turn off the SpanSelector, set `span_selector.active=False`. To
    turn it back on, set `span_selector.active=True`.

    Parameters
    ----------
    ax :  :class:`matplotlib.axes.Axes` object

    onselect : func(min, max), min/max are floats

    direction : "horizontal" or "vertical"
      The axis along which to draw the span selector

    minspan : float, default is None
     If selection is less than *minspan*, do not call *onselect*

    useblit : bool, default is False
      If True, use the backend-dependent blitting features for faster
      canvas updates.

    rectprops : dict, default is None
      Dictionary of :class:`matplotlib.patches.Patch` properties

    onmove_callback : func(min, max), min/max are floats, default is None
      Called on mouse move while the span is being selected

    span_stays : bool, default is False
      If True, the span stays visible after the mouse is released

    button : int or list of ints
      Determines which mouse buttons activate the span selector
        1 = left mouse button\n
        2 = center mouse button (scroll wheel)\n
        3 = right mouse button\n

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import matplotlib.widgets as mwidgets
    >>> fig, ax = plt.subplots()
    >>> ax.plot([1, 2, 3], [10, 50, 100])
    >>> def onselect(vmin, vmax):
            print(vmin, vmax)
    >>> rectprops = dict(facecolor='blue', alpha=0.5)
    >>> span = mwidgets.SpanSelector(ax, onselect, 'horizontal',
                                     rectprops=rectprops)
    >>> fig.show()

    See also: :ref:`sphx_glr_gallery_widgets_span_selector.py`

    """

    def __init__(self, ax, onselect, direction, minspan=None, useblit=False,
                 rectprops=None, onmove_callback=None, span_stays=False,
                 button=None):

        _SelectorWidget.__init__(self, ax, onselect, useblit=useblit,
                                 button=button)

        if rectprops is None:
            rectprops = dict(facecolor='red', alpha=0.5)

        rectprops['animated'] = self.useblit

        if direction not in ['horizontal', 'vertical']:
            raise ValueError("direction must be 'horizontal' or 'vertical'")
        self.direction = direction

        self.rect = None
        self.pressv = None

        self.rectprops = rectprops
        self.onmove_callback = onmove_callback
        self.minspan = minspan
        self.span_stays = span_stays

        # Needed when dragging out of axes
        self.prev = (0, 0)

        # Reset canvas so that `new_axes` connects events.
        self.canvas = None
        self.new_axes(ax)

    def new_axes(self, ax):
        """Set SpanSelector to operate on a new Axes"""
        self.ax = ax
        if self.canvas is not ax.figure.canvas:
            if self.canvas is not None:
                self.disconnect_events()

            self.canvas = ax.figure.canvas
            self.connect_default_events()

        if self.direction == 'horizontal':
            trans = blended_transform_factory(self.ax.transData,
                                              self.ax.transAxes)
            w, h = 0, 1
        else:
            trans = blended_transform_factory(self.ax.transAxes,
                                              self.ax.transData)
            w, h = 1, 0
        self.rect = Rectangle((0, 0), w, h,
                              transform=trans,
                              visible=False,
                              **self.rectprops)
        if self.span_stays:
            self.stay_rect = Rectangle((0, 0), w, h,
                                       transform=trans,
                                       visible=False,
                                       **self.rectprops)
            self.stay_rect.set_animated(False)
            self.ax.add_patch(self.stay_rect)

        self.ax.add_patch(self.rect)
        self.artists = [self.rect]

    def ignore(self, event):
        """return *True* if *event* should be ignored"""
        return _SelectorWidget.ignore(self, event) or not self.visible

    def _press(self, event):
        """on button press event"""
        self.rect.set_visible(self.visible)
        if self.span_stays:
            self.stay_rect.set_visible(False)
            # really force a draw so that the stay rect is not in
            # the blit background
            if self.useblit:
                self.canvas.draw()
        xdata, ydata = self._get_data(event)
        if self.direction == 'horizontal':
            self.pressv = xdata
        else:
            self.pressv = ydata
        return False

    def _release(self, event):
        """on button release event"""
        if self.pressv is None:
            return
        self.buttonDown = False

        self.rect.set_visible(False)

        if self.span_stays:
            self.stay_rect.set_x(self.rect.get_x())
            self.stay_rect.set_y(self.rect.get_y())
            self.stay_rect.set_width(self.rect.get_width())
            self.stay_rect.set_height(self.rect.get_height())
            self.stay_rect.set_visible(True)

        self.canvas.draw_idle()
        vmin = self.pressv
        xdata, ydata = self._get_data(event)
        if self.direction == 'horizontal':
            vmax = xdata or self.prev[0]
        else:
            vmax = ydata or self.prev[1]

        if vmin > vmax:
            vmin, vmax = vmax, vmin
        span = vmax - vmin
        if self.minspan is not None and span < self.minspan:
            return
        self.onselect(vmin, vmax)
        self.pressv = None
        return False

    def _onmove(self, event):
        """on motion notify event"""
        if self.pressv is None:
            return
        x, y = self._get_data(event)
        if x is None:
            return

        self.prev = x, y
        if self.direction == 'horizontal':
            v = x
        else:
            v = y

        minv, maxv = v, self.pressv
        if minv > maxv:
            minv, maxv = maxv, minv
        if self.direction == 'horizontal':
            self.rect.set_x(minv)
            self.rect.set_width(maxv - minv)
        else:
            self.rect.set_y(minv)
            self.rect.set_height(maxv - minv)

        if self.onmove_callback is not None:
            vmin = self.pressv
            xdata, ydata = self._get_data(event)
            if self.direction == 'horizontal':
                vmax = xdata or self.prev[0]
            else:
                vmax = ydata or self.prev[1]

            if vmin > vmax:
                vmin, vmax = vmax, vmin
            self.onmove_callback(vmin, vmax)

        self.update()
        return False


class ToolHandles(object):
    """Control handles for canvas tools.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes`
        Matplotlib axes where tool handles are displayed.
    x, y : 1D arrays
        Coordinates of control handles.
    marker : str
        Shape of marker used to display handle. See `matplotlib.pyplot.plot`.
    marker_props : dict
        Additional marker properties. See :class:`matplotlib.lines.Line2D`.
    """

    def __init__(self, ax, x, y, marker='o', marker_props=None, useblit=True):
        self.ax = ax

        props = dict(marker=marker, markersize=7, mfc='w', ls='none',
                     alpha=0.5, visible=False, label='_nolegend_')
        props.update(marker_props if marker_props is not None else {})
        self._markers = Line2D(x, y, animated=useblit, **props)
        self.ax.add_line(self._markers)
        self.artist = self._markers

    @property
    def x(self):
        return self._markers.get_xdata()

    @property
    def y(self):
        return self._markers.get_ydata()

    def set_data(self, pts, y=None):
        """Set x and y positions of handles"""
        if y is not None:
            x = pts
            pts = np.array([x, y])
        self._markers.set_data(pts)

    def set_visible(self, val):
        self._markers.set_visible(val)

    def set_animated(self, val):
        self._markers.set_animated(val)

    def closest(self, x, y):
        """Return index and pixel distance to closest index."""
        pts = np.transpose((self.x, self.y))
        # Transform data coordinates to pixel coordinates.
        pts = self.ax.transData.transform(pts)
        diff = pts - ((x, y))
        if diff.ndim == 2:
            dist = np.sqrt(np.sum(diff ** 2, axis=1))
            return np.argmin(dist), np.min(dist)
        else:
            return 0, np.sqrt(np.sum(diff ** 2))


class RectangleSelector(_SelectorWidget):
    """
    Select a rectangular region of an axes.

    For the cursor to remain responsive you must keep a reference to
    it.

    Example usage::

        from matplotlib.widgets import  RectangleSelector
        from pylab import *

        def onselect(eclick, erelease):
          'eclick and erelease are matplotlib events at press and release'
          print(' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata))
          print(' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata))
          print(' used button   : ', eclick.button)

        def toggle_selector(event):
            print(' Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.RS.active:
                print(' RectangleSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.RS.active:
                print(' RectangleSelector activated.')
                toggle_selector.RS.set_active(True)

        x = arange(100)/(99.0)
        y = sin(x)
        fig = figure
        ax = subplot(111)
        ax.plot(x,y)

        toggle_selector.RS = RectangleSelector(ax, onselect, drawtype='line')
        connect('key_press_event', toggle_selector)
        show()
    """

    _shape_klass = Rectangle

    def __init__(self, ax, onselect, drawtype='box',
                 minspanx=None, minspany=None, useblit=False,
                 lineprops=None, rectprops=None, spancoords='data',
                 button=None, maxdist=10, marker_props=None,
                 interactive=False, state_modifier_keys=None):

        """
        Create a selector in *ax*.  When a selection is made, clear
        the span and call onselect with::

          onselect(pos_1, pos_2)

        and clear the drawn box/line. The ``pos_1`` and ``pos_2`` are
        arrays of length 2 containing the x- and y-coordinate.

        If *minspanx* is not *None* then events smaller than *minspanx*
        in x direction are ignored (it's the same for y).

        The rectangle is drawn with *rectprops*; default::

          rectprops = dict(facecolor='red', edgecolor = 'black',
                           alpha=0.2, fill=True)

        The line is drawn with *lineprops*; default::

          lineprops = dict(color='black', linestyle='-',
                           linewidth = 2, alpha=0.5)

        Use *drawtype* if you want the mouse to draw a line,
        a box or nothing between click and actual position by setting

        ``drawtype = 'line'``, ``drawtype='box'`` or ``drawtype = 'none'``.
        Drawing a line would result in a line from vertex A to vertex C in
        a rectangle ABCD.

        *spancoords* is one of 'data' or 'pixels'.  If 'data', *minspanx*
        and *minspanx* will be interpreted in the same coordinates as
        the x and y axis. If 'pixels', they are in pixels.

        *button* is a list of integers indicating which mouse buttons should
        be used for rectangle selection.  You can also specify a single
        integer if only a single button is desired.  Default is *None*,
        which does not limit which button can be used.

        Note, typically:
         1 = left mouse button
         2 = center mouse button (scroll wheel)
         3 = right mouse button

        *interactive* will draw a set of handles and allow you interact
        with the widget after it is drawn.

        *state_modifier_keys* are keyboard modifiers that affect the behavior
        of the widget.

        The defaults are:
        dict(move=' ', clear='escape', square='shift', center='ctrl')

        Keyboard modifiers, which:
        'move': Move the existing shape.
        'clear': Clear the current shape.
        'square': Makes the shape square.
        'center': Make the initial point the center of the shape.
        'square' and 'center' can be combined.
        """
        _SelectorWidget.__init__(self, ax, onselect, useblit=useblit,
                                 button=button,
                                 state_modifier_keys=state_modifier_keys)

        self.to_draw = None
        self.visible = True
        self.interactive = interactive

        if drawtype == 'none':
            drawtype = 'line'                        # draw a line but make it
            self.visible = False                     # invisible

        if drawtype == 'box':
            if rectprops is None:
                rectprops = dict(facecolor='red', edgecolor='black',
                                 alpha=0.2, fill=True)
            rectprops['animated'] = self.useblit
            self.rectprops = rectprops
            self.to_draw = self._shape_klass((0, 0), 0, 1, visible=False,
                                             **self.rectprops)
            self.ax.add_patch(self.to_draw)
        if drawtype == 'line':
            if lineprops is None:
                lineprops = dict(color='black', linestyle='-',
                                 linewidth=2, alpha=0.5)
            lineprops['animated'] = self.useblit
            self.lineprops = lineprops
            self.to_draw = Line2D([0, 0], [0, 0], visible=False,
                                  **self.lineprops)
            self.ax.add_line(self.to_draw)

        self.minspanx = minspanx
        self.minspany = minspany

        if spancoords not in ('data', 'pixels'):
            raise ValueError("'spancoords' must be 'data' or 'pixels'")

        self.spancoords = spancoords
        self.drawtype = drawtype

        self.maxdist = maxdist

        if rectprops is None:
            props = dict(mec='r')
        else:
            props = dict(mec=rectprops.get('edgecolor', 'r'))
        self._corner_order = ['NW', 'NE', 'SE', 'SW']
        xc, yc = self.corners
        self._corner_handles = ToolHandles(self.ax, xc, yc, marker_props=props,
                                           useblit=self.useblit)

        self._edge_order = ['W', 'N', 'E', 'S']
        xe, ye = self.edge_centers
        self._edge_handles = ToolHandles(self.ax, xe, ye, marker='s',
                                         marker_props=props,
                                         useblit=self.useblit)

        xc, yc = self.center
        self._center_handle = ToolHandles(self.ax, [xc], [yc], marker='s',
                                          marker_props=props,
                                          useblit=self.useblit)

        self.active_handle = None

        self.artists = [self.to_draw, self._center_handle.artist,
                        self._corner_handles.artist,
                        self._edge_handles.artist]

        if not self.interactive:
            self.artists = [self.to_draw]

        self._extents_on_press = None

    def _press(self, event):
        """on button press event"""
        # make the drawed box/line visible get the click-coordinates,
        # button, ...
        if self.interactive and self.to_draw.get_visible():
            self._set_active_handle(event)
        else:
            self.active_handle = None

        if self.active_handle is None or not self.interactive:
            # Clear previous rectangle before drawing new rectangle.
            self.update()

        self.set_visible(self.visible)

    def _release(self, event):
        """on button release event"""
        if not self.interactive:
            self.to_draw.set_visible(False)

        # update the eventpress and eventrelease with the resulting extents
        x1, x2, y1, y2 = self.extents
        self.eventpress.xdata = x1
        self.eventpress.ydata = y1
        xy1 = self.ax.transData.transform_point([x1, y1])
        self.eventpress.x, self.eventpress.y = xy1

        self.eventrelease.xdata = x2
        self.eventrelease.ydata = y2
        xy2 = self.ax.transData.transform_point([x2, y2])
        self.eventrelease.x, self.eventrelease.y = xy2

        if self.spancoords == 'data':
            xmin, ymin = self.eventpress.xdata, self.eventpress.ydata
            xmax, ymax = self.eventrelease.xdata, self.eventrelease.ydata
            # calculate dimensions of box or line get values in the right
            # order
        elif self.spancoords == 'pixels':
            xmin, ymin = self.eventpress.x, self.eventpress.y
            xmax, ymax = self.eventrelease.x, self.eventrelease.y
        else:
            raise ValueError('spancoords must be "data" or "pixels"')

        if xmin > xmax:
            xmin, xmax = xmax, xmin
        if ymin > ymax:
            ymin, ymax = ymax, ymin

        spanx = xmax - xmin
        spany = ymax - ymin
        xproblems = self.minspanx is not None and spanx < self.minspanx
        yproblems = self.minspany is not None and spany < self.minspany

        # check if drawn distance (if it exists) is not too small in
        # either x or y-direction
        if self.drawtype != 'none' and (xproblems or yproblems):
            for artist in self.artists:
                artist.set_visible(False)
            self.update()
            return

        # call desired function
        self.onselect(self.eventpress, self.eventrelease)
        self.update()

        return False

    def _onmove(self, event):
        """on motion notify event if box/line is wanted"""
        # resize an existing shape
        if self.active_handle and not self.active_handle == 'C':
            x1, x2, y1, y2 = self._extents_on_press
            if self.active_handle in ['E', 'W'] + self._corner_order:
                x2 = event.xdata
            if self.active_handle in ['N', 'S'] + self._corner_order:
                y2 = event.ydata

        # move existing shape
        elif (('move' in self.state or self.active_handle == 'C')
              and self._extents_on_press is not None):
            x1, x2, y1, y2 = self._extents_on_press
            dx = event.xdata - self.eventpress.xdata
            dy = event.ydata - self.eventpress.ydata
            x1 += dx
            x2 += dx
            y1 += dy
            y2 += dy

        # new shape
        else:
            center = [self.eventpress.xdata, self.eventpress.ydata]
            center_pix = [self.eventpress.x, self.eventpress.y]
            dx = (event.xdata - center[0]) / 2.
            dy = (event.ydata - center[1]) / 2.

            # square shape
            if 'square' in self.state:
                dx_pix = abs(event.x - center_pix[0])
                dy_pix = abs(event.y - center_pix[1])
                if not dx_pix:
                    return
                maxd = max(abs(dx_pix), abs(dy_pix))
                if abs(dx_pix) < maxd:
                    dx *= maxd / (abs(dx_pix) + 1e-6)
                if abs(dy_pix) < maxd:
                    dy *= maxd / (abs(dy_pix) + 1e-6)

            # from center
            if 'center' in self.state:
                dx *= 2
                dy *= 2

            # from corner
            else:
                center[0] += dx
                center[1] += dy

            x1, x2, y1, y2 = (center[0] - dx, center[0] + dx,
                              center[1] - dy, center[1] + dy)

        self.extents = x1, x2, y1, y2

    @property
    def _rect_bbox(self):
        if self.drawtype == 'box':
            x0 = self.to_draw.get_x()
            y0 = self.to_draw.get_y()
            width = self.to_draw.get_width()
            height = self.to_draw.get_height()
            return x0, y0, width, height
        else:
            x, y = self.to_draw.get_data()
            x0, x1 = min(x), max(x)
            y0, y1 = min(y), max(y)
            return x0, y0, x1 - x0, y1 - y0

    @property
    def corners(self):
        """Corners of rectangle from lower left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        xc = x0, x0 + width, x0 + width, x0
        yc = y0, y0, y0 + height, y0 + height
        return xc, yc

    @property
    def edge_centers(self):
        """Midpoint of rectangle edges from left, moving clockwise."""
        x0, y0, width, height = self._rect_bbox
        w = width / 2.
        h = height / 2.
        xe = x0, x0 + w, x0 + width, x0 + w
        ye = y0 + h, y0, y0 + h, y0 + height
        return xe, ye

    @property
    def center(self):
        """Center of rectangle"""
        x0, y0, width, height = self._rect_bbox
        return x0 + width / 2., y0 + height / 2.

    @property
    def extents(self):
        """Return (xmin, xmax, ymin, ymax)."""
        x0, y0, width, height = self._rect_bbox
        xmin, xmax = sorted([x0, x0 + width])
        ymin, ymax = sorted([y0, y0 + height])
        return xmin, xmax, ymin, ymax

    @extents.setter
    def extents(self, extents):
        # Update displayed shape
        self.draw_shape(extents)
        # Update displayed handles
        self._corner_handles.set_data(*self.corners)
        self._edge_handles.set_data(*self.edge_centers)
        self._center_handle.set_data(*self.center)
        self.set_visible(self.visible)
        self.update()

    def draw_shape(self, extents):
        x0, x1, y0, y1 = extents
        xmin, xmax = sorted([x0, x1])
        ymin, ymax = sorted([y0, y1])
        xlim = sorted(self.ax.get_xlim())
        ylim = sorted(self.ax.get_ylim())

        xmin = max(xlim[0], xmin)
        ymin = max(ylim[0], ymin)
        xmax = min(xmax, xlim[1])
        ymax = min(ymax, ylim[1])

        if self.drawtype == 'box':
            self.to_draw.set_x(xmin)
            self.to_draw.set_y(ymin)
            self.to_draw.set_width(xmax - xmin)
            self.to_draw.set_height(ymax - ymin)

        elif self.drawtype == 'line':
            self.to_draw.set_data([xmin, xmax], [ymin, ymax])

    def _set_active_handle(self, event):
        """Set active handle based on the location of the mouse event"""
        # Note: event.xdata/ydata in data coordinates, event.x/y in pixels
        c_idx, c_dist = self._corner_handles.closest(event.x, event.y)
        e_idx, e_dist = self._edge_handles.closest(event.x, event.y)
        m_idx, m_dist = self._center_handle.closest(event.x, event.y)

        if 'move' in self.state:
            self.active_handle = 'C'
            self._extents_on_press = self.extents

        # Set active handle as closest handle, if mouse click is close enough.
        elif m_dist < self.maxdist * 2:
            self.active_handle = 'C'
        elif c_dist > self.maxdist and e_dist > self.maxdist:
            self.active_handle = None
            return
        elif c_dist < e_dist:
            self.active_handle = self._corner_order[c_idx]
        else:
            self.active_handle = self._edge_order[e_idx]

        # Save coordinates of rectangle at the start of handle movement.
        x1, x2, y1, y2 = self.extents
        # Switch variables so that only x2 and/or y2 are updated on move.
        if self.active_handle in ['W', 'SW', 'NW']:
            x1, x2 = x2, event.xdata
        if self.active_handle in ['N', 'NW', 'NE']:
            y1, y2 = y2, event.ydata
        self._extents_on_press = x1, x2, y1, y2

    @property
    def geometry(self):
        """
        Returns numpy.ndarray of shape (2,5) containing
        x (``RectangleSelector.geometry[1,:]``) and
        y (``RectangleSelector.geometry[0,:]``)
        coordinates of the four corners of the rectangle starting
        and ending in the top left corner.
        """
        if hasattr(self.to_draw, 'get_verts'):
            xfm = self.ax.transData.inverted()
            y, x = xfm.transform(self.to_draw.get_verts()).T
            return np.array([x, y])
        else:
            return np.array(self.to_draw.get_data())


class EllipseSelector(RectangleSelector):
    """
    Select an elliptical region of an axes.

    For the cursor to remain responsive you must keep a reference to
    it.

    Example usage::

        from matplotlib.widgets import  EllipseSelector
        from pylab import *

        def onselect(eclick, erelease):
          'eclick and erelease are matplotlib events at press and release'
          print(' startposition : (%f, %f)' % (eclick.xdata, eclick.ydata))
          print(' endposition   : (%f, %f)' % (erelease.xdata, erelease.ydata))
          print(' used button   : ', eclick.button)

        def toggle_selector(event):
            print(' Key pressed.')
            if event.key in ['Q', 'q'] and toggle_selector.ES.active:
                print(' EllipseSelector deactivated.')
                toggle_selector.RS.set_active(False)
            if event.key in ['A', 'a'] and not toggle_selector.ES.active:
                print(' EllipseSelector activated.')
                toggle_selector.ES.set_active(True)

        x = arange(100)/(99.0)
        y = sin(x)
        fig = figure
        ax = subplot(111)
        ax.plot(x,y)

        toggle_selector.ES = EllipseSelector(ax, onselect, drawtype='line')
        connect('key_press_event', toggle_selector)
        show()
    """
    _shape_klass = Ellipse

    def draw_shape(self, extents):
        x1, x2, y1, y2 = extents
        xmin, xmax = sorted([x1, x2])
        ymin, ymax = sorted([y1, y2])
        center = [x1 + (x2 - x1) / 2., y1 + (y2 - y1) / 2.]
        a = (xmax - xmin) / 2.
        b = (ymax - ymin) / 2.

        if self.drawtype == 'box':
            self.to_draw.center = center
            self.to_draw.width = 2 * a
            self.to_draw.height = 2 * b
        else:
            rad = np.deg2rad(np.arange(31) * 12)
            x = a * np.cos(rad) + center[0]
            y = b * np.sin(rad) + center[1]
            self.to_draw.set_data(x, y)

    @property
    def _rect_bbox(self):
        if self.drawtype == 'box':
            x, y = self.to_draw.center
            width = self.to_draw.width
            height = self.to_draw.height
            return x - width / 2., y - height / 2., width, height
        else:
            x, y = self.to_draw.get_data()
            x0, x1 = min(x), max(x)
            y0, y1 = min(y), max(y)
            return x0, y0, x1 - x0, y1 - y0


class LassoSelector(_SelectorWidget):
    """
    Selection curve of an arbitrary shape.

    For the selector to remain responsive you must keep a reference to it.

    The selected path can be used in conjunction with `~.Path.contains_point`
    to select data points from an image.

    In contrast to `Lasso`, `LassoSelector` is written with an interface
    similar to `RectangleSelector` and `SpanSelector`, and will continue to
    interact with the axes until disconnected.

    Example usage::

        ax = subplot(111)
        ax.plot(x,y)

        def onselect(verts):
            print(verts)
        lasso = LassoSelector(ax, onselect)

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        The parent axes for the widget.
    onselect : function
        Whenever the lasso is released, the *onselect* function is called and
        passed the vertices of the selected path.
    button : List[Int], optional
        A list of integers indicating which mouse buttons should be used for
        rectangle selection. You can also specify a single integer if only a
        single button is desired.  Default is ``None``, which does not limit
        which button can be used.

        Note, typically:

        - 1 = left mouse button
        - 2 = center mouse button (scroll wheel)
        - 3 = right mouse button

    """

    def __init__(self, ax, onselect=None, useblit=True, lineprops=None,
                 button=None):
        _SelectorWidget.__init__(self, ax, onselect, useblit=useblit,
                                 button=button)

        self.verts = None

        if lineprops is None:
            lineprops = dict()
        if useblit:
            lineprops['animated'] = True
        self.line = Line2D([], [], **lineprops)
        self.line.set_visible(False)
        self.ax.add_line(self.line)
        self.artists = [self.line]

    def onpress(self, event):
        self.press(event)

    def _press(self, event):
        self.verts = [self._get_data(event)]
        self.line.set_visible(True)

    def onrelease(self, event):
        self.release(event)

    def _release(self, event):
        if self.verts is not None:
            self.verts.append(self._get_data(event))
            self.onselect(self.verts)
        self.line.set_data([[], []])
        self.line.set_visible(False)
        self.verts = None

    def _onmove(self, event):
        if self.verts is None:
            return
        self.verts.append(self._get_data(event))

        self.line.set_data(list(zip(*self.verts)))

        self.update()


class PolygonSelector(_SelectorWidget):
    """Select a polygon region of an axes.

    Place vertices with each mouse click, and make the selection by completing
    the polygon (clicking on the first vertex). Hold the *ctrl* key and click
    and drag a vertex to reposition it (the *ctrl* key is not necessary if the
    polygon has already been completed). Hold the *shift* key and click and
    drag anywhere in the axes to move all vertices. Press the *esc* key to
    start a new polygon.

    For the selector to remain responsive you must keep a reference to
    it.

    Parameters
    ----------
    ax : :class:`~matplotlib.axes.Axes`
        The parent axes for the widget.
    onselect : function
        When a polygon is completed or modified after completion,
        the `onselect` function is called and passed a list of the vertices as
        ``(xdata, ydata)`` tuples.
    useblit : bool, optional
    lineprops : dict, optional
        The line for the sides of the polygon is drawn with the properties
        given by `lineprops`. The default is ``dict(color='k', linestyle='-',
        linewidth=2, alpha=0.5)``.
    markerprops : dict, optional
        The markers for the vertices of the polygon are drawn with the
        properties given by `markerprops`. The default is ``dict(marker='o',
        markersize=7, mec='k', mfc='k', alpha=0.5)``.
    vertex_select_radius : float, optional
        A vertex is selected (to complete the polygon or to move a vertex)
        if the mouse click is within `vertex_select_radius` pixels of the
        vertex. The default radius is 15 pixels.

    See Also
    --------
    :ref:`sphx_glr_gallery_widgets_polygon_selector_demo.py`
    """

    def __init__(self, ax, onselect, useblit=False,
                 lineprops=None, markerprops=None, vertex_select_radius=15):
        # The state modifiers 'move', 'square', and 'center' are expected by
        # _SelectorWidget but are not supported by PolygonSelector
        # Note: could not use the existing 'move' state modifier in-place of
        # 'move_all' because _SelectorWidget automatically discards 'move'
        # from the state on button release.
        state_modifier_keys = dict(clear='escape', move_vertex='control',
                                   move_all='shift', move='not-applicable',
                                   square='not-applicable',
                                   center='not-applicable')
        _SelectorWidget.__init__(self, ax, onselect, useblit=useblit,
                                 state_modifier_keys=state_modifier_keys)

        self._xs, self._ys = [0], [0]
        self._polygon_completed = False

        if lineprops is None:
            lineprops = dict(color='k', linestyle='-', linewidth=2, alpha=0.5)
        lineprops['animated'] = self.useblit
        self.line = Line2D(self._xs, self._ys, **lineprops)
        self.ax.add_line(self.line)

        if markerprops is None:
            markerprops = dict(mec='k', mfc=lineprops.get('color', 'k'))
        self._polygon_handles = ToolHandles(self.ax, self._xs, self._ys,
                                            useblit=self.useblit,
                                            marker_props=markerprops)

        self._active_handle_idx = -1
        self.vertex_select_radius = vertex_select_radius

        self.artists = [self.line, self._polygon_handles.artist]
        self.set_visible(True)

    def _press(self, event):
        """Button press event handler"""
        # Check for selection of a tool handle.
        if ((self._polygon_completed or 'move_vertex' in self.state)
                and len(self._xs) > 0):
            h_idx, h_dist = self._polygon_handles.closest(event.x, event.y)
            if h_dist < self.vertex_select_radius:
                self._active_handle_idx = h_idx
        # Save the vertex positions at the time of the press event (needed to
        # support the 'move_all' state modifier).
        self._xs_at_press, self._ys_at_press = self._xs[:], self._ys[:]

    def _release(self, event):
        """Button release event handler"""
        # Release active tool handle.
        if self._active_handle_idx >= 0:
            self._active_handle_idx = -1

        # Complete the polygon.
        elif (len(self._xs) > 3
              and self._xs[-1] == self._xs[0]
              and self._ys[-1] == self._ys[0]):
            self._polygon_completed = True

        # Place new vertex.
        elif (not self._polygon_completed
              and 'move_all' not in self.state
              and 'move_vertex' not in self.state):
            self._xs.insert(-1, event.xdata)
            self._ys.insert(-1, event.ydata)

        if self._polygon_completed:
            self.onselect(self.verts)

    def onmove(self, event):
        """Cursor move event handler and validator"""
        # Method overrides _SelectorWidget.onmove because the polygon selector
        # needs to process the move callback even if there is no button press.
        # _SelectorWidget.onmove include logic to ignore move event if
        # eventpress is None.
        if not self.ignore(event):
            event = self._clean_event(event)
            self._onmove(event)
            return True
        return False

    def _onmove(self, event):
        """Cursor move event handler"""
        # Move the active vertex (ToolHandle).
        if self._active_handle_idx >= 0:
            idx = self._active_handle_idx
            self._xs[idx], self._ys[idx] = event.xdata, event.ydata
            # Also update the end of the polygon line if the first vertex is
            # the active handle and the polygon is completed.
            if idx == 0 and self._polygon_completed:
                self._xs[-1], self._ys[-1] = event.xdata, event.ydata

        # Move all vertices.
        elif 'move_all' in self.state and self.eventpress:
            dx = event.xdata - self.eventpress.xdata
            dy = event.ydata - self.eventpress.ydata
            for k in range(len(self._xs)):
                self._xs[k] = self._xs_at_press[k] + dx
                self._ys[k] = self._ys_at_press[k] + dy

        # Do nothing if completed or waiting for a move.
        elif (self._polygon_completed
              or 'move_vertex' in self.state or 'move_all' in self.state):
            return

        # Position pending vertex.
        else:
            # Calculate distance to the start vertex.
            x0, y0 = self.line.get_transform().transform((self._xs[0],
                                                          self._ys[0]))
            v0_dist = np.sqrt((x0 - event.x) ** 2 + (y0 - event.y) ** 2)
            # Lock on to the start vertex if near it and ready to complete.
            if len(self._xs) > 3 and v0_dist < self.vertex_select_radius:
                self._xs[-1], self._ys[-1] = self._xs[0], self._ys[0]
            else:
                self._xs[-1], self._ys[-1] = event.xdata, event.ydata

        self._draw_polygon()

    def _on_key_press(self, event):
        """Key press event handler"""
        # Remove the pending vertex if entering the 'move_vertex' or
        # 'move_all' mode
        if (not self._polygon_completed
                and ('move_vertex' in self.state or 'move_all' in self.state)):
            self._xs, self._ys = self._xs[:-1], self._ys[:-1]
            self._draw_polygon()

    def _on_key_release(self, event):
        """Key release event handler"""
        # Add back the pending vertex if leaving the 'move_vertex' or
        # 'move_all' mode (by checking the released key)
        if (not self._polygon_completed
                and
                (event.key == self.state_modifier_keys.get('move_vertex')
                 or event.key == self.state_modifier_keys.get('move_all'))):
            self._xs.append(event.xdata)
            self._ys.append(event.ydata)
            self._draw_polygon()
        # Reset the polygon if the released key is the 'clear' key.
        elif event.key == self.state_modifier_keys.get('clear'):
            event = self._clean_event(event)
            self._xs, self._ys = [event.xdata], [event.ydata]
            self._polygon_completed = False
            self.set_visible(True)

    def _draw_polygon(self):
        """Redraw the polygon based on the new vertex positions."""
        self.line.set_data(self._xs, self._ys)
        # Only show one tool handle at the start and end vertex of the polygon
        # if the polygon is completed or the user is locked on to the start
        # vertex.
        if (self._polygon_completed
                or (len(self._xs) > 3
                    and self._xs[-1] == self._xs[0]
                    and self._ys[-1] == self._ys[0])):
            self._polygon_handles.set_data(self._xs[:-1], self._ys[:-1])
        else:
            self._polygon_handles.set_data(self._xs, self._ys)
        self.update()

    @property
    def verts(self):
        """Get the polygon vertices.

        Returns
        -------
        list
            A list of the vertices of the polygon as ``(xdata, ydata)`` tuples.
        """
        return list(zip(self._xs[:-1], self._ys[:-1]))


class Lasso(AxesWidget):
    """Selection curve of an arbitrary shape.

    The selected path can be used in conjunction with
    :func:`~matplotlib.path.Path.contains_point` to select data points
    from an image.

    Unlike :class:`LassoSelector`, this must be initialized with a starting
    point `xy`, and the `Lasso` events are destroyed upon release.

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        The parent axes for the widget.
    xy : array
        Coordinates of the start of the lasso.
    callback : callable
        Whenever the lasso is released, the `callback` function is called and
        passed the vertices of the selected path.
    """

    def __init__(self, ax, xy, callback=None, useblit=True):
        AxesWidget.__init__(self, ax)

        self.useblit = useblit and self.canvas.supports_blit
        if self.useblit:
            self.background = self.canvas.copy_from_bbox(self.ax.bbox)

        x, y = xy
        self.verts = [(x, y)]
        self.line = Line2D([x], [y], linestyle='-', color='black', lw=2)
        self.ax.add_line(self.line)
        self.callback = callback
        self.connect_event('button_release_event', self.onrelease)
        self.connect_event('motion_notify_event', self.onmove)

    def onrelease(self, event):
        if self.ignore(event):
            return
        if self.verts is not None:
            self.verts.append((event.xdata, event.ydata))
            if len(self.verts) > 2:
                self.callback(self.verts)
            self.ax.lines.remove(self.line)
        self.verts = None
        self.disconnect_events()

    def onmove(self, event):
        if self.ignore(event):
            return
        if self.verts is None:
            return
        if event.inaxes != self.ax:
            return
        if event.button != 1:
            return
        self.verts.append((event.xdata, event.ydata))

        self.line.set_data(list(zip(*self.verts)))

        if self.useblit:
            self.canvas.restore_region(self.background)
            self.ax.draw_artist(self.line)
            self.canvas.blit(self.ax.bbox)
        else:
            self.canvas.draw_idle()
