"""
This provides several classes used for blocking interaction with figure
windows:

:class:`BlockingInput`
    creates a callable object to retrieve events in a blocking way for
    interactive sessions

:class:`BlockingKeyMouseInput`
    creates a callable object to retrieve key or mouse clicks in a blocking
    way for interactive sessions.
    Note: Subclass of BlockingInput. Used by waitforbuttonpress

:class:`BlockingMouseInput`
    creates a callable object to retrieve mouse clicks in a blocking way for
    interactive sessions.
    Note: Subclass of BlockingInput.  Used by ginput

:class:`BlockingContourLabeler`
    creates a callable object to retrieve mouse clicks in a blocking way that
    will then be used to place labels on a ContourSet
    Note: Subclass of BlockingMouseInput.  Used by clabel
"""

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import six
import matplotlib.lines as mlines

import logging

_log = logging.getLogger(__name__)


class BlockingInput(object):
    """
    Class that creates a callable object to retrieve events in a
    blocking way.
    """
    def __init__(self, fig, eventslist=()):
        self.fig = fig
        self.eventslist = eventslist

    def on_event(self, event):
        """
        Event handler that will be passed to the current figure to
        retrieve events.
        """
        # Add a new event to list - using a separate function is
        # overkill for the base class, but this is consistent with
        # subclasses
        self.add_event(event)
        _log.info("Event %i", len(self.events))

        # This will extract info from events
        self.post_event()

        # Check if we have enough events already
        if len(self.events) >= self.n and self.n > 0:
            self.fig.canvas.stop_event_loop()

    def post_event(self):
        """For baseclass, do nothing but collect events"""
        pass

    def cleanup(self):
        """Disconnect all callbacks"""
        for cb in self.callbacks:
            self.fig.canvas.mpl_disconnect(cb)

        self.callbacks = []

    def add_event(self, event):
        """For base class, this just appends an event to events."""
        self.events.append(event)

    def pop_event(self, index=-1):
        """
        This removes an event from the event list.  Defaults to
        removing last event, but an index can be supplied.  Note that
        this does not check that there are events, much like the
        normal pop method.  If not events exist, this will throw an
        exception.
        """
        self.events.pop(index)

    def pop(self, index=-1):
        self.pop_event(index)
    pop.__doc__ = pop_event.__doc__

    def __call__(self, n=1, timeout=30):
        """
        Blocking call to retrieve n events
        """

        if not isinstance(n, int):
            raise ValueError("Requires an integer argument")
        self.n = n

        self.events = []
        self.callbacks = []

        if hasattr(self.fig.canvas, "manager"):
            # Ensure that the figure is shown, if we are managing it.
            self.fig.show()

        # connect the events to the on_event function call
        for n in self.eventslist:
            self.callbacks.append(
                self.fig.canvas.mpl_connect(n, self.on_event))

        try:
            # Start event loop
            self.fig.canvas.start_event_loop(timeout=timeout)
        finally:  # Run even on exception like ctrl-c
            # Disconnect the callbacks
            self.cleanup()

        # Return the events in this case
        return self.events


class BlockingMouseInput(BlockingInput):
    """
    Class that creates a callable object to retrieve mouse clicks in a
    blocking way.

    This class will also retrieve keyboard clicks and treat them like
    appropriate mouse clicks (delete and backspace are like mouse button 3,
    enter is like mouse button 2 and all others are like mouse button 1).
    """

    button_add = 1
    button_pop = 3
    button_stop = 2

    def __init__(self, fig, mouse_add=1, mouse_pop=3, mouse_stop=2):
        BlockingInput.__init__(self, fig=fig,
                               eventslist=('button_press_event',
                                           'key_press_event'))
        self.button_add = mouse_add
        self.button_pop = mouse_pop
        self.button_stop = mouse_stop

    def post_event(self):
        """
        This will be called to process events
        """
        if len(self.events) == 0:
            _log.warning("No events yet")
        elif self.events[-1].name == 'key_press_event':
            self.key_event()
        else:
            self.mouse_event()

    def mouse_event(self):
        '''Process a mouse click event'''

        event = self.events[-1]
        button = event.button

        if button == self.button_pop:
            self.mouse_event_pop(event)
        elif button == self.button_stop:
            self.mouse_event_stop(event)
        else:
            self.mouse_event_add(event)

    def key_event(self):
        '''
        Process a key click event.  This maps certain keys to appropriate
        mouse click events.
        '''

        event = self.events[-1]
        if event.key is None:
            # at least in mac os X gtk backend some key returns None.
            return

        key = event.key.lower()

        if key in ['backspace', 'delete']:
            self.mouse_event_pop(event)
        elif key in ['escape', 'enter']:
            # on windows XP and wxAgg, the enter key doesn't seem to register
            self.mouse_event_stop(event)
        else:
            self.mouse_event_add(event)

    def mouse_event_add(self, event):
        """
        Will be called for any event involving a button other than
        button 2 or 3.  This will add a click if it is inside axes.
        """
        if event.inaxes:
            self.add_click(event)
        else:  # If not a valid click, remove from event list
            BlockingInput.pop(self, -1)

    def mouse_event_stop(self, event):
        """
        Will be called for any event involving button 2.
        Button 2 ends blocking input.
        """

        # Remove last event just for cleanliness
        BlockingInput.pop(self, -1)

        # This will exit even if not in infinite mode.  This is
        # consistent with MATLAB and sometimes quite useful, but will
        # require the user to test how many points were actually
        # returned before using data.
        self.fig.canvas.stop_event_loop()

    def mouse_event_pop(self, event):
        """
        Will be called for any event involving button 3.
        Button 3 removes the last click.
        """
        # Remove this last event
        BlockingInput.pop(self, -1)

        # Now remove any existing clicks if possible
        if len(self.events) > 0:
            self.pop(event, -1)

    def add_click(self, event):
        """
        This add the coordinates of an event to the list of clicks
        """
        self.clicks.append((event.xdata, event.ydata))
        _log.info("input %i: %f,%f" %
                       (len(self.clicks), event.xdata, event.ydata))

        # If desired plot up click
        if self.show_clicks:
            line = mlines.Line2D([event.xdata], [event.ydata],
                                 marker='+', color='r')
            event.inaxes.add_line(line)
            self.marks.append(line)
            self.fig.canvas.draw()

    def pop_click(self, event, index=-1):
        """
        This removes a click from the list of clicks.  Defaults to
        removing the last click.
        """
        self.clicks.pop(index)

        if self.show_clicks:

            mark = self.marks.pop(index)
            mark.remove()

            self.fig.canvas.draw()
            # NOTE: I do NOT understand why the above 3 lines does not work
            # for the keyboard backspace event on windows XP wxAgg.
            # maybe event.inaxes here is a COPY of the actual axes?

    def pop(self, event, index=-1):
        """
        This removes a click and the associated event from the object.
        Defaults to removing the last click, but any index can be
        supplied.
        """
        self.pop_click(event, index)
        BlockingInput.pop(self, index)

    def cleanup(self, event=None):
        # clean the figure
        if self.show_clicks:

            for mark in self.marks:
                mark.remove()
            self.marks = []

            self.fig.canvas.draw()

        # Call base class to remove callbacks
        BlockingInput.cleanup(self)

    def __call__(self, n=1, timeout=30, show_clicks=True):
        """
        Blocking call to retrieve n coordinate pairs through mouse
        clicks.
        """
        self.show_clicks = show_clicks
        self.clicks = []
        self.marks = []
        BlockingInput.__call__(self, n=n, timeout=timeout)

        return self.clicks


class BlockingContourLabeler(BlockingMouseInput):
    """
    Class that creates a callable object that uses mouse clicks or key
    clicks on a figure window to place contour labels.
    """
    def __init__(self, cs):
        self.cs = cs
        BlockingMouseInput.__init__(self, fig=cs.ax.figure)

    def add_click(self, event):
        self.button1(event)

    def pop_click(self, event, index=-1):
        self.button3(event)

    def button1(self, event):
        """
        This will be called if an event involving a button other than
        2 or 3 occcurs.  This will add a label to a contour.
        """

        # Shorthand
        if event.inaxes == self.cs.ax:
            self.cs.add_label_near(event.x, event.y, self.inline,
                                   inline_spacing=self.inline_spacing,
                                   transform=False)
            self.fig.canvas.draw()
        else:  # Remove event if not valid
            BlockingInput.pop(self)

    def button3(self, event):
        """
        This will be called if button 3 is clicked.  This will remove
        a label if not in inline mode.  Unfortunately, if one is doing
        inline labels, then there is currently no way to fix the
        broken contour - once humpty-dumpty is broken, he can't be put
        back together.  In inline mode, this does nothing.
        """

        if self.inline:
            pass
        else:
            self.cs.pop_label()
            self.cs.ax.figure.canvas.draw()

    def __call__(self, inline, inline_spacing=5, n=-1, timeout=-1):
        self.inline = inline
        self.inline_spacing = inline_spacing

        BlockingMouseInput.__call__(self, n=n, timeout=timeout,
                                    show_clicks=False)


class BlockingKeyMouseInput(BlockingInput):
    """
    Class that creates a callable object to retrieve a single mouse or
    keyboard click
    """
    def __init__(self, fig):
        BlockingInput.__init__(self, fig=fig, eventslist=(
            'button_press_event', 'key_press_event'))

    def post_event(self):
        """
        Determines if it is a key event
        """
        if len(self.events) == 0:
            _log.warning("No events yet")
        else:
            self.keyormouse = self.events[-1].name == 'key_press_event'

    def __call__(self, timeout=30):
        """
        Blocking call to retrieve a single mouse or key click
        Returns True if key click, False if mouse, or None if timeout
        """
        self.keyormouse = None
        BlockingInput.__call__(self, n=1, timeout=timeout)

        return self.keyormouse
