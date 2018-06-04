"""
Mouse events.


How it works
------------

The renderer has a 2 dimensional grid of mouse event handlers.
(`prompt_toolkit.layout.MouseHandlers`.) When the layout is rendered, the
`Window` class will make sure that this grid will also be filled with
callbacks. For vt100 terminals, mouse events are received through stdin, just
like any other key press. There is a handler among the key bindings that
catches these events and forwards them to such a mouse event handler. It passes
through the `Window` class where the coordinates are translated from absolute
coordinates to coordinates relative to the user control, and there
`UIControl.mouse_handler` is called.
"""
from __future__ import unicode_literals

__all__ = (
    'MouseEventType',
    'MouseEvent'
)


class MouseEventType:
    MOUSE_UP = 'MOUSE_UP'
    MOUSE_DOWN = 'MOUSE_DOWN'
    SCROLL_UP = 'SCROLL_UP'
    SCROLL_DOWN = 'SCROLL_DOWN'


MouseEventTypes = MouseEventType  # Deprecated: plural for backwards compatibility.


class MouseEvent(object):
    """
    Mouse event, sent to `UIControl.mouse_handler`.

    :param position: `Point` instance.
    :param event_type: `MouseEventType`.
    """
    def __init__(self, position, event_type):
        self.position = position
        self.event_type = event_type

    def __repr__(self):
        return 'MouseEvent(%r, %r)' % (self.position, self.event_type)
