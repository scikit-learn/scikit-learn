from __future__ import unicode_literals

from itertools import product
from collections import defaultdict

__all__ = (
    'MouseHandlers',
)


class MouseHandlers(object):
    """
    Two dimentional raster of callbacks for mouse events.
    """
    def __init__(self):
        def dummy_callback(cli, mouse_event):
            """
            :param mouse_event: `MouseEvent` instance.
            """

        # Map (x,y) tuples to handlers.
        self.mouse_handlers = defaultdict(lambda: dummy_callback)

    def set_mouse_handler_for_range(self, x_min, x_max, y_min, y_max, handler=None):
        """
        Set mouse handler for a region.
        """
        for x, y in product(range(x_min, x_max), range(y_min, y_max)):
            self.mouse_handlers[x,y] = handler
