"""
Eventloop for integration with Python3 asyncio.

Note that we can't use "yield from", because the package should be installable
under Python 2.6 as well, and it should contain syntactically valid Python 2.6
code.
"""
from __future__ import unicode_literals

__all__ = (
    'AsyncioTimeout',
)


class AsyncioTimeout(object):
    """
    Call the `timeout` function when the timeout expires.
    Every call of the `reset` method, resets the timeout and starts a new
    timer.
    """
    def __init__(self, timeout, callback, loop):
        self.timeout = timeout
        self.callback = callback
        self.loop = loop

        self.counter = 0
        self.running = True

    def reset(self):
        """
        Reset the timeout. Starts a new timer.
        """
        self.counter += 1
        local_counter = self.counter

        def timer_timeout():
            if self.counter == local_counter and self.running:
                self.callback()

        self.loop.call_later(self.timeout, timer_timeout)

    def stop(self):
        """
        Ignore timeout. Don't call the callback anymore.
        """
        self.running = False
