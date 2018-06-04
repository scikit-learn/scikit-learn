"""
Posix asyncio event loop.
"""
from __future__ import unicode_literals

from ..terminal.vt100_input import InputStream
from .asyncio_base import AsyncioTimeout
from .base import EventLoop, INPUT_TIMEOUT
from .callbacks import EventLoopCallbacks
from .posix_utils import PosixStdinReader

import asyncio
import signal

__all__ = (
    'PosixAsyncioEventLoop',
)


class PosixAsyncioEventLoop(EventLoop):
    def __init__(self, loop=None):
        self.loop = loop or asyncio.get_event_loop()
        self.closed = False

        self._stopped_f = asyncio.Future(loop=self.loop)

    @asyncio.coroutine
    def run_as_coroutine(self, stdin, callbacks):
        """
        The input 'event loop'.
        """
        assert isinstance(callbacks, EventLoopCallbacks)

        # Create reader class.
        stdin_reader = PosixStdinReader(stdin.fileno())

        if self.closed:
            raise Exception('Event loop already closed.')

        inputstream = InputStream(callbacks.feed_key)

        try:
            # Create a new Future every time.
            self._stopped_f = asyncio.Future(loop=self.loop)

            # Handle input timouts
            def timeout_handler():
                """
                When no input has been received for INPUT_TIMEOUT seconds,
                flush the input stream and fire the timeout event.
                """
                inputstream.flush()

                callbacks.input_timeout()

            timeout = AsyncioTimeout(INPUT_TIMEOUT, timeout_handler, self.loop)

            # Catch sigwinch
            def received_winch():
                self.call_from_executor(callbacks.terminal_size_changed)

            self.loop.add_signal_handler(signal.SIGWINCH, received_winch)

            # Read input data.
            def stdin_ready():
                data = stdin_reader.read()
                inputstream.feed(data)
                timeout.reset()

                # Quit when the input stream was closed.
                if stdin_reader.closed:
                    self.stop()

            self.loop.add_reader(stdin.fileno(), stdin_ready)

            # Block this coroutine until stop() has been called.
            for f in self._stopped_f:
                yield f

        finally:
            # Clean up.
            self.loop.remove_reader(stdin.fileno())
            self.loop.remove_signal_handler(signal.SIGWINCH)

            # Don't trigger any timeout events anymore.
            timeout.stop()

    def stop(self):
        # Trigger the 'Stop' future.
        self._stopped_f.set_result(True)

    def close(self):
        # Note: we should not close the asyncio loop itself, because that one
        # was not created here.
        self.closed = True

    def run_in_executor(self, callback):
        self.loop.run_in_executor(None, callback)

    def call_from_executor(self, callback, _max_postpone_until=None):
        """
        Call this function in the main event loop.
        Similar to Twisted's ``callFromThread``.
        """
        self.loop.call_soon_threadsafe(callback)

    def add_reader(self, fd, callback):
        " Start watching the file descriptor for read availability. "
        self.loop.add_reader(fd, callback)

    def remove_reader(self, fd):
        " Stop watching the file descriptor for read availability. "
        self.loop.remove_reader(fd)
