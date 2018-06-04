"""
Abstraction of CLI Input.
"""
from __future__ import unicode_literals

from .utils import DummyContext, is_windows
from abc import ABCMeta, abstractmethod
from six import with_metaclass

import io
import os
import sys

if is_windows():
    from .terminal.win32_input import raw_mode, cooked_mode
else:
    from .terminal.vt100_input import raw_mode, cooked_mode

__all__ = (
    'Input',
    'StdinInput',
    'PipeInput',
)


class Input(with_metaclass(ABCMeta, object)):
    """
    Abstraction for any input.

    An instance of this class can be given to the constructor of a
    :class:`~prompt_toolkit.interface.CommandLineInterface` and will also be
    passed to the :class:`~prompt_toolkit.eventloop.base.EventLoop`.
    """
    @abstractmethod
    def fileno(self):
        """
        Fileno for putting this in an event loop.
        """

    @abstractmethod
    def read(self):
        """
        Return text from the input.
        """

    @abstractmethod
    def raw_mode(self):
        """
        Context manager that turns the input into raw mode.
        """

    @abstractmethod
    def cooked_mode(self):
        """
        Context manager that turns the input into cooked mode.
        """


class StdinInput(Input):
    """
    Simple wrapper around stdin.
    """
    def __init__(self, stdin=None):
        self.stdin = stdin or sys.stdin

        # The input object should be a TTY.
        assert self.stdin.isatty()

        # Test whether the given input object has a file descriptor.
        # (Idle reports stdin to be a TTY, but fileno() is not implemented.)
        try:
            # This should not raise, but can return 0.
            self.stdin.fileno()
        except io.UnsupportedOperation:
            if 'idlelib.run' in sys.modules:
                raise io.UnsupportedOperation(
                    'Stdin is not a terminal. Running from Idle is not supported.')
            else:
                raise io.UnsupportedOperation('Stdin is not a terminal.')

    def __repr__(self):
        return 'StdinInput(stdin=%r)' % (self.stdin,)

    def raw_mode(self):
        return raw_mode(self.stdin.fileno())

    def cooked_mode(self):
        return cooked_mode(self.stdin.fileno())

    def fileno(self):
        return self.stdin.fileno()

    def read(self):
        return self.stdin.read()


class PipeInput(Input):
    """
    Input that is send through a pipe.
    This is useful if we want to send the input programatically into the
    interface, but still use the eventloop.

    Usage::

        input = PipeInput()
        input.send('inputdata')
    """
    def __init__(self):
        self._r, self._w = os.pipe()

    def fileno(self):
        return self._r

    def read(self):
        return os.read(self._r)

    def send_text(self, data):
        " Send text to the input. "
        os.write(self._w, data.encode('utf-8'))

    # Deprecated alias for `send_text`.
    send = send_text

    def raw_mode(self):
        return DummyContext()

    def cooked_mode(self):
        return DummyContext()

    def close(self):
        " Close pipe fds. "
        os.close(self._r)
        os.close(self._w)
        self._r = None
        self._w = None
