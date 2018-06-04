from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from six import with_metaclass

__all__ = (
    'EventLoopCallbacks',
)


class EventLoopCallbacks(with_metaclass(ABCMeta, object)):
    """
    This is the glue between the :class:`~prompt_toolkit.eventloop.base.EventLoop`
    and :class:`~prompt_toolkit.interface.CommandLineInterface`.

    :meth:`~prompt_toolkit.eventloop.base.EventLoop.run` takes an
    :class:`.EventLoopCallbacks` instance and operates on that one, driving the
    interface.
    """
    @abstractmethod
    def terminal_size_changed(self):
        pass

    @abstractmethod
    def input_timeout(self):
        pass

    @abstractmethod
    def feed_key(self, key):
        pass
