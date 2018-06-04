from __future__ import unicode_literals
from abc import ABCMeta, abstractmethod
from six import with_metaclass

__all__ = (
    'EventLoop',
    'INPUT_TIMEOUT',
)


#: When to trigger the `onInputTimeout` event.
INPUT_TIMEOUT = .5


class EventLoop(with_metaclass(ABCMeta, object)):
    """
    Eventloop interface.
    """
    def run(self, stdin, callbacks):
        """
        Run the eventloop until stop() is called. Report all
        input/timeout/terminal-resize events to the callbacks.

        :param stdin: :class:`~prompt_toolkit.input.Input` instance.
        :param callbacks: :class:`~prompt_toolkit.eventloop.callbacks.EventLoopCallbacks` instance.
        """
        raise NotImplementedError("This eventloop doesn't implement synchronous 'run()'.")

    def run_as_coroutine(self, stdin, callbacks):
        """
        Similar to `run`, but this is a coroutine. (For asyncio integration.)
        """
        raise NotImplementedError("This eventloop doesn't implement 'run_as_coroutine()'.")

    @abstractmethod
    def stop(self):
        """
        Stop the `run` call. (Normally called by
        :class:`~prompt_toolkit.interface.CommandLineInterface`, when a result
        is available, or Abort/Quit has been called.)
        """

    @abstractmethod
    def close(self):
        """
        Clean up of resources. Eventloop cannot be reused a second time after
        this call.
        """

    @abstractmethod
    def add_reader(self, fd, callback):
        """
        Start watching the file descriptor for read availability and then call
        the callback.
        """

    @abstractmethod
    def remove_reader(self, fd):
        """
        Stop watching the file descriptor for read availability.
        """

    @abstractmethod
    def run_in_executor(self, callback):
        """
        Run a long running function in a background thread. (This is
        recommended for code that could block the event loop.)
        Similar to Twisted's ``deferToThread``.
        """

    @abstractmethod
    def call_from_executor(self, callback, _max_postpone_until=None):
        """
        Call this function in the main event loop. Similar to Twisted's
        ``callFromThread``.

        :param _max_postpone_until: `None` or `time.time` value. For interal
            use. If the eventloop is saturated, consider this task to be low
            priority and postpone maximum until this timestamp. (For instance,
            repaint is done using low priority.)

            Note: In the past, this used to be a datetime.datetime instance,
                  but apparently, executing `time.time` is more efficient: it
                  does fewer system calls. (It doesn't read /etc/localtime.)
        """
