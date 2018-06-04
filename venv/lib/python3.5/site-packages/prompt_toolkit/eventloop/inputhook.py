"""
Similar to `PyOS_InputHook` of the Python API. Some eventloops can have an
inputhook to allow easy integration with other event loops.

When the eventloop of prompt-toolkit is idle, it can call such a hook. This
hook can call another eventloop that runs for a short while, for instance to
keep a graphical user interface responsive.

It's the responsibility of this hook to exit when there is input ready.
There are two ways to detect when input is ready:

- Call the `input_is_ready` method periodically. Quit when this returns `True`.

- Add the `fileno` as a watch to the external eventloop. Quit when file descriptor
  becomes readable. (But don't read from it.)

  Note that this is not the same as checking for `sys.stdin.fileno()`. The
  eventloop of prompt-toolkit allows thread-based executors, for example for
  asynchronous autocompletion. When the completion for instance is ready, we
  also want prompt-toolkit to gain control again in order to display that.

An alternative to using input hooks, is to create a custom `EventLoop` class that
controls everything.
"""
from __future__ import unicode_literals
import os
import threading
from prompt_toolkit.utils import is_windows
from .select import select_fds

__all__ = (
    'InputHookContext',
)


class InputHookContext(object):
    """
    Given as a parameter to the inputhook.
    """
    def __init__(self, inputhook):
        assert callable(inputhook)

        self.inputhook = inputhook
        self._input_is_ready = None

        self._r, self._w = os.pipe()

    def input_is_ready(self):
        """
        Return True when the input is ready.
        """
        return self._input_is_ready(wait=False)

    def fileno(self):
        """
        File descriptor that will become ready when the event loop needs to go on.
        """
        return self._r

    def call_inputhook(self, input_is_ready_func):
        """
        Call the inputhook. (Called by a prompt-toolkit eventloop.)
        """
        self._input_is_ready = input_is_ready_func

        # Start thread that activates this pipe when there is input to process.
        def thread():
            input_is_ready_func(wait=True)
            os.write(self._w, b'x')

        threading.Thread(target=thread).start()

        # Call inputhook.
        self.inputhook(self)

        # Flush the read end of the pipe.
        try:
            # Before calling 'os.read', call select.select. This is required
            # when the gevent monkey patch has been applied. 'os.read' is never
            # monkey patched and won't be cooperative, so that would block all
            # other select() calls otherwise.
            # See: http://www.gevent.org/gevent.os.html

            # Note: On Windows, this is apparently not an issue.
            #       However, if we would ever want to add a select call, it
            #       should use `windll.kernel32.WaitForMultipleObjects`,
            #       because `select.select` can't wait for a pipe on Windows.
            if not is_windows():
                select_fds([self._r], timeout=None)

            os.read(self._r, 1024)
        except OSError:
            # This happens when the window resizes and a SIGWINCH was received.
            # We get 'Error: [Errno 4] Interrupted system call'
            # Just ignore.
            pass
        self._input_is_ready = None

    def close(self):
        """
        Clean up resources.
        """
        if self._r:
            os.close(self._r)
            os.close(self._w)

        self._r = self._w = None
