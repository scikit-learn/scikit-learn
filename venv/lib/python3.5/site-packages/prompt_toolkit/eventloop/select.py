"""
Selectors for the Posix event loop.
"""
from __future__ import unicode_literals, absolute_import
import sys
import abc
import errno
import select
import six

__all__ = (
    'AutoSelector',
    'PollSelector',
    'SelectSelector',
    'Selector',
    'fd_to_int',
)

def fd_to_int(fd):
    assert isinstance(fd, int) or hasattr(fd, 'fileno')

    if isinstance(fd, int):
        return fd
    else:
        return fd.fileno()


class Selector(six.with_metaclass(abc.ABCMeta, object)):
    @abc.abstractmethod
    def register(self, fd):
        assert isinstance(fd, int)

    @abc.abstractmethod
    def unregister(self, fd):
        assert isinstance(fd, int)

    @abc.abstractmethod
    def select(self, timeout):
        pass

    @abc.abstractmethod
    def close(self):
        pass


class AutoSelector(Selector):
    def __init__(self):
        self._fds = []

        self._select_selector = SelectSelector()
        self._selectors = [self._select_selector]

        # When 'select.poll' exists, create a PollSelector.
        if hasattr(select, 'poll'):
            self._poll_selector = PollSelector()
            self._selectors.append(self._poll_selector)
        else:
            self._poll_selector = None

        # Use of the 'select' module, that was introduced in Python3.4. We don't
        # use it before 3.5 however, because this is the point where this module
        # retries interrupted system calls.
        if sys.version_info >= (3, 5):
            self._py3_selector = Python3Selector()
            self._selectors.append(self._py3_selector)
        else:
            self._py3_selector = None

    def register(self, fd):
        assert isinstance(fd, int)

        self._fds.append(fd)

        for sel in self._selectors:
            sel.register(fd)

    def unregister(self, fd):
        assert isinstance(fd, int)

        self._fds.remove(fd)

        for sel in self._selectors:
            sel.unregister(fd)

    def select(self, timeout):
        # Try Python 3 selector first.
        if self._py3_selector:
            try:
                return self._py3_selector.select(timeout)
            except PermissionError:
                # We had a situation (in pypager) where epoll raised a
                # PermissionError when a local file descriptor was registered,
                # however poll and select worked fine. So, in that case, just
                # try using select below.
                pass

        try:
            # Prefer 'select.select', if we don't have much file descriptors.
            # This is more universal.
            return self._select_selector.select(timeout)
        except ValueError:
            # When we have more than 1024 open file descriptors, we'll always
            # get a "ValueError: filedescriptor out of range in select()" for
            # 'select'. In this case, try, using 'poll' instead.
            if self._poll_selector is not None:
                return self._poll_selector.select(timeout)
            else:
                raise

    def close(self):
        for sel in self._selectors:
            sel.close()


class Python3Selector(Selector):
    """
    Use of the Python3 'selectors' module.

    NOTE: Only use on Python 3.5 or newer!
    """
    def __init__(self):
        assert sys.version_info >= (3, 5)

        import selectors  # Inline import: Python3 only!
        self._sel = selectors.DefaultSelector()

    def register(self, fd):
        assert isinstance(fd, int)
        import selectors  # Inline import: Python3 only!
        self._sel.register(fd, selectors.EVENT_READ, None)

    def unregister(self, fd):
        assert isinstance(fd, int)
        self._sel.unregister(fd)

    def select(self, timeout):
        events = self._sel.select(timeout=timeout)
        return [key.fileobj for key, mask in events]

    def close(self):
        self._sel.close()


class PollSelector(Selector):
    def __init__(self):
        self._poll = select.poll()

    def register(self, fd):
        assert isinstance(fd, int)
        self._poll.register(fd, select.POLLIN)

    def unregister(self, fd):
        assert isinstance(fd, int)

    def select(self, timeout):
        tuples = self._poll.poll(timeout)  # Returns (fd, event) tuples.
        return [t[0] for t in tuples]

    def close(self):
        pass  # XXX


class SelectSelector(Selector):
    """
    Wrapper around select.select.

    When the SIGWINCH signal is handled, other system calls, like select
    are aborted in Python. This wrapper will retry the system call.
    """
    def __init__(self):
        self._fds = []

    def register(self, fd):
        self._fds.append(fd)

    def unregister(self, fd):
        self._fds.remove(fd)

    def select(self, timeout):
        while True:
            try:
                return select.select(self._fds, [], [], timeout)[0]
            except select.error as e:
                # Retry select call when EINTR
                if e.args and e.args[0] == errno.EINTR:
                    continue
                else:
                    raise

    def close(self):
        pass


def select_fds(read_fds, timeout, selector=AutoSelector):
    """
    Wait for a list of file descriptors (`read_fds`) to become ready for
    reading. This chooses the most appropriate select-tool for use in
    prompt-toolkit.
    """
    # Map to ensure that we return the objects that were passed in originally.
    # Whether they are a fd integer or an object that has a fileno().
    # (The 'poll' implementation for instance, returns always integers.)
    fd_map = dict((fd_to_int(fd), fd) for fd in read_fds)

    # Wait, using selector.
    sel = selector()
    try:
        for fd in read_fds:
            sel.register(fd)

        result = sel.select(timeout)

        if result is not None:
            return [fd_map[fd_to_int(fd)] for fd in result]
    finally:
        sel.close()
