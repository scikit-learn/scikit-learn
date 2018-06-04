from __future__ import unicode_literals
import inspect
import os
import signal
import sys
import threading
import weakref

from wcwidth import wcwidth
from six.moves import range


__all__ = (
    'Event',
    'DummyContext',
    'get_cwidth',
    'suspend_to_background_supported',
    'is_conemu_ansi',
    'is_windows',
    'in_main_thread',
    'take_using_weights',
    'test_callable_args',
)


class Event(object):
    """
    Simple event to which event handlers can be attached. For instance::

        class Cls:
            def __init__(self):
                # Define event. The first parameter is the sender.
                self.event = Event(self)

        obj = Cls()

        def handler(sender):
            pass

        # Add event handler by using the += operator.
        obj.event += handler

        # Fire event.
        obj.event()
    """
    def __init__(self, sender, handler=None):
        self.sender = sender
        self._handlers = []

        if handler is not None:
            self += handler

    def __call__(self):
        " Fire event. "
        for handler in self._handlers:
            handler(self.sender)

    def fire(self):
        " Alias for just calling the event. "
        self()

    def __iadd__(self, handler):
        """
        Add another handler to this callback.
        (Handler should be a callable that takes exactly one parameter: the
        sender object.)
        """
        # Test handler.
        assert callable(handler)
        if not test_callable_args(handler, [None]):
            raise TypeError("%r doesn't take exactly one argument." % handler)

        # Add to list of event handlers.
        self._handlers.append(handler)
        return self

    def __isub__(self, handler):
        """
        Remove a handler from this callback.
        """
        self._handlers.remove(handler)
        return self


# Cache of signatures. Improves the performance of `test_callable_args`.
_signatures_cache = weakref.WeakKeyDictionary()


def test_callable_args(func, args):
    """
    Return True when this function can be called with the given arguments.
    """
    assert isinstance(args, (list, tuple))
    signature = getattr(inspect, 'signature', None)

    if signature is not None:
        # For Python 3, use inspect.signature.
        try:
            sig = _signatures_cache[func]
        except KeyError:
            sig = signature(func)
            _signatures_cache[func] = sig

        try:
            sig.bind(*args)
        except TypeError:
            return False
        else:
            return True
    else:
        # For older Python versions, fall back to using getargspec.
        spec = inspect.getargspec(func)

        # Drop the 'self'
        def drop_self(spec):
            args, varargs, varkw, defaults = spec
            if args[0:1] == ['self']:
                args = args[1:]
            return inspect.ArgSpec(args, varargs, varkw, defaults)

        spec = drop_self(spec)

        # When taking *args, always return True.
        if spec.varargs is not None:
            return True

        # Test whether the given amount of args is between the min and max
        # accepted argument counts.
        return len(spec.args) - len(spec.defaults or []) <= len(args) <= len(spec.args)


class DummyContext(object):
    """
    (contextlib.nested is not available on Py3)
    """
    def __enter__(self):
        pass

    def __exit__(self, *a):
        pass


class _CharSizesCache(dict):
    """
    Cache for wcwidth sizes.
    """
    def __missing__(self, string):
        # Note: We use the `max(0, ...` because some non printable control
        #       characters, like e.g. Ctrl-underscore get a -1 wcwidth value.
        #       It can be possible that these characters end up in the input
        #       text.
        if len(string) == 1:
            result = max(0, wcwidth(string))
        else:
            result = sum(max(0, wcwidth(c)) for c in string)

        # Cache for short strings.
        # (It's hard to tell what we can consider short...)
        if len(string) < 256:
            self[string] = result

        return result


_CHAR_SIZES_CACHE = _CharSizesCache()


def get_cwidth(string):
    """
    Return width of a string. Wrapper around ``wcwidth``.
    """
    return _CHAR_SIZES_CACHE[string]


def suspend_to_background_supported():
    """
    Returns `True` when the Python implementation supports
    suspend-to-background. This is typically `False' on Windows systems.
    """
    return hasattr(signal, 'SIGTSTP')


def is_windows():
    """
    True when we are using Windows.
    """
    return sys.platform.startswith('win')  # E.g. 'win32', not 'darwin' or 'linux2'


def is_conemu_ansi():
    """
    True when the ConEmu Windows console is used.
    """
    return is_windows() and os.environ.get('ConEmuANSI', 'OFF') == 'ON'


def in_main_thread():
    """
    True when the current thread is the main thread.
    """
    return threading.current_thread().__class__.__name__ == '_MainThread'


def take_using_weights(items, weights):
    """
    Generator that keeps yielding items from the items list, in proportion to
    their weight. For instance::

        # Getting the first 70 items from this generator should have yielded 10
        # times A, 20 times B and 40 times C, all distributed equally..
        take_using_weights(['A', 'B', 'C'], [5, 10, 20])

    :param items: List of items to take from.
    :param weights: Integers representing the weight. (Numbers have to be
                    integers, not floats.)
    """
    assert isinstance(items, list)
    assert isinstance(weights, list)
    assert all(isinstance(i, int) for i in weights)
    assert len(items) == len(weights)
    assert len(items) > 0

    already_taken = [0 for i in items]
    item_count = len(items)
    max_weight = max(weights)

    i = 0
    while True:
        # Each iteration of this loop, we fill up until by (total_weight/max_weight).
        adding = True
        while adding:
            adding = False

            for item_i, item, weight in zip(range(item_count), items, weights):
                if already_taken[item_i] < i * weight / float(max_weight):
                    yield item
                    already_taken[item_i] += 1
                    adding = True

        i += 1
