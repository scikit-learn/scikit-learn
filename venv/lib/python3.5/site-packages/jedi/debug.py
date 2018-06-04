from jedi._compatibility import encoding, is_py3, u
import os
import time

def _lazy_colorama_init():
    """
    Lazily init colorama if necessary, not to screw up stdout is debug not
    enabled.

    This version of the function does nothing.
    """
    pass

_inited=False

try:
    if os.name == 'nt':
        # Does not work on Windows, as pyreadline and colorama interfere
        raise ImportError
    else:
        # Use colorama for nicer console output.
        from colorama import Fore, init
        from colorama import initialise
        def _lazy_colorama_init():
            """
            Lazily init colorama if necessary, not to screw up stdout is
            debug not enabled.

            This version of the function does init colorama.
            """
            global _inited
            if not _inited:
                # pytest resets the stream at the end - causes troubles. Since
                # after every output the stream is reset automatically we don't
                # need this.
                initialise.atexit_done = True
                try:
                    init(strip=False)
                except Exception:
                    # Colorama fails with initializing under vim and is buggy in
                    # version 0.3.6.
                    pass
            _inited = True

except ImportError:
    class Fore(object):
        RED = ''
        GREEN = ''
        YELLOW = ''
        MAGENTA = ''
        RESET = ''

NOTICE = object()
WARNING = object()
SPEED = object()

enable_speed = False
enable_warning = False
enable_notice = False

# callback, interface: level, str
debug_function = None
_debug_indent = 0
_start_time = time.time()


def reset_time():
    global _start_time, _debug_indent
    _start_time = time.time()
    _debug_indent = 0


def increase_indent(func):
    """Decorator for makin """
    def wrapper(*args, **kwargs):
        global _debug_indent
        _debug_indent += 1
        try:
            return func(*args, **kwargs)
        finally:
            _debug_indent -= 1
    return wrapper


def dbg(message, *args, **kwargs):
    """ Looks at the stack, to see if a debug message should be printed. """
    # Python 2 compatibility, because it doesn't understand default args
    color = kwargs.pop('color', 'GREEN')
    assert color

    if debug_function and enable_notice:
        i = ' ' * _debug_indent
        _lazy_colorama_init()
        debug_function(color, i + 'dbg: ' + message % tuple(u(repr(a)) for a in args))


def warning(message, *args, **kwargs):
    format = kwargs.pop('format', True)
    assert not kwargs

    if debug_function and enable_warning:
        i = ' ' * _debug_indent
        if format:
            message = message % tuple(u(repr(a)) for a in args)
        debug_function('RED', i + 'warning: ' + message)


def speed(name):
    if debug_function and enable_speed:
        now = time.time()
        i = ' ' * _debug_indent
        debug_function('YELLOW', i + 'speed: ' + '%s %s' % (name, now - _start_time))


def print_to_stdout(color, str_out):
    """
    The default debug function that prints to standard out.

    :param str color: A string that is an attribute of ``colorama.Fore``.
    """
    col = getattr(Fore, color)
    _lazy_colorama_init()
    if not is_py3:
        str_out = str_out.encode(encoding, 'replace')
    print(col + str_out + Fore.RESET)


# debug_function = print_to_stdout
