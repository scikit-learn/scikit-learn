#
#   Errors
#

from __future__ import absolute_import

try:
    from __builtin__ import basestring as any_string_type
except ImportError:
    any_string_type = (bytes, str)

import sys
from contextlib import contextmanager

try:
    from threading import local as _threadlocal
except ImportError:
    class _threadlocal(object): pass

threadlocal = _threadlocal()

from ..Utils import open_new_file
from . import DebugFlags
from . import Options


class PyrexError(Exception):
    pass


class PyrexWarning(Exception):
    pass

class CannotSpecialize(PyrexError):
    pass

def context(position):
    source = position[0]
    assert not (isinstance(source, any_string_type)), (
        "Please replace filename strings with Scanning.FileSourceDescriptor instances %r" % source)
    try:
        F = source.get_lines()
    except UnicodeDecodeError:
        # file has an encoding problem
        s = u"[unprintable code]\n"
    else:
        s = u''.join(F[max(0, position[1]-6):position[1]])
        s = u'...\n%s%s^\n' % (s, u' '*(position[2]))
    s = u'%s\n%s%s\n' % (u'-'*60, s, u'-'*60)
    return s

def format_position(position):
    if position:
        return u"%s:%d:%d: " % (position[0].get_error_description(),
                                position[1], position[2])
    return u''

def format_error(message, position):
    if position:
        pos_str = format_position(position)
        cont = context(position)
        message = u'\nError compiling Cython file:\n%s\n%s%s' % (cont, pos_str, message or u'')
    return message

class CompileError(PyrexError):

    def __init__(self, position = None, message = u""):
        self.position = position
        self.message_only = message
        self.formatted_message = format_error(message, position)
        self.reported = False
        Exception.__init__(self, self.formatted_message)
        # Python Exception subclass pickling is broken,
        # see https://bugs.python.org/issue1692335
        self.args = (position, message)

    def __str__(self):
        return self.formatted_message

class CompileWarning(PyrexWarning):

    def __init__(self, position = None, message = ""):
        self.position = position
        Exception.__init__(self, format_position(position) + message)

class InternalError(Exception):
    # If this is ever raised, there is a bug in the compiler.

    def __init__(self, message):
        self.message_only = message
        Exception.__init__(self, u"Internal compiler error: %s"
            % message)

class AbortError(Exception):
    # Throw this to stop the compilation immediately.

    def __init__(self, message):
        self.message_only = message
        Exception.__init__(self, u"Abort error: %s" % message)

class CompilerCrash(CompileError):
    # raised when an unexpected exception occurs in a transform
    def __init__(self, pos, context, message, cause, stacktrace=None):
        if message:
            message = u'\n' + message
        else:
            message = u'\n'
        self.message_only = message
        if context:
            message = u"Compiler crash in %s%s" % (context, message)
        if stacktrace:
            import traceback
            message += (
                u'\n\nCompiler crash traceback from this point on:\n' +
                u''.join(traceback.format_tb(stacktrace)))
        if cause:
            if not stacktrace:
                message += u'\n'
            message += u'%s: %s' % (cause.__class__.__name__, cause)
        CompileError.__init__(self, pos, message)
        # Python Exception subclass pickling is broken,
        # see https://bugs.python.org/issue1692335
        self.args = (pos, context, message, cause, stacktrace)

class NoElementTreeInstalledException(PyrexError):
    """raised when the user enabled options.gdb_debug but no ElementTree
    implementation was found
    """

def open_listing_file(path, echo_to_stderr=True):
    # Begin a new error listing. If path is None, no file
    # is opened, the error counter is just reset.
    if path is not None:
        threadlocal.cython_errors_listing_file = open_new_file(path)
    else:
        threadlocal.cython_errors_listing_file = None
    if echo_to_stderr:
        threadlocal.cython_errors_echo_file = sys.stderr
    else:
        threadlocal.cython_errors_echo_file = None
    threadlocal.cython_errors_count = 0

def close_listing_file():
    if threadlocal.cython_errors_listing_file:
        threadlocal.cython_errors_listing_file.close()
        threadlocal.cython_errors_listing_file = None

def report_error(err, use_stack=True):
    error_stack = threadlocal.cython_errors_stack
    if error_stack and use_stack:
        error_stack[-1].append(err)
    else:
        # See Main.py for why dual reporting occurs. Quick fix for now.
        if err.reported: return
        err.reported = True
        try: line = u"%s\n" % err
        except UnicodeEncodeError:
            # Python <= 2.5 does this for non-ASCII Unicode exceptions
            line = format_error(getattr(err, 'message_only', "[unprintable exception message]"),
                                getattr(err, 'position', None)) + u'\n'
        listing_file = threadlocal.cython_errors_listing_file
        if listing_file:
            try: listing_file.write(line)
            except UnicodeEncodeError:
                listing_file.write(line.encode('ASCII', 'replace'))
        echo_file = threadlocal.cython_errors_echo_file
        if echo_file:
            try: echo_file.write(line)
            except UnicodeEncodeError:
                echo_file.write(line.encode('ASCII', 'replace'))
        threadlocal.cython_errors_count += 1
        if Options.fast_fail:
            raise AbortError("fatal errors")

def error(position, message):
    #print("Errors.error:", repr(position), repr(message)) ###
    if position is None:
        raise InternalError(message)
    err = CompileError(position, message)
    if DebugFlags.debug_exception_on_error: raise Exception(err)  # debug
    report_error(err)
    return err


LEVEL = 1  # warn about all errors level 1 or higher

def _write_file_encode(file, line):
    try:
        file.write(line)
    except UnicodeEncodeError:
        file.write(line.encode('ascii', 'replace'))


def performance_hint(position, message, env):
    if not env.directives['show_performance_hints']:
        return
    warn = CompileWarning(position, message)
    line = "performance hint: %s\n" % warn
    listing_file = threadlocal.cython_errors_listing_file
    if listing_file:
        _write_file_encode(listing_file, line)
    echo_file = threadlocal.cython_errors_echo_file
    if echo_file:
        _write_file_encode(echo_file, line)
    return warn


def message(position, message, level=1):
    if level < LEVEL:
        return
    warn = CompileWarning(position, message)
    line = u"note: %s\n" % warn
    listing_file = threadlocal.cython_errors_listing_file
    if listing_file:
        _write_file_encode(listing_file, line)
    echo_file = threadlocal.cython_errors_echo_file
    if echo_file:
        _write_file_encode(echo_file, line)
    return warn


def warning(position, message, level=0):
    if level < LEVEL:
        return
    if Options.warning_errors and position:
        return error(position, message)
    warn = CompileWarning(position, message)
    line = u"warning: %s\n" % warn
    listing_file = threadlocal.cython_errors_listing_file
    if listing_file:
        _write_file_encode(listing_file, line)
    echo_file = threadlocal.cython_errors_echo_file
    if echo_file:
        _write_file_encode(echo_file, line)
    return warn


def warn_once(position, message, level=0):
    if level < LEVEL:
        return
    warn_once_seen = threadlocal.cython_errors_warn_once_seen
    if message in warn_once_seen:
        return
    warn = CompileWarning(position, message)
    line = u"warning: %s\n" % warn
    listing_file = threadlocal.cython_errors_listing_file
    if listing_file:
        _write_file_encode(listing_file, line)
    echo_file = threadlocal.cython_errors_echo_file
    if echo_file:
        _write_file_encode(echo_file, line)
    warn_once_seen.add(message)
    return warn


# These functions can be used to momentarily suppress errors.

def hold_errors():
    errors = []
    threadlocal.cython_errors_stack.append(errors)
    return errors


def release_errors(ignore=False):
    held_errors = threadlocal.cython_errors_stack.pop()
    if not ignore:
        for err in held_errors:
            report_error(err)


def held_errors():
    return threadlocal.cython_errors_stack[-1]


# same as context manager:

@contextmanager
def local_errors(ignore=False):
    errors = hold_errors()
    try:
        yield errors
    finally:
        release_errors(ignore=ignore)


# Keep all global state in thread local storage to support parallel cythonisation in distutils.

def init_thread():
    threadlocal.cython_errors_count = 0
    threadlocal.cython_errors_listing_file = None
    threadlocal.cython_errors_echo_file = None
    threadlocal.cython_errors_warn_once_seen = set()
    threadlocal.cython_errors_stack = []

def reset():
    threadlocal.cython_errors_warn_once_seen.clear()
    del threadlocal.cython_errors_stack[:]

def get_errors_count():
    return threadlocal.cython_errors_count
