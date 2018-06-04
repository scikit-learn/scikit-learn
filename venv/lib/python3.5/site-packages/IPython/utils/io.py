# encoding: utf-8
"""
IO related utilities.
"""

# Copyright (c) IPython Development Team.
# Distributed under the terms of the Modified BSD License.



import atexit
import os
import sys
import tempfile
import warnings
from warnings import warn

from IPython.utils.decorators import undoc
from .capture import CapturedIO, capture_output
from .py3compat import input

@undoc
class IOStream:

    def __init__(self, stream, fallback=None):
        warn('IOStream is deprecated since IPython 5.0, use sys.{stdin,stdout,stderr} instead',
             DeprecationWarning, stacklevel=2)
        if not hasattr(stream,'write') or not hasattr(stream,'flush'):
            if fallback is not None:
                stream = fallback
            else:
                raise ValueError("fallback required, but not specified")
        self.stream = stream
        self._swrite = stream.write

        # clone all methods not overridden:
        def clone(meth):
            return not hasattr(self, meth) and not meth.startswith('_')
        for meth in filter(clone, dir(stream)):
            try:
                val = getattr(stream, meth)
            except AttributeError:
                pass
            else:
                setattr(self, meth, val)

    def __repr__(self):
        cls = self.__class__
        tpl = '{mod}.{cls}({args})'
        return tpl.format(mod=cls.__module__, cls=cls.__name__, args=self.stream)

    def write(self,data):
        warn('IOStream is deprecated since IPython 5.0, use sys.{stdin,stdout,stderr} instead',
             DeprecationWarning, stacklevel=2)
        try:
            self._swrite(data)
        except:
            try:
                # print handles some unicode issues which may trip a plain
                # write() call.  Emulate write() by using an empty end
                # argument.
                print(data, end='', file=self.stream)
            except:
                # if we get here, something is seriously broken.
                print('ERROR - failed to write data to stream:', self.stream,
                      file=sys.stderr)

    def writelines(self, lines):
        warn('IOStream is deprecated since IPython 5.0, use sys.{stdin,stdout,stderr} instead',
             DeprecationWarning, stacklevel=2)
        if isinstance(lines, str):
            lines = [lines]
        for line in lines:
            self.write(line)

    # This class used to have a writeln method, but regular files and streams
    # in Python don't have this method. We need to keep this completely
    # compatible so we removed it.

    @property
    def closed(self):
        return self.stream.closed

    def close(self):
        pass

# setup stdin/stdout/stderr to sys.stdin/sys.stdout/sys.stderr
devnull = open(os.devnull, 'w')
atexit.register(devnull.close)

# io.std* are deprecated, but don't show our own deprecation warnings
# during initialization of the deprecated API.
with warnings.catch_warnings():
    warnings.simplefilter('ignore', DeprecationWarning)
    stdin = IOStream(sys.stdin, fallback=devnull)
    stdout = IOStream(sys.stdout, fallback=devnull)
    stderr = IOStream(sys.stderr, fallback=devnull)

class Tee(object):
    """A class to duplicate an output stream to stdout/err.

    This works in a manner very similar to the Unix 'tee' command.

    When the object is closed or deleted, it closes the original file given to
    it for duplication.
    """
    # Inspired by:
    # http://mail.python.org/pipermail/python-list/2007-May/442737.html

    def __init__(self, file_or_name, mode="w", channel='stdout'):
        """Construct a new Tee object.

        Parameters
        ----------
        file_or_name : filename or open filehandle (writable)
          File that will be duplicated

        mode : optional, valid mode for open().
          If a filename was give, open with this mode.

        channel : str, one of ['stdout', 'stderr']
        """
        if channel not in ['stdout', 'stderr']:
            raise ValueError('Invalid channel spec %s' % channel)

        if hasattr(file_or_name, 'write') and hasattr(file_or_name, 'seek'):
            self.file = file_or_name
        else:
            self.file = open(file_or_name, mode)
        self.channel = channel
        self.ostream = getattr(sys, channel)
        setattr(sys, channel, self)
        self._closed = False

    def close(self):
        """Close the file and restore the channel."""
        self.flush()
        setattr(sys, self.channel, self.ostream)
        self.file.close()
        self._closed = True

    def write(self, data):
        """Write data to both channels."""
        self.file.write(data)
        self.ostream.write(data)
        self.ostream.flush()

    def flush(self):
        """Flush both channels."""
        self.file.flush()
        self.ostream.flush()

    def __del__(self):
        if not self._closed:
            self.close()


def ask_yes_no(prompt, default=None, interrupt=None):
    """Asks a question and returns a boolean (y/n) answer.

    If default is given (one of 'y','n'), it is used if the user input is
    empty. If interrupt is given (one of 'y','n'), it is used if the user
    presses Ctrl-C. Otherwise the question is repeated until an answer is
    given.

    An EOF is treated as the default answer.  If there is no default, an
    exception is raised to prevent infinite loops.

    Valid answers are: y/yes/n/no (match is not case sensitive)."""

    answers = {'y':True,'n':False,'yes':True,'no':False}
    ans = None
    while ans not in answers.keys():
        try:
            ans = input(prompt+' ').lower()
            if not ans:  # response was an empty string
                ans = default
        except KeyboardInterrupt:
            if interrupt:
                ans = interrupt
            print("\r")
        except EOFError:
            if default in answers.keys():
                ans = default
                print()
            else:
                raise

    return answers[ans]


def temp_pyfile(src, ext='.py'):
    """Make a temporary python file, return filename and filehandle.

    Parameters
    ----------
    src : string or list of strings (no need for ending newlines if list)
      Source code to be written to the file.

    ext : optional, string
      Extension for the generated file.

    Returns
    -------
    (filename, open filehandle)
      It is the caller's responsibility to close the open file and unlink it.
    """
    fname = tempfile.mkstemp(ext)[1]
    f = open(fname,'w')
    f.write(src)
    f.flush()
    return fname, f

def atomic_writing(*args, **kwargs):
    """DEPRECATED: moved to notebook.services.contents.fileio"""
    warn("IPython.utils.io.atomic_writing has moved to notebook.services.contents.fileio", stacklevel=2)
    from notebook.services.contents.fileio import atomic_writing
    return atomic_writing(*args, **kwargs)

def raw_print(*args, **kw):
    """Raw print to sys.__stdout__, otherwise identical interface to print()."""

    print(*args, sep=kw.get('sep', ' '), end=kw.get('end', '\n'),
          file=sys.__stdout__)
    sys.__stdout__.flush()


def raw_print_err(*args, **kw):
    """Raw print to sys.__stderr__, otherwise identical interface to print()."""

    print(*args, sep=kw.get('sep', ' '), end=kw.get('end', '\n'),
          file=sys.__stderr__)
    sys.__stderr__.flush()


# Short aliases for quick debugging, do NOT use these in production code.
rprint = raw_print
rprinte = raw_print_err


def unicode_std_stream(stream='stdout'):
    """DEPRECATED, moved to nbconvert.utils.io"""
    warn("IPython.utils.io.unicode_std_stream has moved to nbconvert.utils.io", stacklevel=2)
    from nbconvert.utils.io import unicode_std_stream
    return unicode_std_stream(stream)
