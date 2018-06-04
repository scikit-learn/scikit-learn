# -*- coding: utf-8 -*-
"""
Pdb debugger class.

Modified from the standard pdb.Pdb class to avoid including readline, so that
the command line completion of other programs which include this isn't
damaged.

In the future, this class will be expanded with improvements over the standard
pdb.

The code in this file is mainly lifted out of cmd.py in Python 2.2, with minor
changes. Licensing should therefore be under the standard Python terms.  For
details on the PSF (Python Software Foundation) standard license, see:

https://docs.python.org/2/license.html
"""

#*****************************************************************************
#
#       This file is licensed under the PSF license.
#
#       Copyright (C) 2001 Python Software Foundation, www.python.org
#       Copyright (C) 2005-2006 Fernando Perez. <fperez@colorado.edu>
#
#
#*****************************************************************************

import bdb
import functools
import inspect
import linecache
import sys
import warnings
import re

from IPython import get_ipython
from IPython.utils import PyColorize
from IPython.utils import coloransi, py3compat
from IPython.core.excolors import exception_colors
from IPython.testing.skipdoctest import skip_doctest


prompt = 'ipdb> '

#We have to check this directly from sys.argv, config struct not yet available
from pdb import Pdb as OldPdb

# Allow the set_trace code to operate outside of an ipython instance, even if
# it does so with some limitations.  The rest of this support is implemented in
# the Tracer constructor.

def make_arrow(pad):
    """generate the leading arrow in front of traceback or debugger"""
    if pad >= 2:
        return '-'*(pad-2) + '> '
    elif pad == 1:
        return '>'
    return ''


def BdbQuit_excepthook(et, ev, tb, excepthook=None):
    """Exception hook which handles `BdbQuit` exceptions.

    All other exceptions are processed using the `excepthook`
    parameter.
    """
    warnings.warn("`BdbQuit_excepthook` is deprecated since version 5.1",
                  DeprecationWarning, stacklevel=2)
    if et==bdb.BdbQuit:
        print('Exiting Debugger.')
    elif excepthook is not None:
        excepthook(et, ev, tb)
    else:
        # Backwards compatibility. Raise deprecation warning?
        BdbQuit_excepthook.excepthook_ori(et,ev,tb)


def BdbQuit_IPython_excepthook(self,et,ev,tb,tb_offset=None):
    warnings.warn(
        "`BdbQuit_IPython_excepthook` is deprecated since version 5.1",
        DeprecationWarning, stacklevel=2)
    print('Exiting Debugger.')


class Tracer(object):
    """
    DEPRECATED

    Class for local debugging, similar to pdb.set_trace.

    Instances of this class, when called, behave like pdb.set_trace, but
    providing IPython's enhanced capabilities.

    This is implemented as a class which must be initialized in your own code
    and not as a standalone function because we need to detect at runtime
    whether IPython is already active or not.  That detection is done in the
    constructor, ensuring that this code plays nicely with a running IPython,
    while functioning acceptably (though with limitations) if outside of it.
    """

    @skip_doctest
    def __init__(self, colors=None):
        """
        DEPRECATED

        Create a local debugger instance.

        Parameters
        ----------

        colors : str, optional
            The name of the color scheme to use, it must be one of IPython's
            valid color schemes.  If not given, the function will default to
            the current IPython scheme when running inside IPython, and to
            'NoColor' otherwise.

        Examples
        --------
        ::

            from IPython.core.debugger import Tracer; debug_here = Tracer()

        Later in your code::

            debug_here()  # -> will open up the debugger at that point.

        Once the debugger activates, you can use all of its regular commands to
        step through code, set breakpoints, etc.  See the pdb documentation
        from the Python standard library for usage details.
        """
        warnings.warn("`Tracer` is deprecated since version 5.1, directly use "
                      "`IPython.core.debugger.Pdb.set_trace()`",
                      DeprecationWarning, stacklevel=2)

        ip = get_ipython()
        if ip is None:
            # Outside of ipython, we set our own exception hook manually
            sys.excepthook = functools.partial(BdbQuit_excepthook,
                                               excepthook=sys.excepthook)
            def_colors = 'NoColor'
        else:
            # In ipython, we use its custom exception handler mechanism
            def_colors = ip.colors
            ip.set_custom_exc((bdb.BdbQuit,), BdbQuit_IPython_excepthook)

        if colors is None:
            colors = def_colors

        # The stdlib debugger internally uses a modified repr from the `repr`
        # module, that limits the length of printed strings to a hardcoded
        # limit of 30 characters.  That much trimming is too aggressive, let's
        # at least raise that limit to 80 chars, which should be enough for
        # most interactive uses.
        try:
            try:
                from reprlib import aRepr  # Py 3
            except ImportError:
                from repr import aRepr  # Py 2
            aRepr.maxstring = 80
        except:
            # This is only a user-facing convenience, so any error we encounter
            # here can be warned about but can be otherwise ignored.  These
            # printouts will tell us about problems if this API changes
            import traceback
            traceback.print_exc()

        self.debugger = Pdb(colors)

    def __call__(self):
        """Starts an interactive debugger at the point where called.

        This is similar to the pdb.set_trace() function from the std lib, but
        using IPython's enhanced debugger."""

        self.debugger.set_trace(sys._getframe().f_back)


RGX_EXTRA_INDENT = re.compile('(?<=\n)\s+')


def strip_indentation(multiline_string):
    return RGX_EXTRA_INDENT.sub('', multiline_string)


def decorate_fn_with_doc(new_fn, old_fn, additional_text=""):
    """Make new_fn have old_fn's doc string. This is particularly useful
    for the ``do_...`` commands that hook into the help system.
    Adapted from from a comp.lang.python posting
    by Duncan Booth."""
    def wrapper(*args, **kw):
        return new_fn(*args, **kw)
    if old_fn.__doc__:
        wrapper.__doc__ = strip_indentation(old_fn.__doc__) + additional_text
    return wrapper


def _file_lines(fname):
    """Return the contents of a named file as a list of lines.

    This function never raises an IOError exception: if the file can't be
    read, it simply returns an empty list."""

    try:
        outfile = open(fname)
    except IOError:
        return []
    else:
        out = outfile.readlines()
        outfile.close()
        return out


class Pdb(OldPdb):
    """Modified Pdb class, does not load readline.

    for a standalone version that uses prompt_toolkit, see
    `IPython.terminal.debugger.TerminalPdb` and
    `IPython.terminal.debugger.set_trace()`
    """

    def __init__(self, color_scheme=None, completekey=None,
                 stdin=None, stdout=None, context=5):

        # Parent constructor:
        try:
            self.context = int(context)
            if self.context <= 0:
                raise ValueError("Context must be a positive integer")
        except (TypeError, ValueError):
                raise ValueError("Context must be a positive integer")

        OldPdb.__init__(self, completekey, stdin, stdout)

        # IPython changes...
        self.shell = get_ipython()

        if self.shell is None:
            save_main = sys.modules['__main__']
            # No IPython instance running, we must create one
            from IPython.terminal.interactiveshell import \
                TerminalInteractiveShell
            self.shell = TerminalInteractiveShell.instance()
            # needed by any code which calls __import__("__main__") after
            # the debugger was entered. See also #9941.
            sys.modules['__main__'] = save_main 

        if color_scheme is not None:
            warnings.warn(
                "The `color_scheme` argument is deprecated since version 5.1",
                DeprecationWarning, stacklevel=2)
        else:
            color_scheme = self.shell.colors

        self.aliases = {}

        # Create color table: we copy the default one from the traceback
        # module and add a few attributes needed for debugging
        self.color_scheme_table = exception_colors()

        # shorthands
        C = coloransi.TermColors
        cst = self.color_scheme_table

        cst['NoColor'].colors.prompt = C.NoColor
        cst['NoColor'].colors.breakpoint_enabled = C.NoColor
        cst['NoColor'].colors.breakpoint_disabled = C.NoColor

        cst['Linux'].colors.prompt = C.Green
        cst['Linux'].colors.breakpoint_enabled = C.LightRed
        cst['Linux'].colors.breakpoint_disabled = C.Red

        cst['LightBG'].colors.prompt = C.Blue
        cst['LightBG'].colors.breakpoint_enabled = C.LightRed
        cst['LightBG'].colors.breakpoint_disabled = C.Red

        cst['Neutral'].colors.prompt = C.Blue
        cst['Neutral'].colors.breakpoint_enabled = C.LightRed
        cst['Neutral'].colors.breakpoint_disabled = C.Red


        # Add a python parser so we can syntax highlight source while
        # debugging.
        self.parser = PyColorize.Parser(style=color_scheme)
        self.set_colors(color_scheme)

        # Set the prompt - the default prompt is '(Pdb)'
        self.prompt = prompt

    def set_colors(self, scheme):
        """Shorthand access to the color table scheme selector method."""
        self.color_scheme_table.set_active_scheme(scheme)
        self.parser.style = scheme

    def interaction(self, frame, traceback):
        try:
            OldPdb.interaction(self, frame, traceback)
        except KeyboardInterrupt:
            sys.stdout.write('\n' + self.shell.get_exception_only())

    def new_do_up(self, arg):
        OldPdb.do_up(self, arg)
    do_u = do_up = decorate_fn_with_doc(new_do_up, OldPdb.do_up)

    def new_do_down(self, arg):
        OldPdb.do_down(self, arg)

    do_d = do_down = decorate_fn_with_doc(new_do_down, OldPdb.do_down)

    def new_do_frame(self, arg):
        OldPdb.do_frame(self, arg)

    def new_do_quit(self, arg):

        if hasattr(self, 'old_all_completions'):
            self.shell.Completer.all_completions=self.old_all_completions

        return OldPdb.do_quit(self, arg)

    do_q = do_quit = decorate_fn_with_doc(new_do_quit, OldPdb.do_quit)

    def new_do_restart(self, arg):
        """Restart command. In the context of ipython this is exactly the same
        thing as 'quit'."""
        self.msg("Restart doesn't make sense here. Using 'quit' instead.")
        return self.do_quit(arg)

    def print_stack_trace(self, context=None):
        if context is None:
            context = self.context
        try:
            context=int(context)
            if context <= 0:
                raise ValueError("Context must be a positive integer")
        except (TypeError, ValueError):
                raise ValueError("Context must be a positive integer")
        try:
            for frame_lineno in self.stack:
                self.print_stack_entry(frame_lineno, context=context)
        except KeyboardInterrupt:
            pass

    def print_stack_entry(self,frame_lineno, prompt_prefix='\n-> ',
                          context=None):
        if context is None:
            context = self.context
        try:
            context=int(context)
            if context <= 0:
                raise ValueError("Context must be a positive integer")
        except (TypeError, ValueError):
                raise ValueError("Context must be a positive integer")
        print(self.format_stack_entry(frame_lineno, '', context))

        # vds: >>
        frame, lineno = frame_lineno
        filename = frame.f_code.co_filename
        self.shell.hooks.synchronize_with_editor(filename, lineno, 0)
        # vds: <<

    def format_stack_entry(self, frame_lineno, lprefix=': ', context=None):
        if context is None:
            context = self.context
        try:
            context=int(context)
            if context <= 0:
                print("Context must be a positive integer")
        except (TypeError, ValueError):
                print("Context must be a positive integer")
        try:
            import reprlib  # Py 3
        except ImportError:
            import repr as reprlib  # Py 2

        ret = []

        Colors = self.color_scheme_table.active_colors
        ColorsNormal = Colors.Normal
        tpl_link = u'%s%%s%s' % (Colors.filenameEm, ColorsNormal)
        tpl_call = u'%s%%s%s%%s%s' % (Colors.vName, Colors.valEm, ColorsNormal)
        tpl_line = u'%%s%s%%s %s%%s' % (Colors.lineno, ColorsNormal)
        tpl_line_em = u'%%s%s%%s %s%%s%s' % (Colors.linenoEm, Colors.line,
                                            ColorsNormal)

        frame, lineno = frame_lineno

        return_value = ''
        if '__return__' in frame.f_locals:
            rv = frame.f_locals['__return__']
            #return_value += '->'
            return_value += reprlib.repr(rv) + '\n'
        ret.append(return_value)

        #s = filename + '(' + `lineno` + ')'
        filename = self.canonic(frame.f_code.co_filename)
        link = tpl_link % py3compat.cast_unicode(filename)

        if frame.f_code.co_name:
            func = frame.f_code.co_name
        else:
            func = "<lambda>"

        call = ''
        if func != '?':
            if '__args__' in frame.f_locals:
                args = reprlib.repr(frame.f_locals['__args__'])
            else:
                args = '()'
            call = tpl_call % (func, args)

        # The level info should be generated in the same format pdb uses, to
        # avoid breaking the pdbtrack functionality of python-mode in *emacs.
        if frame is self.curframe:
            ret.append('> ')
        else:
            ret.append('  ')
        ret.append(u'%s(%s)%s\n' % (link,lineno,call))

        start = lineno - 1 - context//2
        lines = linecache.getlines(filename)
        start = min(start, len(lines) - context)
        start = max(start, 0)
        lines = lines[start : start + context]

        for i,line in enumerate(lines):
            show_arrow = (start + 1 + i == lineno)
            linetpl = (frame is self.curframe or show_arrow) \
                      and tpl_line_em \
                      or tpl_line
            ret.append(self.__format_line(linetpl, filename,
                                          start + 1 + i, line,
                                          arrow = show_arrow) )
        return ''.join(ret)

    def __format_line(self, tpl_line, filename, lineno, line, arrow = False):
        bp_mark = ""
        bp_mark_color = ""

        new_line, err = self.parser.format2(line, 'str')
        if not err:
            line = new_line

        bp = None
        if lineno in self.get_file_breaks(filename):
            bps = self.get_breaks(filename, lineno)
            bp = bps[-1]

        if bp:
            Colors = self.color_scheme_table.active_colors
            bp_mark = str(bp.number)
            bp_mark_color = Colors.breakpoint_enabled
            if not bp.enabled:
                bp_mark_color = Colors.breakpoint_disabled

        numbers_width = 7
        if arrow:
            # This is the line with the error
            pad = numbers_width - len(str(lineno)) - len(bp_mark)
            num = '%s%s' % (make_arrow(pad), str(lineno))
        else:
            num = '%*s' % (numbers_width - len(bp_mark), str(lineno))

        return tpl_line % (bp_mark_color + bp_mark, num, line)


    def print_list_lines(self, filename, first, last):
        """The printing (as opposed to the parsing part of a 'list'
        command."""
        try:
            Colors = self.color_scheme_table.active_colors
            ColorsNormal = Colors.Normal
            tpl_line = '%%s%s%%s %s%%s' % (Colors.lineno, ColorsNormal)
            tpl_line_em = '%%s%s%%s %s%%s%s' % (Colors.linenoEm, Colors.line, ColorsNormal)
            src = []
            if filename == "<string>" and hasattr(self, "_exec_filename"):
                filename = self._exec_filename

            for lineno in range(first, last+1):
                line = linecache.getline(filename, lineno)
                if not line:
                    break

                if lineno == self.curframe.f_lineno:
                    line = self.__format_line(tpl_line_em, filename, lineno, line, arrow = True)
                else:
                    line = self.__format_line(tpl_line, filename, lineno, line, arrow = False)

                src.append(line)
                self.lineno = lineno

            print(''.join(src))

        except KeyboardInterrupt:
            pass

    def do_list(self, arg):
        """Print lines of code from the current stack frame
        """
        self.lastcmd = 'list'
        last = None
        if arg:
            try:
                x = eval(arg, {}, {})
                if type(x) == type(()):
                    first, last = x
                    first = int(first)
                    last = int(last)
                    if last < first:
                        # Assume it's a count
                        last = first + last
                else:
                    first = max(1, int(x) - 5)
            except:
                print('*** Error in argument:', repr(arg))
                return
        elif self.lineno is None:
            first = max(1, self.curframe.f_lineno - 5)
        else:
            first = self.lineno + 1
        if last is None:
            last = first + 10
        self.print_list_lines(self.curframe.f_code.co_filename, first, last)

        # vds: >>
        lineno = first
        filename = self.curframe.f_code.co_filename
        self.shell.hooks.synchronize_with_editor(filename, lineno, 0)
        # vds: <<

    do_l = do_list

    def getsourcelines(self, obj):
        lines, lineno = inspect.findsource(obj)
        if inspect.isframe(obj) and obj.f_globals is obj.f_locals:
            # must be a module frame: do not try to cut a block out of it
            return lines, 1
        elif inspect.ismodule(obj):
            return lines, 1
        return inspect.getblock(lines[lineno:]), lineno+1

    def do_longlist(self, arg):
        """Print lines of code from the current stack frame.

        Shows more lines than 'list' does.
        """
        self.lastcmd = 'longlist'
        try:
            lines, lineno = self.getsourcelines(self.curframe)
        except OSError as err:
            self.error(err)
            return
        last = lineno + len(lines)
        self.print_list_lines(self.curframe.f_code.co_filename, lineno, last)
    do_ll = do_longlist

    def do_debug(self, arg):
        """debug code
        Enter a recursive debugger that steps through the code
        argument (which is an arbitrary expression or statement to be
        executed in the current environment).
        """
        sys.settrace(None)
        globals = self.curframe.f_globals
        locals = self.curframe_locals
        p = self.__class__(completekey=self.completekey,
                           stdin=self.stdin, stdout=self.stdout)
        p.use_rawinput = self.use_rawinput
        p.prompt = "(%s) " % self.prompt.strip()
        self.message("ENTERING RECURSIVE DEBUGGER")
        sys.call_tracing(p.run, (arg, globals, locals))
        self.message("LEAVING RECURSIVE DEBUGGER")
        sys.settrace(self.trace_dispatch)
        self.lastcmd = p.lastcmd

    def do_pdef(self, arg):
        """Print the call signature for any callable object.

        The debugger interface to %pdef"""
        namespaces = [('Locals', self.curframe.f_locals),
                      ('Globals', self.curframe.f_globals)]
        self.shell.find_line_magic('pdef')(arg, namespaces=namespaces)

    def do_pdoc(self, arg):
        """Print the docstring for an object.

        The debugger interface to %pdoc."""
        namespaces = [('Locals', self.curframe.f_locals),
                      ('Globals', self.curframe.f_globals)]
        self.shell.find_line_magic('pdoc')(arg, namespaces=namespaces)

    def do_pfile(self, arg):
        """Print (or run through pager) the file where an object is defined.

        The debugger interface to %pfile.
        """
        namespaces = [('Locals', self.curframe.f_locals),
                      ('Globals', self.curframe.f_globals)]
        self.shell.find_line_magic('pfile')(arg, namespaces=namespaces)

    def do_pinfo(self, arg):
        """Provide detailed information about an object.

        The debugger interface to %pinfo, i.e., obj?."""
        namespaces = [('Locals', self.curframe.f_locals),
                      ('Globals', self.curframe.f_globals)]
        self.shell.find_line_magic('pinfo')(arg, namespaces=namespaces)

    def do_pinfo2(self, arg):
        """Provide extra detailed information about an object.

        The debugger interface to %pinfo2, i.e., obj??."""
        namespaces = [('Locals', self.curframe.f_locals),
                      ('Globals', self.curframe.f_globals)]
        self.shell.find_line_magic('pinfo2')(arg, namespaces=namespaces)

    def do_psource(self, arg):
        """Print (or run through pager) the source code for an object."""
        namespaces = [('Locals', self.curframe.f_locals),
                      ('Globals', self.curframe.f_globals)]
        self.shell.find_line_magic('psource')(arg, namespaces=namespaces)

    def do_where(self, arg):
        """w(here)
        Print a stack trace, with the most recent frame at the bottom.
        An arrow indicates the "current frame", which determines the
        context of most commands. 'bt' is an alias for this command.

        Take a number as argument as an (optional) number of context line to
        print"""
        if arg:
            context = int(arg)
            self.print_stack_trace(context)
        else:
            self.print_stack_trace()

    do_w = do_where


def set_trace(frame=None):
    """
    Start debugging from `frame`.

    If frame is not specified, debugging starts from caller's frame.
    """
    Pdb().set_trace(frame or sys._getframe().f_back)
