"""Implementation of basic magic functions."""


import argparse
import textwrap
import io
import sys
from pprint import pformat

from IPython.core import magic_arguments, page
from IPython.core.error import UsageError
from IPython.core.magic import Magics, magics_class, line_magic, magic_escapes
from IPython.utils.text import format_screen, dedent, indent
from IPython.testing.skipdoctest import skip_doctest
from IPython.utils.ipstruct import Struct
from warnings import warn
from logging import error


class MagicsDisplay(object):
    def __init__(self, magics_manager, ignore=None):
        self.ignore = ignore if ignore else []
        self.magics_manager = magics_manager
    
    def _lsmagic(self):
        """The main implementation of the %lsmagic"""
        mesc = magic_escapes['line']
        cesc = magic_escapes['cell']
        mman = self.magics_manager
        magics = mman.lsmagic()
        out = ['Available line magics:',
               mesc + ('  '+mesc).join(sorted([m for m,v in magics['line'].items() if (v not in self.ignore)])),
               '',
               'Available cell magics:',
               cesc + ('  '+cesc).join(sorted([m for m,v in magics['cell'].items() if (v not in self.ignore)])),
               '',
               mman.auto_status()]
        return '\n'.join(out)

    def _repr_pretty_(self, p, cycle):
        p.text(self._lsmagic())
    
    def __str__(self):
        return self._lsmagic()
    
    def _jsonable(self):
        """turn magics dict into jsonable dict of the same structure
        
        replaces object instances with their class names as strings
        """
        magic_dict = {}
        mman = self.magics_manager
        magics = mman.lsmagic()
        for key, subdict in magics.items():
            d = {}
            magic_dict[key] = d
            for name, obj in subdict.items():
                try:
                    classname = obj.__self__.__class__.__name__
                except AttributeError:
                    classname = 'Other'
                
                d[name] = classname
        return magic_dict
        
    def _repr_json_(self):
        return self._jsonable()


@magics_class
class BasicMagics(Magics):
    """Magics that provide central IPython functionality.

    These are various magics that don't fit into specific categories but that
    are all part of the base 'IPython experience'."""

    @magic_arguments.magic_arguments()
    @magic_arguments.argument(
        '-l', '--line', action='store_true',
        help="""Create a line magic alias."""
    )
    @magic_arguments.argument(
        '-c', '--cell', action='store_true',
        help="""Create a cell magic alias."""
    )
    @magic_arguments.argument(
        'name',
        help="""Name of the magic to be created."""
    )
    @magic_arguments.argument(
        'target',
        help="""Name of the existing line or cell magic."""
    )
    @magic_arguments.argument(
        '-p', '--params', default=None,
        help="""Parameters passed to the magic function."""
    )
    @line_magic
    def alias_magic(self, line=''):
        """Create an alias for an existing line or cell magic.

        Examples
        --------
        ::

          In [1]: %alias_magic t timeit
          Created `%t` as an alias for `%timeit`.
          Created `%%t` as an alias for `%%timeit`.

          In [2]: %t -n1 pass
          1 loops, best of 3: 954 ns per loop

          In [3]: %%t -n1
             ...: pass
             ...:
          1 loops, best of 3: 954 ns per loop

          In [4]: %alias_magic --cell whereami pwd
          UsageError: Cell magic function `%%pwd` not found.
          In [5]: %alias_magic --line whereami pwd
          Created `%whereami` as an alias for `%pwd`.

          In [6]: %whereami
          Out[6]: u'/home/testuser'
          
          In [7]: %alias_magic h history -p "-l 30" --line
          Created `%h` as an alias for `%history -l 30`.
        """

        args = magic_arguments.parse_argstring(self.alias_magic, line)
        shell = self.shell
        mman = self.shell.magics_manager
        escs = ''.join(magic_escapes.values())

        target = args.target.lstrip(escs)
        name = args.name.lstrip(escs)

        params = args.params
        if (params and
                ((params.startswith('"') and params.endswith('"'))
                or (params.startswith("'") and params.endswith("'")))):
            params = params[1:-1]

        # Find the requested magics.
        m_line = shell.find_magic(target, 'line')
        m_cell = shell.find_magic(target, 'cell')
        if args.line and m_line is None:
            raise UsageError('Line magic function `%s%s` not found.' %
                             (magic_escapes['line'], target))
        if args.cell and m_cell is None:
            raise UsageError('Cell magic function `%s%s` not found.' %
                             (magic_escapes['cell'], target))

        # If --line and --cell are not specified, default to the ones
        # that are available.
        if not args.line and not args.cell:
            if not m_line and not m_cell:
                raise UsageError(
                    'No line or cell magic with name `%s` found.' % target
                )
            args.line = bool(m_line)
            args.cell = bool(m_cell)

        params_str = "" if params is None else " " + params

        if args.line:
            mman.register_alias(name, target, 'line', params)
            print('Created `%s%s` as an alias for `%s%s%s`.' % (
                magic_escapes['line'], name,
                magic_escapes['line'], target, params_str))

        if args.cell:
            mman.register_alias(name, target, 'cell', params)
            print('Created `%s%s` as an alias for `%s%s%s`.' % (
                magic_escapes['cell'], name,
                magic_escapes['cell'], target, params_str))

    @line_magic
    def lsmagic(self, parameter_s=''):
        """List currently available magic functions."""
        return MagicsDisplay(self.shell.magics_manager, ignore=[self.pip])

    def _magic_docs(self, brief=False, rest=False):
        """Return docstrings from magic functions."""
        mman = self.shell.magics_manager
        docs = mman.lsmagic_docs(brief, missing='No documentation')

        if rest:
            format_string = '**%s%s**::\n\n%s\n\n'
        else:
            format_string = '%s%s:\n%s\n'

        return ''.join(
            [format_string % (magic_escapes['line'], fname,
                              indent(dedent(fndoc)))
             for fname, fndoc in sorted(docs['line'].items())]
            +
            [format_string % (magic_escapes['cell'], fname,
                              indent(dedent(fndoc)))
             for fname, fndoc in sorted(docs['cell'].items())]
        )

    @line_magic
    def magic(self, parameter_s=''):
        """Print information about the magic function system.

        Supported formats: -latex, -brief, -rest
        """

        mode = ''
        try:
            mode = parameter_s.split()[0][1:]
        except IndexError:
            pass

        brief = (mode == 'brief')
        rest = (mode == 'rest')
        magic_docs = self._magic_docs(brief, rest)

        if mode == 'latex':
            print(self.format_latex(magic_docs))
            return
        else:
            magic_docs = format_screen(magic_docs)

        out = ["""
IPython's 'magic' functions
===========================

The magic function system provides a series of functions which allow you to
control the behavior of IPython itself, plus a lot of system-type
features. There are two kinds of magics, line-oriented and cell-oriented.

Line magics are prefixed with the % character and work much like OS
command-line calls: they get as an argument the rest of the line, where
arguments are passed without parentheses or quotes.  For example, this will
time the given statement::

        %timeit range(1000)

Cell magics are prefixed with a double %%, and they are functions that get as
an argument not only the rest of the line, but also the lines below it in a
separate argument.  These magics are called with two arguments: the rest of the
call line and the body of the cell, consisting of the lines below the first.
For example::

        %%timeit x = numpy.random.randn((100, 100))
        numpy.linalg.svd(x)

will time the execution of the numpy svd routine, running the assignment of x
as part of the setup phase, which is not timed.

In a line-oriented client (the terminal or Qt console IPython), starting a new
input with %% will automatically enter cell mode, and IPython will continue
reading input until a blank line is given.  In the notebook, simply type the
whole cell as one entity, but keep in mind that the %% escape can only be at
the very start of the cell.

NOTE: If you have 'automagic' enabled (via the command line option or with the
%automagic function), you don't need to type in the % explicitly for line
magics; cell magics always require an explicit '%%' escape.  By default,
IPython ships with automagic on, so you should only rarely need the % escape.

Example: typing '%cd mydir' (without the quotes) changes your working directory
to 'mydir', if it exists.

For a list of the available magic functions, use %lsmagic. For a description
of any of them, type %magic_name?, e.g. '%cd?'.

Currently the magic system has the following functions:""",
       magic_docs,
       "Summary of magic functions (from %slsmagic):" % magic_escapes['line'],
       str(self.lsmagic()),
       ]
        page.page('\n'.join(out))


    @line_magic
    def page(self, parameter_s=''):
        """Pretty print the object and display it through a pager.

        %page [options] OBJECT

        If no object is given, use _ (last output).

        Options:

          -r: page str(object), don't pretty-print it."""

        # After a function contributed by Olivier Aubert, slightly modified.

        # Process options/args
        opts, args = self.parse_options(parameter_s, 'r')
        raw = 'r' in opts

        oname = args and args or '_'
        info = self.shell._ofind(oname)
        if info['found']:
            txt = (raw and str or pformat)( info['obj'] )
            page.page(txt)
        else:
            print('Object `%s` not found' % oname)

    @line_magic
    def profile(self, parameter_s=''):
        """DEPRECATED since IPython 2.0.

        Raise `UsageError`. To profile code use the :magic:`prun` magic.
        

        See Also
        --------
        prun : run code using the Python profiler (:magic:`prun`)
        """
        raise UsageError("The `%profile` magic has been deprecated since IPython 2.0. "
            "and removed in IPython 6.0. Please use the value of `get_ipython().profile` instead "
            "to see current profile in use. Perhaps you meant to use `%prun` to profile code?")

    @line_magic
    def pprint(self, parameter_s=''):
        """Toggle pretty printing on/off."""
        ptformatter = self.shell.display_formatter.formatters['text/plain']
        ptformatter.pprint = bool(1 - ptformatter.pprint)
        print('Pretty printing has been turned',
              ['OFF','ON'][ptformatter.pprint])

    @line_magic
    def colors(self, parameter_s=''):
        """Switch color scheme for prompts, info system and exception handlers.

        Currently implemented schemes: NoColor, Linux, LightBG.

        Color scheme names are not case-sensitive.

        Examples
        --------
        To get a plain black and white terminal::

          %colors nocolor
        """
        def color_switch_err(name):
            warn('Error changing %s color schemes.\n%s' %
                 (name, sys.exc_info()[1]), stacklevel=2)


        new_scheme = parameter_s.strip()
        if not new_scheme:
            raise UsageError(
                "%colors: you must specify a color scheme. See '%colors?'")
        # local shortcut
        shell = self.shell

        # Set shell colour scheme
        try:
            shell.colors = new_scheme
            shell.refresh_style()
        except:
            color_switch_err('shell')

        # Set exception colors
        try:
            shell.InteractiveTB.set_colors(scheme = new_scheme)
            shell.SyntaxTB.set_colors(scheme = new_scheme)
        except:
            color_switch_err('exception')

        # Set info (for 'object?') colors
        if shell.color_info:
            try:
                shell.inspector.set_active_scheme(new_scheme)
            except:
                color_switch_err('object inspector')
        else:
            shell.inspector.set_active_scheme('NoColor')

    @line_magic
    def xmode(self, parameter_s=''):
        """Switch modes for the exception handlers.

        Valid modes: Plain, Context and Verbose.

        If called without arguments, acts as a toggle."""

        def xmode_switch_err(name):
            warn('Error changing %s exception modes.\n%s' %
                 (name,sys.exc_info()[1]))

        shell = self.shell
        new_mode = parameter_s.strip().capitalize()
        try:
            shell.InteractiveTB.set_mode(mode=new_mode)
            print('Exception reporting mode:',shell.InteractiveTB.mode)
        except:
            xmode_switch_err('user')

    @line_magic
    def pip(self, args=''):
        """
        Intercept usage of ``pip`` in IPython and direct user to run command outside of IPython.
        """
        print(textwrap.dedent('''
        The following command must be run outside of the IPython shell:

            $ pip {args}

        The Python package manager (pip) can only be used from outside of IPython.
        Please reissue the `pip` command in a separate terminal or command prompt.

        See the Python documentation for more information on how to install packages:

            https://docs.python.org/3/installing/'''.format(args=args)))

    @line_magic
    def quickref(self, arg):
        """ Show a quick reference sheet """
        from IPython.core.usage import quick_reference
        qr = quick_reference + self._magic_docs(brief=True)
        page.page(qr)

    @line_magic
    def doctest_mode(self, parameter_s=''):
        """Toggle doctest mode on and off.

        This mode is intended to make IPython behave as much as possible like a
        plain Python shell, from the perspective of how its prompts, exceptions
        and output look.  This makes it easy to copy and paste parts of a
        session into doctests.  It does so by:

        - Changing the prompts to the classic ``>>>`` ones.
        - Changing the exception reporting mode to 'Plain'.
        - Disabling pretty-printing of output.

        Note that IPython also supports the pasting of code snippets that have
        leading '>>>' and '...' prompts in them.  This means that you can paste
        doctests from files or docstrings (even if they have leading
        whitespace), and the code will execute correctly.  You can then use
        '%history -t' to see the translated history; this will give you the
        input after removal of all the leading prompts and whitespace, which
        can be pasted back into an editor.

        With these features, you can switch into this mode easily whenever you
        need to do testing and changes to doctests, without having to leave
        your existing IPython session.
        """

        # Shorthands
        shell = self.shell
        meta = shell.meta
        disp_formatter = self.shell.display_formatter
        ptformatter = disp_formatter.formatters['text/plain']
        # dstore is a data store kept in the instance metadata bag to track any
        # changes we make, so we can undo them later.
        dstore = meta.setdefault('doctest_mode',Struct())
        save_dstore = dstore.setdefault

        # save a few values we'll need to recover later
        mode = save_dstore('mode',False)
        save_dstore('rc_pprint',ptformatter.pprint)
        save_dstore('xmode',shell.InteractiveTB.mode)
        save_dstore('rc_separate_out',shell.separate_out)
        save_dstore('rc_separate_out2',shell.separate_out2)
        save_dstore('rc_separate_in',shell.separate_in)
        save_dstore('rc_active_types',disp_formatter.active_types)

        if not mode:
            # turn on

            # Prompt separators like plain python
            shell.separate_in = ''
            shell.separate_out = ''
            shell.separate_out2 = ''


            ptformatter.pprint = False
            disp_formatter.active_types = ['text/plain']

            shell.magic('xmode Plain')
        else:
            # turn off
            shell.separate_in = dstore.rc_separate_in

            shell.separate_out = dstore.rc_separate_out
            shell.separate_out2 = dstore.rc_separate_out2

            ptformatter.pprint = dstore.rc_pprint
            disp_formatter.active_types = dstore.rc_active_types

            shell.magic('xmode ' + dstore.xmode)

        # mode here is the state before we switch; switch_doctest_mode takes
        # the mode we're switching to.
        shell.switch_doctest_mode(not mode)

        # Store new mode and inform
        dstore.mode = bool(not mode)
        mode_label = ['OFF','ON'][dstore.mode]
        print('Doctest mode is:', mode_label)

    @line_magic
    def gui(self, parameter_s=''):
        """Enable or disable IPython GUI event loop integration.

        %gui [GUINAME]

        This magic replaces IPython's threaded shells that were activated
        using the (pylab/wthread/etc.) command line flags.  GUI toolkits
        can now be enabled at runtime and keyboard
        interrupts should work without any problems.  The following toolkits
        are supported:  wxPython, PyQt4, PyGTK, Tk and Cocoa (OSX)::

            %gui wx      # enable wxPython event loop integration
            %gui qt4|qt  # enable PyQt4 event loop integration
            %gui qt5     # enable PyQt5 event loop integration
            %gui gtk     # enable PyGTK event loop integration
            %gui gtk3    # enable Gtk3 event loop integration
            %gui tk      # enable Tk event loop integration
            %gui osx     # enable Cocoa event loop integration
                         # (requires %matplotlib 1.1)
            %gui         # disable all event loop integration

        WARNING:  after any of these has been called you can simply create
        an application object, but DO NOT start the event loop yourself, as
        we have already handled that.
        """
        opts, arg = self.parse_options(parameter_s, '')
        if arg=='': arg = None
        try:
            return self.shell.enable_gui(arg)
        except Exception as e:
            # print simple error message, rather than traceback if we can't
            # hook up the GUI
            error(str(e))

    @skip_doctest
    @line_magic
    def precision(self, s=''):
        """Set floating point precision for pretty printing.

        Can set either integer precision or a format string.

        If numpy has been imported and precision is an int,
        numpy display precision will also be set, via ``numpy.set_printoptions``.

        If no argument is given, defaults will be restored.

        Examples
        --------
        ::

            In [1]: from math import pi

            In [2]: %precision 3
            Out[2]: u'%.3f'

            In [3]: pi
            Out[3]: 3.142

            In [4]: %precision %i
            Out[4]: u'%i'

            In [5]: pi
            Out[5]: 3

            In [6]: %precision %e
            Out[6]: u'%e'

            In [7]: pi**10
            Out[7]: 9.364805e+04

            In [8]: %precision
            Out[8]: u'%r'

            In [9]: pi**10
            Out[9]: 93648.047476082982
        """
        ptformatter = self.shell.display_formatter.formatters['text/plain']
        ptformatter.float_precision = s
        return ptformatter.float_format

    @magic_arguments.magic_arguments()
    @magic_arguments.argument(
        '-e', '--export', action='store_true', default=False,
        help=argparse.SUPPRESS
    )
    @magic_arguments.argument(
        'filename', type=str,
        help='Notebook name or filename'
    )
    @line_magic
    def notebook(self, s):
        """Export and convert IPython notebooks.

        This function can export the current IPython history to a notebook file.
        For example, to export the history to "foo.ipynb" do "%notebook foo.ipynb".

        The -e or --export flag is deprecated in IPython 5.2, and will be
        removed in the future.
        """
        args = magic_arguments.parse_argstring(self.notebook, s)

        from nbformat import write, v4

        cells = []
        hist = list(self.shell.history_manager.get_range())
        if(len(hist)<=1):
            raise ValueError('History is empty, cannot export')
        for session, execution_count, source in hist[:-1]:
            cells.append(v4.new_code_cell(
                execution_count=execution_count,
                source=source
            ))
        nb = v4.new_notebook(cells=cells)
        with io.open(args.filename, 'w', encoding='utf-8') as f:
            write(nb, f, version=4)
