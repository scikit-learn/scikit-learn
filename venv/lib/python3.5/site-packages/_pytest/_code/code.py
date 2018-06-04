from __future__ import absolute_import, division, print_function
import inspect
import sys
import traceback
from inspect import CO_VARARGS, CO_VARKEYWORDS

import attr
import re
from weakref import ref
from _pytest.compat import _PY2, _PY3, PY35, safe_str

import py
builtin_repr = repr

if _PY3:
    from traceback import format_exception_only
else:
    from ._py2traceback import format_exception_only


class Code(object):
    """ wrapper around Python code objects """

    def __init__(self, rawcode):
        if not hasattr(rawcode, "co_filename"):
            rawcode = getrawcode(rawcode)
        try:
            self.filename = rawcode.co_filename
            self.firstlineno = rawcode.co_firstlineno - 1
            self.name = rawcode.co_name
        except AttributeError:
            raise TypeError("not a code object: %r" % (rawcode,))
        self.raw = rawcode

    def __eq__(self, other):
        return self.raw == other.raw

    __hash__ = None

    def __ne__(self, other):
        return not self == other

    @property
    def path(self):
        """ return a path object pointing to source code (note that it
        might not point to an actually existing file). """
        try:
            p = py.path.local(self.raw.co_filename)
            # maybe don't try this checking
            if not p.check():
                raise OSError("py.path check failed.")
        except OSError:
            # XXX maybe try harder like the weird logic
            # in the standard lib [linecache.updatecache] does?
            p = self.raw.co_filename

        return p

    @property
    def fullsource(self):
        """ return a _pytest._code.Source object for the full source file of the code
        """
        from _pytest._code import source
        full, _ = source.findsource(self.raw)
        return full

    def source(self):
        """ return a _pytest._code.Source object for the code object's source only
        """
        # return source only for that part of code
        import _pytest._code
        return _pytest._code.Source(self.raw)

    def getargs(self, var=False):
        """ return a tuple with the argument names for the code object

            if 'var' is set True also return the names of the variable and
            keyword arguments when present
        """
        # handfull shortcut for getting args
        raw = self.raw
        argcount = raw.co_argcount
        if var:
            argcount += raw.co_flags & CO_VARARGS
            argcount += raw.co_flags & CO_VARKEYWORDS
        return raw.co_varnames[:argcount]


class Frame(object):
    """Wrapper around a Python frame holding f_locals and f_globals
    in which expressions can be evaluated."""

    def __init__(self, frame):
        self.lineno = frame.f_lineno - 1
        self.f_globals = frame.f_globals
        self.f_locals = frame.f_locals
        self.raw = frame
        self.code = Code(frame.f_code)

    @property
    def statement(self):
        """ statement this frame is at """
        import _pytest._code
        if self.code.fullsource is None:
            return _pytest._code.Source("")
        return self.code.fullsource.getstatement(self.lineno)

    def eval(self, code, **vars):
        """ evaluate 'code' in the frame

            'vars' are optional additional local variables

            returns the result of the evaluation
        """
        f_locals = self.f_locals.copy()
        f_locals.update(vars)
        return eval(code, self.f_globals, f_locals)

    def exec_(self, code, **vars):
        """ exec 'code' in the frame

            'vars' are optiona; additional local variables
        """
        f_locals = self.f_locals.copy()
        f_locals.update(vars)
        py.builtin.exec_(code, self.f_globals, f_locals)

    def repr(self, object):
        """ return a 'safe' (non-recursive, one-line) string repr for 'object'
        """
        return py.io.saferepr(object)

    def is_true(self, object):
        return object

    def getargs(self, var=False):
        """ return a list of tuples (name, value) for all arguments

            if 'var' is set True also include the variable and keyword
            arguments when present
        """
        retval = []
        for arg in self.code.getargs(var):
            try:
                retval.append((arg, self.f_locals[arg]))
            except KeyError:
                pass     # this can occur when using Psyco
        return retval


class TracebackEntry(object):
    """ a single entry in a traceback """

    _repr_style = None
    exprinfo = None

    def __init__(self, rawentry, excinfo=None):
        self._excinfo = excinfo
        self._rawentry = rawentry
        self.lineno = rawentry.tb_lineno - 1

    def set_repr_style(self, mode):
        assert mode in ("short", "long")
        self._repr_style = mode

    @property
    def frame(self):
        import _pytest._code
        return _pytest._code.Frame(self._rawentry.tb_frame)

    @property
    def relline(self):
        return self.lineno - self.frame.code.firstlineno

    def __repr__(self):
        return "<TracebackEntry %s:%d>" % (self.frame.code.path, self.lineno + 1)

    @property
    def statement(self):
        """ _pytest._code.Source object for the current statement """
        source = self.frame.code.fullsource
        return source.getstatement(self.lineno)

    @property
    def path(self):
        """ path to the source code """
        return self.frame.code.path

    def getlocals(self):
        return self.frame.f_locals
    locals = property(getlocals, None, None, "locals of underlaying frame")

    def getfirstlinesource(self):
        # on Jython this firstlineno can be -1 apparently
        return max(self.frame.code.firstlineno, 0)

    def getsource(self, astcache=None):
        """ return failing source code. """
        # we use the passed in astcache to not reparse asttrees
        # within exception info printing
        from _pytest._code.source import getstatementrange_ast
        source = self.frame.code.fullsource
        if source is None:
            return None
        key = astnode = None
        if astcache is not None:
            key = self.frame.code.path
            if key is not None:
                astnode = astcache.get(key, None)
        start = self.getfirstlinesource()
        try:
            astnode, _, end = getstatementrange_ast(self.lineno, source,
                                                    astnode=astnode)
        except SyntaxError:
            end = self.lineno + 1
        else:
            if key is not None:
                astcache[key] = astnode
        return source[start:end]

    source = property(getsource)

    def ishidden(self):
        """ return True if the current frame has a var __tracebackhide__
            resolving to True

            If __tracebackhide__ is a callable, it gets called with the
            ExceptionInfo instance and can decide whether to hide the traceback.

            mostly for internal use
        """
        try:
            tbh = self.frame.f_locals['__tracebackhide__']
        except KeyError:
            try:
                tbh = self.frame.f_globals['__tracebackhide__']
            except KeyError:
                return False

        if callable(tbh):
            return tbh(None if self._excinfo is None else self._excinfo())
        else:
            return tbh

    def __str__(self):
        try:
            fn = str(self.path)
        except py.error.Error:
            fn = '???'
        name = self.frame.code.name
        try:
            line = str(self.statement).lstrip()
        except KeyboardInterrupt:
            raise
        except:  # noqa
            line = "???"
        return "  File %r:%d in %s\n  %s\n" % (fn, self.lineno + 1, name, line)

    def name(self):
        return self.frame.code.raw.co_name
    name = property(name, None, None, "co_name of underlaying code")


class Traceback(list):
    """ Traceback objects encapsulate and offer higher level
        access to Traceback entries.
    """
    Entry = TracebackEntry

    def __init__(self, tb, excinfo=None):
        """ initialize from given python traceback object and ExceptionInfo """
        self._excinfo = excinfo
        if hasattr(tb, 'tb_next'):
            def f(cur):
                while cur is not None:
                    yield self.Entry(cur, excinfo=excinfo)
                    cur = cur.tb_next
            list.__init__(self, f(tb))
        else:
            list.__init__(self, tb)

    def cut(self, path=None, lineno=None, firstlineno=None, excludepath=None):
        """ return a Traceback instance wrapping part of this Traceback

            by provding any combination of path, lineno and firstlineno, the
            first frame to start the to-be-returned traceback is determined

            this allows cutting the first part of a Traceback instance e.g.
            for formatting reasons (removing some uninteresting bits that deal
            with handling of the exception/traceback)
        """
        for x in self:
            code = x.frame.code
            codepath = code.path
            if ((path is None or codepath == path) and
                (excludepath is None or not hasattr(codepath, 'relto') or
                 not codepath.relto(excludepath)) and
                (lineno is None or x.lineno == lineno) and
                    (firstlineno is None or x.frame.code.firstlineno == firstlineno)):
                return Traceback(x._rawentry, self._excinfo)
        return self

    def __getitem__(self, key):
        val = super(Traceback, self).__getitem__(key)
        if isinstance(key, type(slice(0))):
            val = self.__class__(val)
        return val

    def filter(self, fn=lambda x: not x.ishidden()):
        """ return a Traceback instance with certain items removed

            fn is a function that gets a single argument, a TracebackEntry
            instance, and should return True when the item should be added
            to the Traceback, False when not

            by default this removes all the TracebackEntries which are hidden
            (see ishidden() above)
        """
        return Traceback(filter(fn, self), self._excinfo)

    def getcrashentry(self):
        """ return last non-hidden traceback entry that lead
        to the exception of a traceback.
        """
        for i in range(-1, -len(self) - 1, -1):
            entry = self[i]
            if not entry.ishidden():
                return entry
        return self[-1]

    def recursionindex(self):
        """ return the index of the frame/TracebackEntry where recursion
            originates if appropriate, None if no recursion occurred
        """
        cache = {}
        for i, entry in enumerate(self):
            # id for the code.raw is needed to work around
            # the strange metaprogramming in the decorator lib from pypi
            # which generates code objects that have hash/value equality
            # XXX needs a test
            key = entry.frame.code.path, id(entry.frame.code.raw), entry.lineno
            # print "checking for recursion at", key
            values = cache.setdefault(key, [])
            if values:
                f = entry.frame
                loc = f.f_locals
                for otherloc in values:
                    if f.is_true(f.eval(co_equal,
                                        __recursioncache_locals_1=loc,
                                        __recursioncache_locals_2=otherloc)):
                        return i
            values.append(entry.frame.f_locals)
        return None


co_equal = compile('__recursioncache_locals_1 == __recursioncache_locals_2',
                   '?', 'eval')


class ExceptionInfo(object):
    """ wraps sys.exc_info() objects and offers
        help for navigating the traceback.
    """
    _striptext = ''
    _assert_start_repr = "AssertionError(u\'assert " if _PY2 else "AssertionError(\'assert "

    def __init__(self, tup=None, exprinfo=None):
        import _pytest._code
        if tup is None:
            tup = sys.exc_info()
            if exprinfo is None and isinstance(tup[1], AssertionError):
                exprinfo = getattr(tup[1], 'msg', None)
                if exprinfo is None:
                    exprinfo = py.io.saferepr(tup[1])
                if exprinfo and exprinfo.startswith(self._assert_start_repr):
                    self._striptext = 'AssertionError: '
        self._excinfo = tup
        #: the exception class
        self.type = tup[0]
        #: the exception instance
        self.value = tup[1]
        #: the exception raw traceback
        self.tb = tup[2]
        #: the exception type name
        self.typename = self.type.__name__
        #: the exception traceback (_pytest._code.Traceback instance)
        self.traceback = _pytest._code.Traceback(self.tb, excinfo=ref(self))

    def __repr__(self):
        return "<ExceptionInfo %s tblen=%d>" % (self.typename, len(self.traceback))

    def exconly(self, tryshort=False):
        """ return the exception as a string

            when 'tryshort' resolves to True, and the exception is a
            _pytest._code._AssertionError, only the actual exception part of
            the exception representation is returned (so 'AssertionError: ' is
            removed from the beginning)
        """
        lines = format_exception_only(self.type, self.value)
        text = ''.join(lines)
        text = text.rstrip()
        if tryshort:
            if text.startswith(self._striptext):
                text = text[len(self._striptext):]
        return text

    def errisinstance(self, exc):
        """ return True if the exception is an instance of exc """
        return isinstance(self.value, exc)

    def _getreprcrash(self):
        exconly = self.exconly(tryshort=True)
        entry = self.traceback.getcrashentry()
        path, lineno = entry.frame.code.raw.co_filename, entry.lineno
        return ReprFileLocation(path, lineno + 1, exconly)

    def getrepr(self, showlocals=False, style="long",
                abspath=False, tbfilter=True, funcargs=False):
        """ return str()able representation of this exception info.
            showlocals: show locals per traceback entry
            style: long|short|no|native traceback style
            tbfilter: hide entries (where __tracebackhide__ is true)

            in case of style==native, tbfilter and showlocals is ignored.
        """
        if style == 'native':
            return ReprExceptionInfo(ReprTracebackNative(
                traceback.format_exception(
                    self.type,
                    self.value,
                    self.traceback[0]._rawentry,
                )), self._getreprcrash())

        fmt = FormattedExcinfo(showlocals=showlocals, style=style,
                               abspath=abspath, tbfilter=tbfilter, funcargs=funcargs)
        return fmt.repr_excinfo(self)

    def __str__(self):
        entry = self.traceback[-1]
        loc = ReprFileLocation(entry.path, entry.lineno + 1, self.exconly())
        return str(loc)

    def __unicode__(self):
        entry = self.traceback[-1]
        loc = ReprFileLocation(entry.path, entry.lineno + 1, self.exconly())
        return unicode(loc)

    def match(self, regexp):
        """
        Match the regular expression 'regexp' on the string representation of
        the exception. If it matches then True is returned (so that it is
        possible to write 'assert excinfo.match()'). If it doesn't match an
        AssertionError is raised.
        """
        __tracebackhide__ = True
        if not re.search(regexp, str(self.value)):
            assert 0, "Pattern '{0!s}' not found in '{1!s}'".format(
                regexp, self.value)
        return True


@attr.s
class FormattedExcinfo(object):
    """ presenting information about failing Functions and Generators. """
    # for traceback entries
    flow_marker = ">"
    fail_marker = "E"

    showlocals = attr.ib(default=False)
    style = attr.ib(default="long")
    abspath = attr.ib(default=True)
    tbfilter = attr.ib(default=True)
    funcargs = attr.ib(default=False)
    astcache = attr.ib(default=attr.Factory(dict), init=False, repr=False)

    def _getindent(self, source):
        # figure out indent for given source
        try:
            s = str(source.getstatement(len(source) - 1))
        except KeyboardInterrupt:
            raise
        except:  # noqa
            try:
                s = str(source[-1])
            except KeyboardInterrupt:
                raise
            except:  # noqa
                return 0
        return 4 + (len(s) - len(s.lstrip()))

    def _getentrysource(self, entry):
        source = entry.getsource(self.astcache)
        if source is not None:
            source = source.deindent()
        return source

    def _saferepr(self, obj):
        return py.io.saferepr(obj)

    def repr_args(self, entry):
        if self.funcargs:
            args = []
            for argname, argvalue in entry.frame.getargs(var=True):
                args.append((argname, self._saferepr(argvalue)))
            return ReprFuncArgs(args)

    def get_source(self, source, line_index=-1, excinfo=None, short=False):
        """ return formatted and marked up source lines. """
        import _pytest._code
        lines = []
        if source is None or line_index >= len(source.lines):
            source = _pytest._code.Source("???")
            line_index = 0
        if line_index < 0:
            line_index += len(source)
        space_prefix = "    "
        if short:
            lines.append(space_prefix + source.lines[line_index].strip())
        else:
            for line in source.lines[:line_index]:
                lines.append(space_prefix + line)
            lines.append(self.flow_marker + "   " + source.lines[line_index])
            for line in source.lines[line_index + 1:]:
                lines.append(space_prefix + line)
        if excinfo is not None:
            indent = 4 if short else self._getindent(source)
            lines.extend(self.get_exconly(excinfo, indent=indent, markall=True))
        return lines

    def get_exconly(self, excinfo, indent=4, markall=False):
        lines = []
        indent = " " * indent
        # get the real exception information out
        exlines = excinfo.exconly(tryshort=True).split('\n')
        failindent = self.fail_marker + indent[1:]
        for line in exlines:
            lines.append(failindent + line)
            if not markall:
                failindent = indent
        return lines

    def repr_locals(self, locals):
        if self.showlocals:
            lines = []
            keys = [loc for loc in locals if loc[0] != "@"]
            keys.sort()
            for name in keys:
                value = locals[name]
                if name == '__builtins__':
                    lines.append("__builtins__ = <builtins>")
                else:
                    # This formatting could all be handled by the
                    # _repr() function, which is only reprlib.Repr in
                    # disguise, so is very configurable.
                    str_repr = self._saferepr(value)
                    # if len(str_repr) < 70 or not isinstance(value,
                    #                            (list, tuple, dict)):
                    lines.append("%-10s = %s" % (name, str_repr))
                    # else:
                    #    self._line("%-10s =\\" % (name,))
                    #    # XXX
                    #    pprint.pprint(value, stream=self.excinfowriter)
            return ReprLocals(lines)

    def repr_traceback_entry(self, entry, excinfo=None):
        import _pytest._code
        source = self._getentrysource(entry)
        if source is None:
            source = _pytest._code.Source("???")
            line_index = 0
        else:
            # entry.getfirstlinesource() can be -1, should be 0 on jython
            line_index = entry.lineno - max(entry.getfirstlinesource(), 0)

        lines = []
        style = entry._repr_style
        if style is None:
            style = self.style
        if style in ("short", "long"):
            short = style == "short"
            reprargs = self.repr_args(entry) if not short else None
            s = self.get_source(source, line_index, excinfo, short=short)
            lines.extend(s)
            if short:
                message = "in %s" % (entry.name)
            else:
                message = excinfo and excinfo.typename or ""
            path = self._makepath(entry.path)
            filelocrepr = ReprFileLocation(path, entry.lineno + 1, message)
            localsrepr = None
            if not short:
                localsrepr = self.repr_locals(entry.locals)
            return ReprEntry(lines, reprargs, localsrepr, filelocrepr, style)
        if excinfo:
            lines.extend(self.get_exconly(excinfo, indent=4))
        return ReprEntry(lines, None, None, None, style)

    def _makepath(self, path):
        if not self.abspath:
            try:
                np = py.path.local().bestrelpath(path)
            except OSError:
                return path
            if len(np) < len(str(path)):
                path = np
        return path

    def repr_traceback(self, excinfo):
        traceback = excinfo.traceback
        if self.tbfilter:
            traceback = traceback.filter()

        if is_recursion_error(excinfo):
            traceback, extraline = self._truncate_recursive_traceback(traceback)
        else:
            extraline = None

        last = traceback[-1]
        entries = []
        for index, entry in enumerate(traceback):
            einfo = (last == entry) and excinfo or None
            reprentry = self.repr_traceback_entry(entry, einfo)
            entries.append(reprentry)
        return ReprTraceback(entries, extraline, style=self.style)

    def _truncate_recursive_traceback(self, traceback):
        """
        Truncate the given recursive traceback trying to find the starting point
        of the recursion.

        The detection is done by going through each traceback entry and finding the
        point in which the locals of the frame are equal to the locals of a previous frame (see ``recursionindex()``.

        Handle the situation where the recursion process might raise an exception (for example
        comparing numpy arrays using equality raises a TypeError), in which case we do our best to
        warn the user of the error and show a limited traceback.
        """
        try:
            recursionindex = traceback.recursionindex()
        except Exception as e:
            max_frames = 10
            extraline = (
                '!!! Recursion error detected, but an error occurred locating the origin of recursion.\n'
                '  The following exception happened when comparing locals in the stack frame:\n'
                '    {exc_type}: {exc_msg}\n'
                '  Displaying first and last {max_frames} stack frames out of {total}.'
            ).format(exc_type=type(e).__name__, exc_msg=safe_str(e), max_frames=max_frames, total=len(traceback))
            traceback = traceback[:max_frames] + traceback[-max_frames:]
        else:
            if recursionindex is not None:
                extraline = "!!! Recursion detected (same locals & position)"
                traceback = traceback[:recursionindex + 1]
            else:
                extraline = None

        return traceback, extraline

    def repr_excinfo(self, excinfo):
        if _PY2:
            reprtraceback = self.repr_traceback(excinfo)
            reprcrash = excinfo._getreprcrash()

            return ReprExceptionInfo(reprtraceback, reprcrash)
        else:
            repr_chain = []
            e = excinfo.value
            descr = None
            while e is not None:
                if excinfo:
                    reprtraceback = self.repr_traceback(excinfo)
                    reprcrash = excinfo._getreprcrash()
                else:
                    # fallback to native repr if the exception doesn't have a traceback:
                    # ExceptionInfo objects require a full traceback to work
                    reprtraceback = ReprTracebackNative(traceback.format_exception(type(e), e, None))
                    reprcrash = None

                repr_chain += [(reprtraceback, reprcrash, descr)]
                if e.__cause__ is not None:
                    e = e.__cause__
                    excinfo = ExceptionInfo((type(e), e, e.__traceback__)) if e.__traceback__ else None
                    descr = 'The above exception was the direct cause of the following exception:'
                elif (e.__context__ is not None and not e.__suppress_context__):
                    e = e.__context__
                    excinfo = ExceptionInfo((type(e), e, e.__traceback__)) if e.__traceback__ else None
                    descr = 'During handling of the above exception, another exception occurred:'
                else:
                    e = None
            repr_chain.reverse()
            return ExceptionChainRepr(repr_chain)


class TerminalRepr(object):
    def __str__(self):
        s = self.__unicode__()
        if _PY2:
            s = s.encode('utf-8')
        return s

    def __unicode__(self):
        # FYI this is called from pytest-xdist's serialization of exception
        # information.
        io = py.io.TextIO()
        tw = py.io.TerminalWriter(file=io)
        self.toterminal(tw)
        return io.getvalue().strip()

    def __repr__(self):
        return "<%s instance at %0x>" % (self.__class__, id(self))


class ExceptionRepr(TerminalRepr):
    def __init__(self):
        self.sections = []

    def addsection(self, name, content, sep="-"):
        self.sections.append((name, content, sep))

    def toterminal(self, tw):
        for name, content, sep in self.sections:
            tw.sep(sep, name)
            tw.line(content)


class ExceptionChainRepr(ExceptionRepr):
    def __init__(self, chain):
        super(ExceptionChainRepr, self).__init__()
        self.chain = chain
        # reprcrash and reprtraceback of the outermost (the newest) exception
        # in the chain
        self.reprtraceback = chain[-1][0]
        self.reprcrash = chain[-1][1]

    def toterminal(self, tw):
        for element in self.chain:
            element[0].toterminal(tw)
            if element[2] is not None:
                tw.line("")
                tw.line(element[2], yellow=True)
        super(ExceptionChainRepr, self).toterminal(tw)


class ReprExceptionInfo(ExceptionRepr):
    def __init__(self, reprtraceback, reprcrash):
        super(ReprExceptionInfo, self).__init__()
        self.reprtraceback = reprtraceback
        self.reprcrash = reprcrash

    def toterminal(self, tw):
        self.reprtraceback.toterminal(tw)
        super(ReprExceptionInfo, self).toterminal(tw)


class ReprTraceback(TerminalRepr):
    entrysep = "_ "

    def __init__(self, reprentries, extraline, style):
        self.reprentries = reprentries
        self.extraline = extraline
        self.style = style

    def toterminal(self, tw):
        # the entries might have different styles
        for i, entry in enumerate(self.reprentries):
            if entry.style == "long":
                tw.line("")
            entry.toterminal(tw)
            if i < len(self.reprentries) - 1:
                next_entry = self.reprentries[i + 1]
                if entry.style == "long" or \
                   entry.style == "short" and next_entry.style == "long":
                    tw.sep(self.entrysep)

        if self.extraline:
            tw.line(self.extraline)


class ReprTracebackNative(ReprTraceback):
    def __init__(self, tblines):
        self.style = "native"
        self.reprentries = [ReprEntryNative(tblines)]
        self.extraline = None


class ReprEntryNative(TerminalRepr):
    style = "native"

    def __init__(self, tblines):
        self.lines = tblines

    def toterminal(self, tw):
        tw.write("".join(self.lines))


class ReprEntry(TerminalRepr):
    localssep = "_ "

    def __init__(self, lines, reprfuncargs, reprlocals, filelocrepr, style):
        self.lines = lines
        self.reprfuncargs = reprfuncargs
        self.reprlocals = reprlocals
        self.reprfileloc = filelocrepr
        self.style = style

    def toterminal(self, tw):
        if self.style == "short":
            self.reprfileloc.toterminal(tw)
            for line in self.lines:
                red = line.startswith("E   ")
                tw.line(line, bold=True, red=red)
            # tw.line("")
            return
        if self.reprfuncargs:
            self.reprfuncargs.toterminal(tw)
        for line in self.lines:
            red = line.startswith("E   ")
            tw.line(line, bold=True, red=red)
        if self.reprlocals:
            # tw.sep(self.localssep, "Locals")
            tw.line("")
            self.reprlocals.toterminal(tw)
        if self.reprfileloc:
            if self.lines:
                tw.line("")
            self.reprfileloc.toterminal(tw)

    def __str__(self):
        return "%s\n%s\n%s" % ("\n".join(self.lines),
                               self.reprlocals,
                               self.reprfileloc)


class ReprFileLocation(TerminalRepr):
    def __init__(self, path, lineno, message):
        self.path = str(path)
        self.lineno = lineno
        self.message = message

    def toterminal(self, tw):
        # filename and lineno output for each entry,
        # using an output format that most editors unterstand
        msg = self.message
        i = msg.find("\n")
        if i != -1:
            msg = msg[:i]
        tw.write(self.path, bold=True, red=True)
        tw.line(":%s: %s" % (self.lineno, msg))


class ReprLocals(TerminalRepr):
    def __init__(self, lines):
        self.lines = lines

    def toterminal(self, tw):
        for line in self.lines:
            tw.line(line)


class ReprFuncArgs(TerminalRepr):
    def __init__(self, args):
        self.args = args

    def toterminal(self, tw):
        if self.args:
            linesofar = ""
            for name, value in self.args:
                ns = "%s = %s" % (safe_str(name), safe_str(value))
                if len(ns) + len(linesofar) + 2 > tw.fullwidth:
                    if linesofar:
                        tw.line(linesofar)
                    linesofar = ns
                else:
                    if linesofar:
                        linesofar += ", " + ns
                    else:
                        linesofar = ns
            if linesofar:
                tw.line(linesofar)
            tw.line("")


def getrawcode(obj, trycall=True):
    """ return code object for given function. """
    try:
        return obj.__code__
    except AttributeError:
        obj = getattr(obj, 'im_func', obj)
        obj = getattr(obj, 'func_code', obj)
        obj = getattr(obj, 'f_code', obj)
        obj = getattr(obj, '__code__', obj)
        if trycall and not hasattr(obj, 'co_firstlineno'):
            if hasattr(obj, '__call__') and not inspect.isclass(obj):
                x = getrawcode(obj.__call__, trycall=False)
                if hasattr(x, 'co_firstlineno'):
                    return x
        return obj


if PY35:  # RecursionError introduced in 3.5
    def is_recursion_error(excinfo):
        return excinfo.errisinstance(RecursionError)  # noqa
else:
    def is_recursion_error(excinfo):
        if not excinfo.errisinstance(RuntimeError):
            return False
        try:
            return "maximum recursion depth exceeded" in str(excinfo.value)
        except UnicodeError:
            return False
