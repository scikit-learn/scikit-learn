"""Rewrite assertion AST to produce nice error messages"""
from __future__ import absolute_import, division, print_function
import ast
import errno
import itertools
import imp
import marshal
import os
import re
import six
import struct
import sys
import types

import py
from _pytest.assertion import util


# pytest caches rewritten pycs in __pycache__.
if hasattr(imp, "get_tag"):
    PYTEST_TAG = imp.get_tag() + "-PYTEST"
else:
    if hasattr(sys, "pypy_version_info"):
        impl = "pypy"
    elif sys.platform == "java":
        impl = "jython"
    else:
        impl = "cpython"
    ver = sys.version_info
    PYTEST_TAG = "%s-%s%s-PYTEST" % (impl, ver[0], ver[1])
    del ver, impl

PYC_EXT = ".py" + (__debug__ and "c" or "o")
PYC_TAIL = "." + PYTEST_TAG + PYC_EXT

ASCII_IS_DEFAULT_ENCODING = sys.version_info[0] < 3

if sys.version_info >= (3, 5):
    ast_Call = ast.Call
else:
    def ast_Call(a, b, c):
        return ast.Call(a, b, c, None, None)


class AssertionRewritingHook(object):
    """PEP302 Import hook which rewrites asserts."""

    def __init__(self, config):
        self.config = config
        self.fnpats = config.getini("python_files")
        self.session = None
        self.modules = {}
        self._rewritten_names = set()
        self._register_with_pkg_resources()
        self._must_rewrite = set()

    def set_session(self, session):
        self.session = session

    def find_module(self, name, path=None):
        state = self.config._assertstate
        state.trace("find_module called for: %s" % name)
        names = name.rsplit(".", 1)
        lastname = names[-1]
        pth = None
        if path is not None:
            # Starting with Python 3.3, path is a _NamespacePath(), which
            # causes problems if not converted to list.
            path = list(path)
            if len(path) == 1:
                pth = path[0]
        if pth is None:
            try:
                fd, fn, desc = imp.find_module(lastname, path)
            except ImportError:
                return None
            if fd is not None:
                fd.close()
            tp = desc[2]
            if tp == imp.PY_COMPILED:
                if hasattr(imp, "source_from_cache"):
                    try:
                        fn = imp.source_from_cache(fn)
                    except ValueError:
                        # Python 3 doesn't like orphaned but still-importable
                        # .pyc files.
                        fn = fn[:-1]
                else:
                    fn = fn[:-1]
            elif tp != imp.PY_SOURCE:
                # Don't know what this is.
                return None
        else:
            fn = os.path.join(pth, name.rpartition(".")[2] + ".py")

        fn_pypath = py.path.local(fn)
        if not self._should_rewrite(name, fn_pypath, state):
            return None

        self._rewritten_names.add(name)

        # The requested module looks like a test file, so rewrite it. This is
        # the most magical part of the process: load the source, rewrite the
        # asserts, and load the rewritten source. We also cache the rewritten
        # module code in a special pyc. We must be aware of the possibility of
        # concurrent pytest processes rewriting and loading pycs. To avoid
        # tricky race conditions, we maintain the following invariant: The
        # cached pyc is always a complete, valid pyc. Operations on it must be
        # atomic. POSIX's atomic rename comes in handy.
        write = not sys.dont_write_bytecode
        cache_dir = os.path.join(fn_pypath.dirname, "__pycache__")
        if write:
            try:
                os.mkdir(cache_dir)
            except OSError:
                e = sys.exc_info()[1].errno
                if e == errno.EEXIST:
                    # Either the __pycache__ directory already exists (the
                    # common case) or it's blocked by a non-dir node. In the
                    # latter case, we'll ignore it in _write_pyc.
                    pass
                elif e in [errno.ENOENT, errno.ENOTDIR]:
                    # One of the path components was not a directory, likely
                    # because we're in a zip file.
                    write = False
                elif e in [errno.EACCES, errno.EROFS, errno.EPERM]:
                    state.trace("read only directory: %r" % fn_pypath.dirname)
                    write = False
                else:
                    raise
        cache_name = fn_pypath.basename[:-3] + PYC_TAIL
        pyc = os.path.join(cache_dir, cache_name)
        # Notice that even if we're in a read-only directory, I'm going
        # to check for a cached pyc. This may not be optimal...
        co = _read_pyc(fn_pypath, pyc, state.trace)
        if co is None:
            state.trace("rewriting %r" % (fn,))
            source_stat, co = _rewrite_test(self.config, fn_pypath)
            if co is None:
                # Probably a SyntaxError in the test.
                return None
            if write:
                _make_rewritten_pyc(state, source_stat, pyc, co)
        else:
            state.trace("found cached rewritten pyc for %r" % (fn,))
        self.modules[name] = co, pyc
        return self

    def _should_rewrite(self, name, fn_pypath, state):
        # always rewrite conftest files
        fn = str(fn_pypath)
        if fn_pypath.basename == 'conftest.py':
            state.trace("rewriting conftest file: %r" % (fn,))
            return True

        if self.session is not None:
            if self.session.isinitpath(fn):
                state.trace("matched test file (was specified on cmdline): %r" %
                            (fn,))
                return True

        # modules not passed explicitly on the command line are only
        # rewritten if they match the naming convention for test files
        for pat in self.fnpats:
            if fn_pypath.fnmatch(pat):
                state.trace("matched test file %r" % (fn,))
                return True

        for marked in self._must_rewrite:
            if name == marked or name.startswith(marked + '.'):
                state.trace("matched marked file %r (from %r)" % (name, marked))
                return True

        return False

    def mark_rewrite(self, *names):
        """Mark import names as needing to be rewritten.

        The named module or package as well as any nested modules will
        be rewritten on import.
        """
        already_imported = (set(names)
                            .intersection(sys.modules)
                            .difference(self._rewritten_names))
        for name in already_imported:
            if not AssertionRewriter.is_rewrite_disabled(
                    sys.modules[name].__doc__ or ""):
                self._warn_already_imported(name)
        self._must_rewrite.update(names)

    def _warn_already_imported(self, name):
        self.config.warn(
            'P1',
            'Module already imported so cannot be rewritten: %s' % name)

    def load_module(self, name):
        # If there is an existing module object named 'fullname' in
        # sys.modules, the loader must use that existing module. (Otherwise,
        # the reload() builtin will not work correctly.)
        if name in sys.modules:
            return sys.modules[name]

        co, pyc = self.modules.pop(name)
        # I wish I could just call imp.load_compiled here, but __file__ has to
        # be set properly. In Python 3.2+, this all would be handled correctly
        # by load_compiled.
        mod = sys.modules[name] = imp.new_module(name)
        try:
            mod.__file__ = co.co_filename
            # Normally, this attribute is 3.2+.
            mod.__cached__ = pyc
            mod.__loader__ = self
            py.builtin.exec_(co, mod.__dict__)
        except:  # noqa
            if name in sys.modules:
                del sys.modules[name]
            raise
        return sys.modules[name]

    def is_package(self, name):
        try:
            fd, fn, desc = imp.find_module(name)
        except ImportError:
            return False
        if fd is not None:
            fd.close()
        tp = desc[2]
        return tp == imp.PKG_DIRECTORY

    @classmethod
    def _register_with_pkg_resources(cls):
        """
        Ensure package resources can be loaded from this loader. May be called
        multiple times, as the operation is idempotent.
        """
        try:
            import pkg_resources
            # access an attribute in case a deferred importer is present
            pkg_resources.__name__
        except ImportError:
            return

        # Since pytest tests are always located in the file system, the
        #  DefaultProvider is appropriate.
        pkg_resources.register_loader_type(cls, pkg_resources.DefaultProvider)

    def get_data(self, pathname):
        """Optional PEP302 get_data API.
        """
        with open(pathname, 'rb') as f:
            return f.read()


def _write_pyc(state, co, source_stat, pyc):
    # Technically, we don't have to have the same pyc format as
    # (C)Python, since these "pycs" should never be seen by builtin
    # import. However, there's little reason deviate, and I hope
    # sometime to be able to use imp.load_compiled to load them. (See
    # the comment in load_module above.)
    try:
        fp = open(pyc, "wb")
    except IOError:
        err = sys.exc_info()[1].errno
        state.trace("error writing pyc file at %s: errno=%s" % (pyc, err))
        # we ignore any failure to write the cache file
        # there are many reasons, permission-denied, __pycache__ being a
        # file etc.
        return False
    try:
        fp.write(imp.get_magic())
        mtime = int(source_stat.mtime)
        size = source_stat.size & 0xFFFFFFFF
        fp.write(struct.pack("<ll", mtime, size))
        marshal.dump(co, fp)
    finally:
        fp.close()
    return True


RN = "\r\n".encode("utf-8")
N = "\n".encode("utf-8")

cookie_re = re.compile(r"^[ \t\f]*#.*coding[:=][ \t]*[-\w.]+")
BOM_UTF8 = '\xef\xbb\xbf'


def _rewrite_test(config, fn):
    """Try to read and rewrite *fn* and return the code object."""
    state = config._assertstate
    try:
        stat = fn.stat()
        source = fn.read("rb")
    except EnvironmentError:
        return None, None
    if ASCII_IS_DEFAULT_ENCODING:
        # ASCII is the default encoding in Python 2. Without a coding
        # declaration, Python 2 will complain about any bytes in the file
        # outside the ASCII range. Sadly, this behavior does not extend to
        # compile() or ast.parse(), which prefer to interpret the bytes as
        # latin-1. (At least they properly handle explicit coding cookies.) To
        # preserve this error behavior, we could force ast.parse() to use ASCII
        # as the encoding by inserting a coding cookie. Unfortunately, that
        # messes up line numbers. Thus, we have to check ourselves if anything
        # is outside the ASCII range in the case no encoding is explicitly
        # declared. For more context, see issue #269. Yay for Python 3 which
        # gets this right.
        end1 = source.find("\n")
        end2 = source.find("\n", end1 + 1)
        if (not source.startswith(BOM_UTF8) and
            cookie_re.match(source[0:end1]) is None and
                cookie_re.match(source[end1 + 1:end2]) is None):
            if hasattr(state, "_indecode"):
                # encodings imported us again, so don't rewrite.
                return None, None
            state._indecode = True
            try:
                try:
                    source.decode("ascii")
                except UnicodeDecodeError:
                    # Let it fail in real import.
                    return None, None
            finally:
                del state._indecode
    try:
        tree = ast.parse(source)
    except SyntaxError:
        # Let this pop up again in the real import.
        state.trace("failed to parse: %r" % (fn,))
        return None, None
    rewrite_asserts(tree, fn, config)
    try:
        co = compile(tree, fn.strpath, "exec", dont_inherit=True)
    except SyntaxError:
        # It's possible that this error is from some bug in the
        # assertion rewriting, but I don't know of a fast way to tell.
        state.trace("failed to compile: %r" % (fn,))
        return None, None
    return stat, co


def _make_rewritten_pyc(state, source_stat, pyc, co):
    """Try to dump rewritten code to *pyc*."""
    if sys.platform.startswith("win"):
        # Windows grants exclusive access to open files and doesn't have atomic
        # rename, so just write into the final file.
        _write_pyc(state, co, source_stat, pyc)
    else:
        # When not on windows, assume rename is atomic. Dump the code object
        # into a file specific to this process and atomically replace it.
        proc_pyc = pyc + "." + str(os.getpid())
        if _write_pyc(state, co, source_stat, proc_pyc):
            os.rename(proc_pyc, pyc)


def _read_pyc(source, pyc, trace=lambda x: None):
    """Possibly read a pytest pyc containing rewritten code.

    Return rewritten code if successful or None if not.
    """
    try:
        fp = open(pyc, "rb")
    except IOError:
        return None
    with fp:
        try:
            mtime = int(source.mtime())
            size = source.size()
            data = fp.read(12)
        except EnvironmentError as e:
            trace('_read_pyc(%s): EnvironmentError %s' % (source, e))
            return None
        # Check for invalid or out of date pyc file.
        if (len(data) != 12 or data[:4] != imp.get_magic() or
                struct.unpack("<ll", data[4:]) != (mtime, size)):
            trace('_read_pyc(%s): invalid or out of date pyc' % source)
            return None
        try:
            co = marshal.load(fp)
        except Exception as e:
            trace('_read_pyc(%s): marshal.load error %s' % (source, e))
            return None
        if not isinstance(co, types.CodeType):
            trace('_read_pyc(%s): not a code object' % source)
            return None
        return co


def rewrite_asserts(mod, module_path=None, config=None):
    """Rewrite the assert statements in mod."""
    AssertionRewriter(module_path, config).run(mod)


def _saferepr(obj):
    """Get a safe repr of an object for assertion error messages.

    The assertion formatting (util.format_explanation()) requires
    newlines to be escaped since they are a special character for it.
    Normally assertion.util.format_explanation() does this but for a
    custom repr it is possible to contain one of the special escape
    sequences, especially '\n{' and '\n}' are likely to be present in
    JSON reprs.

    """
    repr = py.io.saferepr(obj)
    if isinstance(repr, six.text_type):
        t = six.text_type
    else:
        t = six.binary_type
    return repr.replace(t("\n"), t("\\n"))


from _pytest.assertion.util import format_explanation as _format_explanation  # noqa


def _format_assertmsg(obj):
    """Format the custom assertion message given.

    For strings this simply replaces newlines with '\n~' so that
    util.format_explanation() will preserve them instead of escaping
    newlines.  For other objects py.io.saferepr() is used first.

    """
    # reprlib appears to have a bug which means that if a string
    # contains a newline it gets escaped, however if an object has a
    # .__repr__() which contains newlines it does not get escaped.
    # However in either case we want to preserve the newline.
    if isinstance(obj, six.text_type) or isinstance(obj, six.binary_type):
        s = obj
        is_repr = False
    else:
        s = py.io.saferepr(obj)
        is_repr = True
    if isinstance(s, six.text_type):
        t = six.text_type
    else:
        t = six.binary_type
    s = s.replace(t("\n"), t("\n~")).replace(t("%"), t("%%"))
    if is_repr:
        s = s.replace(t("\\n"), t("\n~"))
    return s


def _should_repr_global_name(obj):
    return not hasattr(obj, "__name__") and not callable(obj)


def _format_boolop(explanations, is_or):
    explanation = "(" + (is_or and " or " or " and ").join(explanations) + ")"
    if isinstance(explanation, six.text_type):
        t = six.text_type
    else:
        t = six.binary_type
    return explanation.replace(t('%'), t('%%'))


def _call_reprcompare(ops, results, expls, each_obj):
    for i, res, expl in zip(range(len(ops)), results, expls):
        try:
            done = not res
        except Exception:
            done = True
        if done:
            break
    if util._reprcompare is not None:
        custom = util._reprcompare(ops[i], each_obj[i], each_obj[i + 1])
        if custom is not None:
            return custom
    return expl


unary_map = {
    ast.Not: "not %s",
    ast.Invert: "~%s",
    ast.USub: "-%s",
    ast.UAdd: "+%s"
}

binop_map = {
    ast.BitOr: "|",
    ast.BitXor: "^",
    ast.BitAnd: "&",
    ast.LShift: "<<",
    ast.RShift: ">>",
    ast.Add: "+",
    ast.Sub: "-",
    ast.Mult: "*",
    ast.Div: "/",
    ast.FloorDiv: "//",
    ast.Mod: "%%",  # escaped for string formatting
    ast.Eq: "==",
    ast.NotEq: "!=",
    ast.Lt: "<",
    ast.LtE: "<=",
    ast.Gt: ">",
    ast.GtE: ">=",
    ast.Pow: "**",
    ast.Is: "is",
    ast.IsNot: "is not",
    ast.In: "in",
    ast.NotIn: "not in"
}
# Python 3.5+ compatibility
try:
    binop_map[ast.MatMult] = "@"
except AttributeError:
    pass

# Python 3.4+ compatibility
if hasattr(ast, "NameConstant"):
    _NameConstant = ast.NameConstant
else:
    def _NameConstant(c):
        return ast.Name(str(c), ast.Load())


def set_location(node, lineno, col_offset):
    """Set node location information recursively."""
    def _fix(node, lineno, col_offset):
        if "lineno" in node._attributes:
            node.lineno = lineno
        if "col_offset" in node._attributes:
            node.col_offset = col_offset
        for child in ast.iter_child_nodes(node):
            _fix(child, lineno, col_offset)
    _fix(node, lineno, col_offset)
    return node


class AssertionRewriter(ast.NodeVisitor):
    """Assertion rewriting implementation.

    The main entrypoint is to call .run() with an ast.Module instance,
    this will then find all the assert statements and rewrite them to
    provide intermediate values and a detailed assertion error.  See
    http://pybites.blogspot.be/2011/07/behind-scenes-of-pytests-new-assertion.html
    for an overview of how this works.

    The entry point here is .run() which will iterate over all the
    statements in an ast.Module and for each ast.Assert statement it
    finds call .visit() with it.  Then .visit_Assert() takes over and
    is responsible for creating new ast statements to replace the
    original assert statement: it rewrites the test of an assertion
    to provide intermediate values and replace it with an if statement
    which raises an assertion error with a detailed explanation in
    case the expression is false.

    For this .visit_Assert() uses the visitor pattern to visit all the
    AST nodes of the ast.Assert.test field, each visit call returning
    an AST node and the corresponding explanation string.  During this
    state is kept in several instance attributes:

    :statements: All the AST statements which will replace the assert
       statement.

    :variables: This is populated by .variable() with each variable
       used by the statements so that they can all be set to None at
       the end of the statements.

    :variable_counter: Counter to create new unique variables needed
       by statements.  Variables are created using .variable() and
       have the form of "@py_assert0".

    :on_failure: The AST statements which will be executed if the
       assertion test fails.  This is the code which will construct
       the failure message and raises the AssertionError.

    :explanation_specifiers: A dict filled by .explanation_param()
       with %-formatting placeholders and their corresponding
       expressions to use in the building of an assertion message.
       This is used by .pop_format_context() to build a message.

    :stack: A stack of the explanation_specifiers dicts maintained by
       .push_format_context() and .pop_format_context() which allows
       to build another %-formatted string while already building one.

    This state is reset on every new assert statement visited and used
    by the other visitors.

    """

    def __init__(self, module_path, config):
        super(AssertionRewriter, self).__init__()
        self.module_path = module_path
        self.config = config

    def run(self, mod):
        """Find all assert statements in *mod* and rewrite them."""
        if not mod.body:
            # Nothing to do.
            return
        # Insert some special imports at the top of the module but after any
        # docstrings and __future__ imports.
        aliases = [ast.alias(py.builtin.builtins.__name__, "@py_builtins"),
                   ast.alias("_pytest.assertion.rewrite", "@pytest_ar")]
        doc = getattr(mod, "docstring", None)
        expect_docstring = doc is None
        if doc is not None and self.is_rewrite_disabled(doc):
            return
        pos = 0
        lineno = 1
        for item in mod.body:
            if (expect_docstring and isinstance(item, ast.Expr) and
                    isinstance(item.value, ast.Str)):
                doc = item.value.s
                if self.is_rewrite_disabled(doc):
                    return
                expect_docstring = False
            elif (not isinstance(item, ast.ImportFrom) or item.level > 0 or
                  item.module != "__future__"):
                lineno = item.lineno
                break
            pos += 1
        else:
            lineno = item.lineno
        imports = [ast.Import([alias], lineno=lineno, col_offset=0)
                   for alias in aliases]
        mod.body[pos:pos] = imports
        # Collect asserts.
        nodes = [mod]
        while nodes:
            node = nodes.pop()
            for name, field in ast.iter_fields(node):
                if isinstance(field, list):
                    new = []
                    for i, child in enumerate(field):
                        if isinstance(child, ast.Assert):
                            # Transform assert.
                            new.extend(self.visit(child))
                        else:
                            new.append(child)
                            if isinstance(child, ast.AST):
                                nodes.append(child)
                    setattr(node, name, new)
                elif (isinstance(field, ast.AST) and
                      # Don't recurse into expressions as they can't contain
                      # asserts.
                      not isinstance(field, ast.expr)):
                    nodes.append(field)

    @staticmethod
    def is_rewrite_disabled(docstring):
        return "PYTEST_DONT_REWRITE" in docstring

    def variable(self):
        """Get a new variable."""
        # Use a character invalid in python identifiers to avoid clashing.
        name = "@py_assert" + str(next(self.variable_counter))
        self.variables.append(name)
        return name

    def assign(self, expr):
        """Give *expr* a name."""
        name = self.variable()
        self.statements.append(ast.Assign([ast.Name(name, ast.Store())], expr))
        return ast.Name(name, ast.Load())

    def display(self, expr):
        """Call py.io.saferepr on the expression."""
        return self.helper("saferepr", expr)

    def helper(self, name, *args):
        """Call a helper in this module."""
        py_name = ast.Name("@pytest_ar", ast.Load())
        attr = ast.Attribute(py_name, "_" + name, ast.Load())
        return ast_Call(attr, list(args), [])

    def builtin(self, name):
        """Return the builtin called *name*."""
        builtin_name = ast.Name("@py_builtins", ast.Load())
        return ast.Attribute(builtin_name, name, ast.Load())

    def explanation_param(self, expr):
        """Return a new named %-formatting placeholder for expr.

        This creates a %-formatting placeholder for expr in the
        current formatting context, e.g. ``%(py0)s``.  The placeholder
        and expr are placed in the current format context so that it
        can be used on the next call to .pop_format_context().

        """
        specifier = "py" + str(next(self.variable_counter))
        self.explanation_specifiers[specifier] = expr
        return "%(" + specifier + ")s"

    def push_format_context(self):
        """Create a new formatting context.

        The format context is used for when an explanation wants to
        have a variable value formatted in the assertion message.  In
        this case the value required can be added using
        .explanation_param().  Finally .pop_format_context() is used
        to format a string of %-formatted values as added by
        .explanation_param().

        """
        self.explanation_specifiers = {}
        self.stack.append(self.explanation_specifiers)

    def pop_format_context(self, expl_expr):
        """Format the %-formatted string with current format context.

        The expl_expr should be an ast.Str instance constructed from
        the %-placeholders created by .explanation_param().  This will
        add the required code to format said string to .on_failure and
        return the ast.Name instance of the formatted string.

        """
        current = self.stack.pop()
        if self.stack:
            self.explanation_specifiers = self.stack[-1]
        keys = [ast.Str(key) for key in current.keys()]
        format_dict = ast.Dict(keys, list(current.values()))
        form = ast.BinOp(expl_expr, ast.Mod(), format_dict)
        name = "@py_format" + str(next(self.variable_counter))
        self.on_failure.append(ast.Assign([ast.Name(name, ast.Store())], form))
        return ast.Name(name, ast.Load())

    def generic_visit(self, node):
        """Handle expressions we don't have custom code for."""
        assert isinstance(node, ast.expr)
        res = self.assign(node)
        return res, self.explanation_param(self.display(res))

    def visit_Assert(self, assert_):
        """Return the AST statements to replace the ast.Assert instance.

        This rewrites the test of an assertion to provide
        intermediate values and replace it with an if statement which
        raises an assertion error with a detailed explanation in case
        the expression is false.

        """
        if isinstance(assert_.test, ast.Tuple) and self.config is not None:
            fslocation = (self.module_path, assert_.lineno)
            self.config.warn('R1', 'assertion is always true, perhaps '
                             'remove parentheses?', fslocation=fslocation)
        self.statements = []
        self.variables = []
        self.variable_counter = itertools.count()
        self.stack = []
        self.on_failure = []
        self.push_format_context()
        # Rewrite assert into a bunch of statements.
        top_condition, explanation = self.visit(assert_.test)
        # Create failure message.
        body = self.on_failure
        negation = ast.UnaryOp(ast.Not(), top_condition)
        self.statements.append(ast.If(negation, body, []))
        if assert_.msg:
            assertmsg = self.helper('format_assertmsg', assert_.msg)
            explanation = "\n>assert " + explanation
        else:
            assertmsg = ast.Str("")
            explanation = "assert " + explanation
        template = ast.BinOp(assertmsg, ast.Add(), ast.Str(explanation))
        msg = self.pop_format_context(template)
        fmt = self.helper("format_explanation", msg)
        err_name = ast.Name("AssertionError", ast.Load())
        exc = ast_Call(err_name, [fmt], [])
        if sys.version_info[0] >= 3:
            raise_ = ast.Raise(exc, None)
        else:
            raise_ = ast.Raise(exc, None, None)
        body.append(raise_)
        # Clear temporary variables by setting them to None.
        if self.variables:
            variables = [ast.Name(name, ast.Store())
                         for name in self.variables]
            clear = ast.Assign(variables, _NameConstant(None))
            self.statements.append(clear)
        # Fix line numbers.
        for stmt in self.statements:
            set_location(stmt, assert_.lineno, assert_.col_offset)
        return self.statements

    def visit_Name(self, name):
        # Display the repr of the name if it's a local variable or
        # _should_repr_global_name() thinks it's acceptable.
        locs = ast_Call(self.builtin("locals"), [], [])
        inlocs = ast.Compare(ast.Str(name.id), [ast.In()], [locs])
        dorepr = self.helper("should_repr_global_name", name)
        test = ast.BoolOp(ast.Or(), [inlocs, dorepr])
        expr = ast.IfExp(test, self.display(name), ast.Str(name.id))
        return name, self.explanation_param(expr)

    def visit_BoolOp(self, boolop):
        res_var = self.variable()
        expl_list = self.assign(ast.List([], ast.Load()))
        app = ast.Attribute(expl_list, "append", ast.Load())
        is_or = int(isinstance(boolop.op, ast.Or))
        body = save = self.statements
        fail_save = self.on_failure
        levels = len(boolop.values) - 1
        self.push_format_context()
        # Process each operand, short-circuting if needed.
        for i, v in enumerate(boolop.values):
            if i:
                fail_inner = []
                # cond is set in a prior loop iteration below
                self.on_failure.append(ast.If(cond, fail_inner, []))  # noqa
                self.on_failure = fail_inner
            self.push_format_context()
            res, expl = self.visit(v)
            body.append(ast.Assign([ast.Name(res_var, ast.Store())], res))
            expl_format = self.pop_format_context(ast.Str(expl))
            call = ast_Call(app, [expl_format], [])
            self.on_failure.append(ast.Expr(call))
            if i < levels:
                cond = res
                if is_or:
                    cond = ast.UnaryOp(ast.Not(), cond)
                inner = []
                self.statements.append(ast.If(cond, inner, []))
                self.statements = body = inner
        self.statements = save
        self.on_failure = fail_save
        expl_template = self.helper("format_boolop", expl_list, ast.Num(is_or))
        expl = self.pop_format_context(expl_template)
        return ast.Name(res_var, ast.Load()), self.explanation_param(expl)

    def visit_UnaryOp(self, unary):
        pattern = unary_map[unary.op.__class__]
        operand_res, operand_expl = self.visit(unary.operand)
        res = self.assign(ast.UnaryOp(unary.op, operand_res))
        return res, pattern % (operand_expl,)

    def visit_BinOp(self, binop):
        symbol = binop_map[binop.op.__class__]
        left_expr, left_expl = self.visit(binop.left)
        right_expr, right_expl = self.visit(binop.right)
        explanation = "(%s %s %s)" % (left_expl, symbol, right_expl)
        res = self.assign(ast.BinOp(left_expr, binop.op, right_expr))
        return res, explanation

    def visit_Call_35(self, call):
        """
        visit `ast.Call` nodes on Python3.5 and after
        """
        new_func, func_expl = self.visit(call.func)
        arg_expls = []
        new_args = []
        new_kwargs = []
        for arg in call.args:
            res, expl = self.visit(arg)
            arg_expls.append(expl)
            new_args.append(res)
        for keyword in call.keywords:
            res, expl = self.visit(keyword.value)
            new_kwargs.append(ast.keyword(keyword.arg, res))
            if keyword.arg:
                arg_expls.append(keyword.arg + "=" + expl)
            else:  # **args have `arg` keywords with an .arg of None
                arg_expls.append("**" + expl)

        expl = "%s(%s)" % (func_expl, ', '.join(arg_expls))
        new_call = ast.Call(new_func, new_args, new_kwargs)
        res = self.assign(new_call)
        res_expl = self.explanation_param(self.display(res))
        outer_expl = "%s\n{%s = %s\n}" % (res_expl, res_expl, expl)
        return res, outer_expl

    def visit_Starred(self, starred):
        # From Python 3.5, a Starred node can appear in a function call
        res, expl = self.visit(starred.value)
        return starred, '*' + expl

    def visit_Call_legacy(self, call):
        """
        visit `ast.Call nodes on 3.4 and below`
        """
        new_func, func_expl = self.visit(call.func)
        arg_expls = []
        new_args = []
        new_kwargs = []
        new_star = new_kwarg = None
        for arg in call.args:
            res, expl = self.visit(arg)
            new_args.append(res)
            arg_expls.append(expl)
        for keyword in call.keywords:
            res, expl = self.visit(keyword.value)
            new_kwargs.append(ast.keyword(keyword.arg, res))
            arg_expls.append(keyword.arg + "=" + expl)
        if call.starargs:
            new_star, expl = self.visit(call.starargs)
            arg_expls.append("*" + expl)
        if call.kwargs:
            new_kwarg, expl = self.visit(call.kwargs)
            arg_expls.append("**" + expl)
        expl = "%s(%s)" % (func_expl, ', '.join(arg_expls))
        new_call = ast.Call(new_func, new_args, new_kwargs,
                            new_star, new_kwarg)
        res = self.assign(new_call)
        res_expl = self.explanation_param(self.display(res))
        outer_expl = "%s\n{%s = %s\n}" % (res_expl, res_expl, expl)
        return res, outer_expl

    # ast.Call signature changed on 3.5,
    # conditionally change  which methods is named
    # visit_Call depending on Python version
    if sys.version_info >= (3, 5):
        visit_Call = visit_Call_35
    else:
        visit_Call = visit_Call_legacy

    def visit_Attribute(self, attr):
        if not isinstance(attr.ctx, ast.Load):
            return self.generic_visit(attr)
        value, value_expl = self.visit(attr.value)
        res = self.assign(ast.Attribute(value, attr.attr, ast.Load()))
        res_expl = self.explanation_param(self.display(res))
        pat = "%s\n{%s = %s.%s\n}"
        expl = pat % (res_expl, res_expl, value_expl, attr.attr)
        return res, expl

    def visit_Compare(self, comp):
        self.push_format_context()
        left_res, left_expl = self.visit(comp.left)
        if isinstance(comp.left, (ast.Compare, ast.BoolOp)):
            left_expl = "({0})".format(left_expl)
        res_variables = [self.variable() for i in range(len(comp.ops))]
        load_names = [ast.Name(v, ast.Load()) for v in res_variables]
        store_names = [ast.Name(v, ast.Store()) for v in res_variables]
        it = zip(range(len(comp.ops)), comp.ops, comp.comparators)
        expls = []
        syms = []
        results = [left_res]
        for i, op, next_operand in it:
            next_res, next_expl = self.visit(next_operand)
            if isinstance(next_operand, (ast.Compare, ast.BoolOp)):
                next_expl = "({0})".format(next_expl)
            results.append(next_res)
            sym = binop_map[op.__class__]
            syms.append(ast.Str(sym))
            expl = "%s %s %s" % (left_expl, sym, next_expl)
            expls.append(ast.Str(expl))
            res_expr = ast.Compare(left_res, [op], [next_res])
            self.statements.append(ast.Assign([store_names[i]], res_expr))
            left_res, left_expl = next_res, next_expl
        # Use pytest.assertion.util._reprcompare if that's available.
        expl_call = self.helper("call_reprcompare",
                                ast.Tuple(syms, ast.Load()),
                                ast.Tuple(load_names, ast.Load()),
                                ast.Tuple(expls, ast.Load()),
                                ast.Tuple(results, ast.Load()))
        if len(comp.ops) > 1:
            res = ast.BoolOp(ast.And(), load_names)
        else:
            res = load_names[0]
        return res, self.explanation_param(self.pop_format_context(expl_call))
