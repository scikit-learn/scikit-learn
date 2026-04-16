"""
Generator for C/C++.
"""

# Serge Guelton: The licensing terms are not set in the source package, but
# pypi[1] says the software is under the MIT license, so I reproduce it here
# [1] http://pypi.python.org/pypi/cgen
#
# Copyright (C) 2008 Andreas Kloeckner
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

from textwrap import dedent
from pythran.config import cfg
from pythran.tables import pythran_ward
from pythran.spec import signatures_to_string
from pythran.utils import quote_cxxstring

__copyright__ = "Copyright (C) 2008 Andreas Kloeckner"


class Nop(object):

    def generate(self, with_semicolon=True):
        yield ''


class Declarator(object):
    def generate(self, with_semicolon=True):
        tp_lines, tp_decl = self.get_decl_pair()
        tp_lines = list(tp_lines)
        for line in tp_lines[:-1]:
            yield line
        sc = ";" if with_semicolon else ""
        if tp_decl is None:
            yield "%s%s" % (tp_lines[-1], sc)
        else:
            yield "%s %s%s" % (tp_lines[-1], tp_decl, sc)

    def get_decl_pair(self):
        """Return a tuple ``(type_lines, rhs)``.

        *type_lines* is a non-empty list of lines (most often just a
        single one) describing the type of this declarator. *rhs* is the right-
        hand side that actually contains the function/array/constness notation
        making up the bulk of the declarator syntax.
        """

    def inline(self):
        """Return the declarator as a single line."""
        tp_lines, tp_decl = self.get_decl_pair()
        tp_lines = " ".join(tp_lines)
        if tp_decl is None:
            return tp_lines
        else:
            return "%s %s" % (tp_lines, tp_decl)


class Value(Declarator):
    """A simple declarator: *typename* and *name* are given as strings."""

    def __init__(self, typename, name):
        self.typename = typename
        self.name = name

    def get_decl_pair(self):
        return [self.typename], self.name


class NestedDeclarator(Declarator):
    def __init__(self, subdecl):
        self.subdecl = subdecl

    @property
    def name(self):
        return self.subdecl.name

    def get_decl_pair(self):
        return self.subdecl.get_decl_pair()


class DeclSpecifier(NestedDeclarator):
    def __init__(self, subdecl, spec, sep=' '):
        NestedDeclarator.__init__(self, subdecl)
        self.spec = spec
        self.sep = sep

    def get_decl_pair(self):
        def add_spec(sub_it):
            it = iter(sub_it)
            try:
                yield "%s%s%s" % (self.spec, self.sep, next(it))
            except StopIteration:
                pass

            for line in it:
                yield line

        sub_tp, sub_decl = self.subdecl.get_decl_pair()
        return add_spec(sub_tp), sub_decl


class Typedef(DeclSpecifier):
    def __init__(self, subdecl):
        DeclSpecifier.__init__(self, subdecl, "typedef")


class FunctionDeclaration(NestedDeclarator):
    def __init__(self, subdecl, arg_decls, *attributes):
        NestedDeclarator.__init__(self, subdecl)
        self.inline = True
        self.arg_decls = arg_decls
        self.attributes = attributes

    def get_decl_pair(self):
        sub_tp, sub_decl = self.subdecl.get_decl_pair()
        if self.inline:
            sub_tp = ['inline'] + sub_tp
        return sub_tp, ("%s(%s) %s" % (
            sub_decl,
            ", ".join(ad.inline() for ad in self.arg_decls),
            " ".join(self.attributes)))


class Struct(Declarator):

    """
    A structure declarator.

    Attributes
    ----------
    tpname : str
        Name of the structure. (None for unnamed struct)
    fields : [Declarator]
        Content of the structure.
    inherit : str
        Parent class of current structure.
    """

    def __init__(self, tpname, fields, inherit=None):
        """Initialize the structure declarator.  """
        self.tpname = tpname
        self.fields = fields
        self.inherit = inherit

    def get_decl_pair(self):
        """ See Declarator.get_decl_pair."""
        def get_tp():
            """ Iterator generating lines for struct definition. """
            decl = "struct "
            if self.tpname is not None:
                decl += self.tpname
                if self.inherit is not None:
                    decl += " : " + self.inherit
            yield decl
            yield "{"
            for f in self.fields:
                for f_line in f.generate():
                    yield "  " + f_line
            yield "} "
        return get_tp(), ""


# template --------------------------------------------------------------------
class Template(NestedDeclarator):
    def __init__(self, template_spec, subdecl):
        super(Template, self).__init__(subdecl)
        self.template_spec = template_spec

    def generate(self, with_semicolon=False):
        yield "template <%s>" % ", ".join(self.template_spec)
        for i in self.subdecl.generate(with_semicolon):
            yield i
        if not isinstance(self.subdecl, (Template, FunctionDeclaration)):
            yield ";"


# control flow/statement stuff ------------------------------------------------
class ExceptHandler(object):
    def __init__(self, name, body, alias=None):
        self.name = name
        self.body = body
        self.alias = alias

    def generate(self):
        if self.name is None:
            yield "catch(...)"
        else:
            yield "catch (pythonic::types::%s const& %s)" % (self.name,
                                                             self.alias or '')
        for line in self.body.generate():
            yield line


class TryExcept(object):
    def __init__(self, try_, except_):
        self.try_ = try_
        self.except_ = except_

    def generate(self):
        yield "try"

        for line in self.try_.generate():
            yield line

        for exception in self.except_:
            for line in exception.generate():
                yield "  " + line


class If(object):
    def __init__(self, condition, then_, else_=None):
        self.condition = condition
        self.then_ = then_
        self.else_ = else_

    def generate(self):
        yield "if (%s)" % self.condition

        for line in self.then_.generate():
            yield line

        if self.else_ is not None:
            yield "else"
            for line in self.else_.generate():
                yield line


class Loop(object):
    def __init__(self, body):
        self.body = body

    def generate(self):
        yield self.intro_line()
        for line in self.body.generate():
            yield line


class While(Loop):
    def __init__(self, condition, body):
        super(While, self).__init__(body)
        self.condition = condition

    def intro_line(self):
        return "while (%s)" % self.condition


class For(Loop):
    def __init__(self, start, condition, update, body):
        super(For, self).__init__(body)
        self.start = start
        self.condition = condition
        self.update = update

    def intro_line(self):
        return "for (%s; %s; %s)" % (self.start, self.condition, self.update)


class AutoFor(Loop):
    def __init__(self, target, iter_, body):
        super(AutoFor, self).__init__(body)
        self.target = target
        self.iter = iter_

    def intro_line(self):
        if self.target == '_':
            return "for (PYTHRAN_UNUSED auto&& {0}: {1})".format(self.target, self.iter)
        else:
            return "for (auto&& {0}: {1})".format(self.target, self.iter)


# simple statements -----------------------------------------------------------
class Define(object):
    def __init__(self, symbol, value):
        self.symbol = symbol
        self.value = value

    def generate(self):
        yield "#define %s %s" % (self.symbol, self.value)


class Include(object):
    def __init__(self, filename, system=True):
        self.filename = filename
        self.system = system

    def generate(self):
        if self.system:
            yield "#include <%s>" % self.filename
        else:
            yield "#include \"%s\"" % self.filename


class Label(object):
    def __init__(self, label):
        self.label = label

    def generate(self):
        yield self.label + ':;'

class LineInfo(object):
    def __init__(self, filepath, lineno):
        self.filepath = filepath
        self.lineno = lineno

    def generate(self):
        if self.filename and self.lineno:
            yield '#line {} {}'.format(self.lineno, self.filepath)


class Statement(object):
    def __init__(self, text):
        self.text = text

    def generate(self):
        yield self.text + ";"


class AnnotatedStatement(object):
    def __init__(self, stmt, annotations):
        self.stmt = stmt
        self.annotations = annotations

    def generate(self):
        for directive in self.annotations:
            pragma = "#pragma " + directive.s
            yield pragma.format(*directive.deps)
        for s in self.stmt.generate():
            yield s

class StatementWithComments(object):
    def __init__(self, stmt, comment):
        self.stmt = stmt
        self.comment = comment

    def generate(self):
        yield '// {}'.format(self.comment)
        for s in self.stmt.generate():
            yield s


class InstrumentedStatement(object):
    def __init__(self, stmt,instrumentation):
        self.stmt = stmt
        self.instrumentation = instrumentation

    def generate(self):
        yield self.instrumentation
        for s in self.stmt.generate():
            yield s


class ReturnStatement(Statement):
    def generate(self):
        yield "return " + self.text + ";"


class EmptyStatement(Statement):
    def __init__(self):
        Statement.__init__(self, "")


class Assign(object):
    def __init__(self, lvalue, rvalue):
        self.lvalue = lvalue
        self.rvalue = rvalue

    def generate(self):
        yield "%s = %s;" % (self.lvalue, self.rvalue)


class Line(object):
    def __init__(self, text=""):
        self.text = text

    def generate(self):
        yield self.text


# initializers ----------------------------------------------------------------
class FunctionBody(object):
    def __init__(self, fdecl, body):
        """Initialize a function definition. *fdecl* is expected to be
        a :class:`FunctionDeclaration` instance, while *body* is a
        :class:`Block`.
        """
        self.fdecl = fdecl
        self.body = body

    def generate(self):
        for f_line in self.fdecl.generate(with_semicolon=False):
            yield f_line
        for b_line in self.body.generate():
            yield b_line


# block -----------------------------------------------------------------------
class Block(object):
    def __init__(self, contents=None):
        if contents is None:
            contents = []
        self.contents = contents

    def generate(self):
        yield "{"
        for item in self.contents:
            for item_line in item.generate():
                yield "  " + item_line
        yield "}"


class Module(Block):
    def generate(self):
        for c in self.contents:
            for line in c.generate():
                yield line


class Namespace(Block):
    def __init__(self, name, contents=None):
        Block.__init__(self, contents)
        self.name = name

    def generate(self):
        yield "namespace " + self.name
        yield "{"
        for item in self.contents:
            for item_line in item.generate():
                yield "  " + item_line
        yield "}"


# copy-pasted from codepy.bpl, which is a real mess...
# the original code was under MIT License
# cf. http://pypi.python.org/pypi/codepy
# so I reproduce it here
#
# Copyright (C) 2008 Andreas Kloeckner
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

class PythonModule(object):
    '''
    Wraps the creation of a Pythran module wrapped a Python native Module
    '''
    def __init__(self, name, docstrings, metadata):
        '''
        Builds an empty PythonModule
        '''
        self.name = name
        self.preamble = []
        self.includes = []
        self.functions = {}
        self.global_vars = []
        self.implems = []
        self.capsules = []
        self.ufuncs = {}
        self.python_implems = []
        self.wrappers = []
        self.docstrings = docstrings

        self.metadata = metadata
        moduledoc = self.docstring(self.docstrings.get(None, ""))
        self.metadata['moduledoc'] = moduledoc

    def docstring(self, doc):
        return self.splitstring(quote_cxxstring(dedent(doc)))

    def splitstring(self, doc):
        return '"{}"'.format('\\n""'.join(doc.split('\\n')))

    def add_to_preamble(self, *pa):
        self.preamble.extend(pa)

    def add_to_includes(self, *incl):
        self.includes.extend(incl)

    def add_pyfunction(self, func, name, types, signature):
        self.add_function_to(self.python_implems, func, name, types, signature)

    def add_capsule(self, func, ptrname, sig):
        self.capsules.append((ptrname, sig))
        self.implems.append(func)

    def add_ufunc(self, func, funcname, funcobject, functypes, signature):
        self.ufuncs.setdefault(funcname, []).append(
                            (funcobject, functypes, signature)
                            )
        self.implems.append(func)

    def add_function(self, func, name, types, signature):
        self.add_function_to(self.implems, func, name, types, signature)

    def add_function_to(self, to, func, name, ctypes, signature):
        """
        Add a function to be exposed. *func* is expected to be a
        :class:`cgen.FunctionBody`.

        Because a function can have several signatures exported,
        this method actually creates a wrapper for each specialization
        and a global wrapper that checks the argument types and
        runs the correct candidate, if any
        """
        to.append(func)

        args_unboxing = []  # turns PyObject to c++ object
        args_checks = []  # check if the above conversion is valid
        wrapper_name = pythran_ward + 'wrap_' + func.fdecl.name

        for i, t in enumerate(ctypes):
            args_unboxing.append('from_python<{}>(args_obj[{}])'.format(t, i))
            args_checks.append('is_convertible<{}>(args_obj[{}])'.format(t, i))
        arg_decls = func.fdecl.arg_decls[:len(ctypes)]
        keywords = "".join('"{}", '.format(arg.name) for arg in arg_decls)
        wrapper = dedent('''
            static PyObject *
            {wname}(PyObject *self, PyObject *args, PyObject *kw)
            {{
                PyObject* args_obj[{size}+1];
                {silent_warning}
                char const* keywords[] = {{{keywords} nullptr}};
                if(! PyArg_ParseTupleAndKeywords(args, kw, "{fmt}",
                                                 (char**)keywords {objs}))
                    return nullptr;
                if({checks})
                    return to_python({name}({args}));
                else {{
                    return nullptr;
                }}
            }}''')

        self.wrappers.append(
            wrapper.format(name=func.fdecl.name,
                           silent_warning= '' if ctypes else '(void)args_obj;',
                           size=len(ctypes),
                           fmt="O" * len(ctypes),
                           objs=''.join(', &args_obj[%d]' % i
                                        for i in range(len(ctypes))),
                           args=', '.join(args_unboxing),
                           checks=' && '.join(args_checks) or '1',
                           wname=wrapper_name,
                           keywords=keywords,
                           )
        )

        func_descriptor = wrapper_name, ctypes, signature
        self.functions.setdefault(name, []).append(func_descriptor)

    def add_global_var(self, name, init):
        self.global_vars.append((name, 'to_python({})'.format(init)))

    def __str__(self):
        """Generate (i.e. yield) the source code of the
        module line-by-line.
        """
        themethods = []
        theextraobjects = []
        theoverloads = []
        for vname, vinit in self.global_vars:
            theextraobjects.append(
                'PyModule_AddObject(theModule, "{0}", {1});'.format(vname, vinit))

        for fname, overloads in self.functions.items():
            tryall = []
            signatures = []
            for overload, ctypes, signature in overloads:
                try_ = dedent("""
                    if(PyObject* obj = {name}(self, args, kw))
                        return obj;
                    PyErr_Clear();
                    """.format(name=overload))
                tryall.append(try_)
                signatures.append(signature)

            candidates = signatures_to_string(fname, signatures)

            wrapper_name = pythran_ward + 'wrapall_' + fname

            candidate = dedent('''
            static PyObject *
            {wname}(PyObject *self, PyObject *args, PyObject *kw)
            {{
                return pythonic::handle_python_exception([self, args, kw]()
                -> PyObject* {{
                {tryall}
                return pythonic::python::raise_invalid_argument(
                               "{name}", {candidates}, args, kw);
                }});
            }}
            '''.format(name=fname,
                       tryall="\n".join(tryall),
                       candidates=self.splitstring(
                           candidates.replace('\n', '\\n')
                       ),
                       wname=wrapper_name))

            fdoc = self.docstring(self.docstrings.get(fname, ''))
            themethod = dedent('''{{
                "{name}",
                (PyCFunction){wname},
                METH_VARARGS | METH_KEYWORDS,
                {doc}}}'''.format(name=fname,
                                  wname=wrapper_name,
                                  doc=fdoc))
            themethods.append(themethod)
            theoverloads.append(candidate)

        for ptrname, sig in self.capsules:
            capsule = '''
            PyModule_AddObject(theModule, "{ptrname}",
                               PyCapsule_New((void*)&{ptrname}, "{sig}", NULL)
            );'''.format(ptrname=ptrname,
                         sig=sig)
            theextraobjects.append(capsule)

        for fname, overloads in self.ufuncs.items():
            fdoc = self.docstring(self.docstrings.get(fname, ''))
            funcs = []
            types = []
            for wrapper_name, wrapper_types, overload in overloads:
                funcs.append("pythonic::types::ufunc_wrapper<{}, {}>".format(
                    wrapper_name, ", ".join(wrapper_types[-1:] + wrapper_types[:-1])))
                types.extend(overload)

            ufunc = '''
            {{
            static PyUFuncGenericFunction funcs [] = {{{funcs}}};
            static char types[] = {{{types}}};
            PyModule_AddObject(
                theModule,
                "{name}",
                PyUFunc_FromFuncAndData(
                    funcs,
                    NULL,
                    types,
                    {noverloads}, {ninputs}, {noutputs},
                    PyUFunc_None,
                    "{name}",
                    {doc},
                    0
                )
            );
            }}
            '''.format(name=fname,
                       funcs=", ".join(['reinterpret_cast<PyUFuncGenericFunction>(&{})'.format(f) for f in
                                        funcs]),
                       types=", ".join(types),
                       noverloads=len(overloads),
                       ninputs=len(types) // len(funcs) - 1,
                       noutputs=1,
                       doc=fdoc
                       )
            theextraobjects.append(ufunc)

        methods = dedent('''
            static PyMethodDef Methods[] = {{
                {methods}
                {{NULL, NULL, 0, NULL}}
            }};
            '''.format(methods="".join(m + "," for m in themethods)))

        module = dedent('''
            static struct PyModuleDef moduledef = {{
              PyModuleDef_HEAD_INIT,
              "{name}",            /* m_name */
              {moduledoc},         /* m_doc */
              -1,                  /* m_size */
              Methods,             /* m_methods */
              NULL,                /* m_reload */
              NULL,                /* m_traverse */
              NULL,                /* m_clear */
              NULL,                /* m_free */
            }};
            PyMODINIT_FUNC
            PyInit_{name}(void)
            #ifndef _WIN32
            __attribute__ ((visibility("default")))
            #if defined(GNUC) && !defined(__clang__)
            __attribute__ ((externally_visible))
            #endif
            #endif
            ;
            PyMODINIT_FUNC
            PyInit_{name}(void) {{
                import_array();
                {import_umath}
                PyObject* theModule = PyModule_Create(&moduledef);
                if(! theModule)
                    return theModule;
                {freethreading}
                PyObject * theDoc = Py_BuildValue("(ss)",
                                                  "{version}",
                                                  "{hash}");
                if(! theDoc)
                    return theModule;
                PyModule_AddObject(theModule,
                                   "__pythran__",
                                   theDoc);

                {extraobjects}
                return theModule;
            }}
            '''.format(name=self.name,
                       import_umath="import_umath();" if self.ufuncs else "",
                       extraobjects='\n'.join(theextraobjects),
                       freethreading="""
                            #ifdef Py_GIL_DISABLED
                                PyUnstable_Module_SetGIL(theModule, Py_MOD_GIL_NOT_USED);
                            #endif""" if cfg.getboolean("backend", "freethreading_compatible") else "",
                       **self.metadata))

        body = (self.preamble +
                self.includes +
                self.implems +
                [Line('#ifdef ENABLE_PYTHON_MODULE')] +
                self.python_implems +
                [Line(code) for code in self.wrappers + theoverloads] +
                [Line(methods), Line(module), Line('#endif')])

        return "\n".join(Module(body).generate())


class CompilationUnit(object):

    def __init__(self, body):
        self.body = body

    def __str__(self):
        return '\n'.join('\n'.join(s.generate()) for s in self.body)
