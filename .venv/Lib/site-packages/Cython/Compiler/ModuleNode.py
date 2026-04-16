#
#   Module parse tree node
#


import cython
cython.declare(Naming=object, Options=object, PyrexTypes=object, TypeSlots=object,
               error=object, warning=object, py_object_type=object, UtilityCode=object,
               EncodedString=object, re=object)

from collections import defaultdict
import json
import operator
import os
import pathlib
import re
import sys
from typing import Sequence

from .PyrexTypes import CPtrType
from . import Future
from . import Annotate
from . import Code
from . import Naming
from . import Nodes
from . import Options
from . import TypeSlots
from . import PyrexTypes
from . import Pythran

from .Errors import error, warning, CompileError, format_position
from .PyrexTypes import py_object_type, get_all_subtypes
from ..Utils import open_new_file, replace_suffix, decode_filename, build_hex_version, is_cython_generated_file
from .Code import UtilityCode, IncludeCode, TempitaUtilityCode
from .StringEncoding import EncodedString, bytes_literal, encoded_string_or_bytes_literal
from .Pythran import has_np_pythran


def replace_suffix_encoded(path, newsuf):
    # calls replace suffix and returns a EncodedString or BytesLiteral with the encoding set
    newpath = replace_suffix(path, newsuf)
    return as_encoded_filename(newpath)

def as_encoded_filename(path):
    # wraps the path with either EncodedString or BytesLiteral (depending on its input type)
    # and sets the encoding to the file system encoding
    return encoded_string_or_bytes_literal(path, sys.getfilesystemencoding())


def check_c_declarations_pxd(module_node):
    module_node.scope.check_c_classes_pxd()
    return module_node


def check_c_declarations(module_node):
    module_node.scope.check_c_classes()
    module_node.scope.check_c_functions()
    return module_node


def generate_c_code_config(env, options):
    if Options.annotate or options.annotate:
        emit_linenums = False
    else:
        emit_linenums = options.emit_linenums

    if hasattr(options, "emit_code_comments"):
        print('Warning: option emit_code_comments is deprecated. '
              'Instead, use compiler directive emit_code_comments.')

    return Code.CCodeConfig(
        emit_linenums=emit_linenums,
        emit_code_comments=env.directives['emit_code_comments'],
        c_line_in_traceback=options.c_line_in_traceback)

# The code required to generate one comparison from another.
# The keys are (from, to).
# The comparison operator always goes first, with equality possibly second.
# The first value specifies if the comparison is inverted. The second is the
# logic op to use, and the third is if the equality is inverted or not.
TOTAL_ORDERING = {
    # a > b from (not a < b) and (a != b)
    ('__lt__', '__gt__'): (True, '&&', True),
    # a <= b from (a < b) or (a == b)
    ('__lt__', '__le__'): (False, '||', False),
    # a >= b from (not a < b).
    ('__lt__', '__ge__'): (True, '', None),

    # a >= b from (not a <= b) or (a == b)
    ('__le__', '__ge__'): (True, '||', False),
    # a < b, from (a <= b) and (a != b)
    ('__le__', '__lt__'): (False, '&&', True),
    # a > b from (not a <= b)
    ('__le__', '__gt__'): (True, '', None),

    # a < b from (not a > b) and (a != b)
    ('__gt__', '__lt__'): (True, '&&', True),
    # a >= b from (a > b) or (a == b)
    ('__gt__', '__ge__'): (False, '||', False),
    # a <= b from (not a > b)
    ('__gt__', '__le__'): (True, '', None),

    # Return a <= b from (not a >= b) or (a == b)
    ('__ge__', '__le__'): (True, '||', False),
    # a > b from (a >= b) and (a != b)
    ('__ge__', '__gt__'): (False, '&&', True),
    # a < b from (not a >= b)
    ('__ge__', '__lt__'): (True, '', None),
}

class SharedUtilityExporter:
    """
    Class responsible for generating code that imports and exports shared utility functions.

    Mark the positions where the functions should be called with `call_import_code()`/`call_export_code()`.
    The function calls and import/export functions are generated when `generate_exporting_functions()`
    is called. This approach is needed because the list of the shared functions is only known in the later
    stages of compilation.
    """
    def __init__(self, pos, mod_init_subfunction, scope):
        self.in_shared_utility_module = bool(scope.context.shared_c_file_path)
        self.using_shared_utility_module = bool(scope.context.shared_utility_qualified_name)
        self.pos = pos
        self.scope = scope
        self.import_code = mod_init_subfunction("Shared function import code")
        self.export_code = mod_init_subfunction("Shared function export code")

    def has_shared_exports(self, shared_func_definitions: Sequence[Code.SharedFunctionDecl]) -> bool:
        return bool(self.in_shared_utility_module and shared_func_definitions)

    def has_shared_imports(self, shared_func_definitions: Sequence[Code.SharedFunctionDecl]) -> bool:
        return bool(self.using_shared_utility_module and shared_func_definitions)

    def call_import_code(self, code):
        self.import_code.set_call_code(code)

    def call_export_code(self, code):
        self.export_code.set_call_code(code)

    def _generate_c_shared_function_export_code(self, code, shared_function_definitions: Sequence[Code.SharedFunctionDecl]):
        # We use the function cname also as exported name.
        exports = [
            (f"{shared_func_def.ret}({shared_func_def.params})", shared_func_def.name, shared_func_def.name)
            for shared_func_def in shared_function_definitions
        ]
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("FunctionExport", "ImportExport.c"))

        _generate_export_code(code, self.pos, exports, "__Pyx_ExportFunction", "void (*{name})(void)")

    def _generate_c_shared_function_import_code_for_module(self, code, function_definitions: Sequence[Code.SharedFunctionDecl]):
        # We use the function cname also as exported name.

        imports = [
            (f"{shared_func_def.ret}({shared_func_def.params})", shared_func_def.name, shared_func_def.name)
            for shared_func_def in function_definitions
        ]
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("FunctionImport", "ImportExport.c"))

        shared_utility_qualified_name = EncodedString(self.scope.context.shared_utility_qualified_name)
        import_func = f"__Pyx_ImportFunction_{Naming.cyversion}"
        _generate_import_code(
            code, self.pos, imports, shared_utility_qualified_name, import_func, "void (**{name})(void)")

    def _generate_exports(self, shared_utility_functions: Sequence[Code.SharedFunctionDecl]):
        if self.has_shared_exports(shared_utility_functions):
            with self.export_code as inner_code:
                self._generate_c_shared_function_export_code(
                    inner_code,
                    shared_utility_functions
                )

    def _generate_imports(self, shared_utility_functions: Sequence[Code.SharedFunctionDecl]):
        if self.has_shared_imports(shared_utility_functions):
            with self.import_code as inner_code:
                self._generate_c_shared_function_import_code_for_module(
                    inner_code,
                    shared_utility_functions
                )

    def generate_exporting_functions(self, code):
        shared_utility_functions = code.globalstate.shared_utility_functions
        code.enter_cfunc_scope(self.scope)
        self._generate_exports(shared_utility_functions)
        self._generate_imports(shared_utility_functions)
        code.exit_cfunc_scope()


class ModuleNode(Nodes.Node, Nodes.BlockNode):
    #  doc       string or None
    #  body      StatListNode
    #
    #  referenced_modules   [ModuleScope]
    #  full_module_name     string
    #
    #  scope                The module scope.
    #  compilation_source   A CompilationSource (see Main)
    #  directives           Top-level compiler directives

    child_attrs = ["body"]
    directives = None
    # internal - used in merging
    pxd_stats = None
    utility_code_stats = None

    @property
    def local_scope(self):
        # Make the module node (and its init function) look like a FuncDefNode.
        return self.scope

    def merge_in(self, tree, scope, stage):
        # Merges in the contents of another tree, and possibly scope. With the
        # current implementation below, this must be done right prior
        # to code generation.
        # Stage is one of "pxd" or "utility" to indicate pxd file or utility
        # code. This helps define the order.
        #
        # Note: This way of doing it seems strange -- I believe the
        # right concept is to split ModuleNode into a ModuleNode and a
        # CodeGenerator, and tell that CodeGenerator to generate code
        # from multiple sources.
        assert isinstance(self.body, Nodes.StatListNode)
        assert stage in ('pxd', 'utility')

        if self.pxd_stats is None:
            self.pxd_stats = Nodes.StatListNode(self.body.pos, stats=[])
            self.utility_code_stats = Nodes.StatListNode(self.body.pos, stats=[])
            self.body.stats.insert(0, self.pxd_stats)
            self.body.stats.insert(0, self.utility_code_stats)

        if scope.directives != self.scope.directives:
            # merged in nodes should keep their original compiler directives
            # (for example inline cdef functions)
            tree = Nodes.CompilerDirectivesNode(tree.pos, body=tree, directives=scope.directives)

        target_stats = self.pxd_stats if stage == "pxd" else self.utility_code_stats
        if isinstance(tree, Nodes.StatListNode):
            target_stats.stats.extend(tree.stats)
        else:
            target_stats.stats.append(tree)

        self.scope.utility_code_list.extend(scope.utility_code_list)

        for inc in scope.c_includes.values():
            self.scope.process_include(inc)

        def extend_if_not_in(L1, L2):
            for x in L2:
                if x not in L1:
                    L1.append(x)

        extend_if_not_in(self.scope.included_files, scope.included_files)

    def merge_scope(self, scope, internalise_c_class_entries=True):
        # Ensure that we don't generate import code for these entries!
        for entry in scope.c_class_entries:
            entry.type.module_name = self.full_module_name
            entry.type.scope.directives["internal"] = internalise_c_class_entries

        self.scope.merge_in(scope)

    def with_compiler_directives(self):
        # When merging a utility code module into the user code we need to preserve
        # the original compiler directives. This returns the body of the module node,
        # wrapped in its set of directives.
        body = Nodes.CompilerDirectivesNode(self.pos, directives=self.directives, body=self.body)
        return body

    def analyse_declarations(self, env):
        if has_np_pythran(env):
            Pythran.include_pythran_generic(env)
        if self.directives:
            env.old_style_globals = self.directives['old_style_globals']
        if not Options.docstrings:
            env.doc = self.doc = None
        elif Options.embed_pos_in_docstring:
            env.doc = EncodedString('File: %s (starting at line %s)' % Nodes.relative_position(self.pos))
            if self.doc is not None:
                env.doc = EncodedString(env.doc + '\n' + self.doc)
                env.doc.encoding = self.doc.encoding
        else:
            env.doc = self.doc
        env.directives = self.directives

        self.body.analyse_declarations(env)

        cy_pymutex_type = PyrexTypes.get_cy_pymutex_type()
        if env.find_shared_usages_of_type(cy_pymutex_type):
            # Be very suspicious of cython locks that are shared.
            # They have the potential to cause ABI issues.
            self.scope.use_utility_code(
                UtilityCode.load_cached(
                    "CythonPyMutexPublicCheck", "Synchronization.c"
                ))

    def prepare_utility_code(self):
        # prepare any utility code that must be created before code generation
        # specifically: CythonUtilityCode
        env = self.scope
        if env.has_import_star:
            self.create_import_star_conversion_utility_code(env)
        for name, entry in sorted(env.entries.items()):
            if (entry.create_wrapper and entry.scope is env
                    and entry.is_type and (entry.type.is_enum or entry.type.is_cpp_enum)):
                entry.type.create_type_wrapper(env)

    def process_implementation(self, options, result):
        env = self.scope
        env.return_type = PyrexTypes.c_void_type
        self.referenced_modules = []
        self.find_referenced_modules(env, self.referenced_modules, {})
        self.sort_cdef_classes(env)
        self.generate_c_code(env, options, result)
        self.generate_h_code(env, options, result)
        self.generate_api_code(env, options, result)

    def has_imported_c_functions(self):
        for module in self.referenced_modules:
            for entry in module.cfunc_entries:
                if entry.defined_in_pxd:
                    return 1
        return 0

    def assure_safe_target(self, path, allow_failed=False):
        # Check for a common gotcha for new users: naming your .pyx file after the .c file you want to wrap
        if not is_cython_generated_file(path, allow_failed=allow_failed, if_not_found=True):
            # Raising a fatal CompileError instead of calling error() to prevent castrating an existing file.
            raise CompileError(
                self.pos, 'The output file already exists and does not look like it was generated by Cython: "%s"' %
                          os.path.basename(path))

    def generate_h_code(self, env, options, result):
        def h_entries(entries, api=0, pxd=0):
            return [entry for entry in entries
                    if ((entry.visibility == 'public') or
                        (api and entry.api) or
                        (pxd and entry.defined_in_pxd))]
        h_types = h_entries(env.type_entries, api=1)
        h_vars = h_entries(env.var_entries)
        h_funcs = h_entries(env.cfunc_entries)
        h_extension_types = h_entries(env.c_class_entries)

        if h_types or h_vars or h_funcs or h_extension_types:
            result.h_file = replace_suffix_encoded(result.c_file, ".h")
            self.assure_safe_target(result.h_file)

            h_code_writer = Code.CCodeWriter()
            c_code_config = generate_c_code_config(env, options)
            globalstate = Code.GlobalState(h_code_writer, self, c_code_config)
            globalstate.initialize_main_h_code()  # in-case utility code is used in the header
            h_code_start = globalstate.parts['h_code']
            h_code_main = globalstate.parts['type_declarations']
            h_code_end = globalstate.parts['end']
            if options.generate_pxi:
                result.i_file = replace_suffix_encoded(result.c_file, ".pxi")
                i_code = Code.PyrexCodeWriter(result.i_file)
            else:
                i_code = None

            h_code_start.put_generated_by()
            h_guard = self.api_name(Naming.h_guard_prefix, env)
            h_code_start.put_h_guard(h_guard)
            h_code_start.putln("")
            h_code_start.putln('#include "Python.h"')
            self.generate_type_header_code(h_types, h_code_main)
            if options.capi_reexport_cincludes:
                self.generate_includes(env, [], h_code_main)
            h_code_main.putln("")
            api_guard = self.api_name(Naming.api_guard_prefix, env)
            h_code_main.putln("#ifndef %s" % api_guard)
            h_code_main.putln("")
            self.generate_extern_c_macro_definition(h_code_main, env.is_cpp())
            h_code_main.putln("")
            self.generate_dl_import_macro(h_code_main)
            if h_extension_types:
                h_code_main.putln("")
                for entry in h_extension_types:
                    self.generate_cclass_header_code(entry.type, h_code_main)
                    if i_code:
                        self.generate_cclass_include_code(entry.type, i_code)
                    globalstate.use_entry_utility_code(entry)
            if h_funcs:
                h_code_main.putln("")
                for entry in h_funcs:
                    self.generate_public_declaration(entry, h_code_main, i_code)
                    globalstate.use_entry_utility_code(entry)
            if h_vars:
                h_code_main.putln("")
                for entry in h_vars:
                    self.generate_public_declaration(entry, h_code_main, i_code)
                    globalstate.use_entry_utility_code(entry)
            h_code_main.putln("")
            h_code_main.putln("#endif /* !%s */" % api_guard)
            h_code_main.putln("")
            h_code_main.putln("/* WARNING: the interface of the module init function changed in CPython 3.5. */")
            h_code_main.putln("/* It now returns a PyModuleDef instance instead of a PyModule instance. */")
            h_code_main.putln("")
            py3_mod_func_name = self.mod_init_func_cname('PyInit', env)
            warning_string = EncodedString('Use PyImport_AppendInittab(%s, %s) instead of calling %s directly.' % (
                env.module_name.as_c_string_literal(), py3_mod_func_name, py3_mod_func_name))
            h_code_main.putln('/* WARNING: %s from Python 3.5 */' % warning_string.rstrip('.'))
            h_code_main.putln("PyMODINIT_FUNC %s(void);" % py3_mod_func_name)
            h_code_main.putln("")
            h_code_main.putln("#if PY_VERSION_HEX >= 0x03050000 "
                "&& (defined(__GNUC__) || defined(__clang__) || defined(_MSC_VER) "
                "|| (defined(__cplusplus) && __cplusplus >= 201402L))")
            h_code_main.putln("#if defined(__cplusplus) && __cplusplus >= 201402L")
            h_code_main.putln("[[deprecated(%s)]] inline" % warning_string.as_c_string_literal())
            h_code_main.putln("#elif defined(__GNUC__) || defined(__clang__)")
            h_code_main.putln('__attribute__ ((__deprecated__(%s), __unused__)) __inline__' % (
                warning_string.as_c_string_literal()))
            h_code_main.putln("#elif defined(_MSC_VER)")
            h_code_main.putln('__declspec(deprecated(%s)) __inline' % (
                warning_string.as_c_string_literal()))
            h_code_main.putln('#endif')
            h_code_main.putln("static PyObject* __PYX_WARN_IF_%s_INIT_CALLED(PyObject* res) {" % py3_mod_func_name)
            h_code_main.putln("return res;")
            h_code_main.putln("}")
            # Function call is converted to warning macro; uncalled (pointer) is not
            h_code_main.putln('#define %s() __PYX_WARN_IF_%s_INIT_CALLED(%s())' % (
                py3_mod_func_name, py3_mod_func_name, py3_mod_func_name))
            h_code_main.putln('#endif')

            h_code_end.putln("")
            h_code_end.putln("#endif /* !%s */" % h_guard)

            with open_new_file(result.h_file) as f:
                h_code_writer.copyto(f)

    def generate_public_declaration(self, entry, h_code, i_code):
        h_code.putln("%s %s;" % (
            Naming.extern_c_macro,
            entry.type.declaration_code(entry.cname)))
        if i_code:
            i_code.putln("cdef extern %s" % (
                entry.type.declaration_code(entry.cname, pyrex=1)))

    def api_name(self, prefix, env):
        api_name = self.punycode_module_name(prefix, env.qualified_name)
        return api_name.replace(".", "__")

    def generate_api_code(self, env, options, result):
        def api_entries(entries, pxd=0):
            return [entry for entry in entries
                    if entry.api or (pxd and entry.defined_in_pxd)]
        api_vars = api_entries(env.var_entries)
        api_funcs = api_entries(env.cfunc_entries)
        api_extension_types = api_entries(env.c_class_entries)

        if not (api_vars or api_funcs or api_extension_types):
            return

        result.api_file = replace_suffix_encoded(result.c_file, "_api.h")
        self.assure_safe_target(result.api_file)

        h_code = Code.CCodeWriter()
        c_code_config = generate_c_code_config(env, options)
        globalstate = Code.GlobalState(h_code, self, c_code_config)
        globalstate.initialize_main_h_code()  # in-case utility code is used in the header
        h_code.put_generated_by()
        api_guard = self.api_name(Naming.api_guard_prefix, env)
        h_code.put_h_guard(api_guard)
        # Work around https://bugs.python.org/issue4709
        h_code.putln('#ifdef __MINGW64__')
        h_code.putln('#define MS_WIN64')
        h_code.putln('#endif')

        def put_utility_code(name, src_file, include_requires=True):
            proto, impl = UtilityCode.load_as_string(name, src_file, include_requires=include_requires)
            if proto:
                h_code.put(proto)
            if impl:
                h_code.put(impl)

        h_code.putln('#include "Python.h"')
        if result.h_file:
            h_filename = os.path.basename(result.h_file)
            h_filename = as_encoded_filename(h_filename)
            h_code.putln('#include %s' % h_filename.as_c_string_literal())
        if api_extension_types:
            h_code.putln("")
            for entry in api_extension_types:
                type = entry.type
                h_code.putln("static PyTypeObject *%s = 0;" % type.typeptr_cname)
                h_code.putln("#define %s (*%s)" % (
                    type.typeobj_cname, type.typeptr_cname))
                h_code.globalstate.use_entry_utility_code(entry)
        if api_funcs:
            h_code.putln("")
            for entry in api_funcs:
                type = CPtrType(entry.type)
                cname = env.mangle(Naming.func_prefix_api, entry.name)
                h_code.putln("static %s = 0;" % type.declaration_code(cname))
                h_code.putln("#define %s %s" % (entry.name, cname))
                h_code.globalstate.use_entry_utility_code(entry)
        if api_vars:
            h_code.putln("")
            for entry in api_vars:
                type = CPtrType(entry.type)
                cname = env.mangle(Naming.varptr_prefix_api, entry.name)
                h_code.putln("static %s = 0;" % type.declaration_code(cname))
                h_code.putln("#define %s (*%s)" % (entry.name, cname))
                h_code.globalstate.use_entry_utility_code(entry)
        if api_vars:
            put_utility_code("VoidPtrImport", "ImportExport.c")
        if api_funcs:
            put_utility_code("FunctionImport", "ImportExport.c")
        if api_extension_types:
            put_utility_code("TypeImport", "ImportExport.c")
        h_code.putln("")
        h_code.putln("static int %s(void) {" % self.api_name("import", env))
        h_code.putln("PyObject *module = 0;")
        h_code.putln('module = PyImport_ImportModule(%s);' % env.qualified_name.as_c_string_literal())
        h_code.putln("if (!module) goto bad;")
        for entry in api_funcs:
            cname = env.mangle(Naming.func_prefix_api, entry.name)
            sig = entry.type.signature_string()
            h_code.putln(
                'if (__Pyx_ImportFunction_%s(module, %s, (void (**)(void))&%s, "%s") < 0) goto bad;'
                % (Naming.cyversion, entry.name.as_c_string_literal(), cname, sig))
        for entry in api_vars:
            cname = env.mangle(Naming.varptr_prefix_api, entry.name)
            sig = entry.type.empty_declaration_code()
            h_code.putln(
                'if (__Pyx_ImportVoidPtr_%s(module, %s, (void **)&%s, "%s") < 0) goto bad;'
                % (Naming.cyversion, entry.name.as_c_string_literal(), cname, sig))
        with ModuleImportGenerator(h_code, imported_modules={env.qualified_name: 'module'}) as import_generator:
            for entry in api_extension_types:
                self.generate_type_import_call(entry.type, h_code, import_generator, error_code="goto bad;", is_api=True)
        h_code.putln("Py_DECREF(module); module = 0;")
        h_code.putln("return 0;")
        h_code.putln("bad:")
        h_code.putln("Py_XDECREF(module);")
        h_code.putln("return -1;")
        h_code.putln("}")
        h_code.putln("")
        h_code.putln("#endif /* !%s */" % api_guard)

        f = open_new_file(result.api_file)
        try:
            h_code.copyto(f)
        finally:
            f.close()

    def generate_cclass_header_code(self, type, h_code):
        h_code.putln("%s %s %s;" % (
            Naming.extern_c_macro,
            PyrexTypes.public_decl("PyTypeObject", "DL_IMPORT"),
            type.typeobj_cname))

    def generate_cclass_include_code(self, type, i_code):
        i_code.putln("cdef extern class %s.%s:" % (
            type.module_name, type.name))
        i_code.indent()
        var_entries = type.scope.var_entries
        if var_entries:
            for entry in var_entries:
                i_code.putln("cdef %s" % (
                    entry.type.declaration_code(entry.cname, pyrex=1)))
        else:
            i_code.putln("pass")
        i_code.dedent()

    def generate_c_code(self, env, options, result):
        self.assure_safe_target(result.c_file, allow_failed=True)
        modules = self.referenced_modules

        if Options.annotate or options.annotate:
            show_entire_c_code = Options.annotate == "fullc" or options.annotate == "fullc"
            rootwriter = Annotate.AnnotationCCodeWriter(
                show_entire_c_code=show_entire_c_code,
                source_desc=self.compilation_source.source_desc,
            )
        else:
            rootwriter = Code.CCodeWriter()

        c_code_config = generate_c_code_config(env, options)

        globalstate = Code.GlobalState(
            rootwriter, self,
            code_config=c_code_config,
            common_utility_include_dir=options.common_utility_include_dir,
        )
        globalstate.initialize_main_c_code()
        h_code = globalstate['h_code']

        globalstate.module_pos = self.pos
        globalstate.directives = self.directives

        self.generate_module_preamble(env, options, modules, result.embedded_metadata, h_code)

        globalstate.use_utility_code(refnanny_utility_code)

        code = globalstate['before_global_var']
        code.putln('#define __Pyx_MODULE_NAME %s' %
                   self.full_module_name.as_c_string_literal())
        module_is_main = self.is_main_module_flag_cname()
        code.putln("extern int %s;" % module_is_main)
        code.putln("int %s = 0;" % module_is_main)
        code.putln("")
        code.putln("/* Implementation of %s */" % env.qualified_name.as_c_string_literal())

        code = globalstate['late_includes']
        self.generate_includes(env, modules, code, early=False)

        code = globalstate['module_code']

        self.generate_cached_builtins_decls(env, code)

        # generate normal variable and function definitions
        self.generate_lambda_definitions(env, code)
        self.generate_variable_definitions(env, code)
        self.body.generate_function_definitions(env, code)

        # generate extension types and methods
        code = globalstate['module_exttypes']
        self.generate_typeobj_definitions(env, code)
        self.generate_method_table(env, code)
        if env.has_import_star:
            self.generate_import_star(env, code)

        # initialise the macro to reduce the code size of one-time functionality
        globalstate['module_state'].put_code_here(
            UtilityCode.load("SmallCodeConfig", "ModuleSetupCode.c"))

        self.generate_module_state_start(env, globalstate['module_state'])
        self.generate_module_state_clear(env, globalstate['module_state_clear'])
        self.generate_module_state_traverse(env, globalstate['module_state_traverse'])

        shared_utility_exporter = SharedUtilityExporter(
            self.pos,
            self.mod_init_subfunction(self.pos, self.scope, globalstate['init_module']),
            self.scope
        )

        # init_globals is inserted before this
        self.generate_module_init_func(
            modules[:-1], shared_utility_exporter, env, globalstate['init_module']
        )
        self.generate_module_cleanup_func(env, globalstate['cleanup_module'])
        if Options.embed:
            self.generate_main_method(env, globalstate['main_method'])
        self.generate_filename_table(globalstate['filename_table'])

        self.generate_declarations_for_modules(env, modules, globalstate)
        h_code.write('\n')

        for utilcode in env.utility_code_list[:]:
            globalstate.use_utility_code(utilcode)

        shared_utility_exporter.generate_exporting_functions(code)

        globalstate.finalize_main_c_code()

        self.generate_module_state_end(env, modules, globalstate)

        f = open_new_file(result.c_file)
        try:
            rootwriter.copyto(f)
        finally:
            f.close()
        result.c_file_generated = 1
        if options.gdb_debug:
            self._serialize_lineno_map(env, rootwriter)
        if Options.annotate or options.annotate:
            self._generate_annotations(rootwriter, result, options)

    def _generate_annotations(self, rootwriter, result, options):
        self.annotate(rootwriter)

        coverage_xml_filename = Options.annotate_coverage_xml or options.annotate_coverage_xml
        if coverage_xml_filename and os.path.exists(coverage_xml_filename):
            import xml.etree.ElementTree as ET
            coverage_xml = ET.parse(coverage_xml_filename).getroot()
            for el in coverage_xml.iter():
                el.tail = None  # save some memory
        else:
            coverage_xml = None

        rootwriter.save_annotation(result.main_source_file, result.c_file, coverage_xml=coverage_xml)

        # if we included files, additionally generate one annotation file for each
        if not self.scope.included_files:
            return

        search_include_file = self.scope.context.search_include_directories
        target_dir = os.path.abspath(os.path.dirname(result.c_file))
        for included_file in self.scope.included_files:
            target_file = os.path.abspath(os.path.join(target_dir, included_file))
            target_file_dir = os.path.dirname(target_file)
            if not target_file_dir.startswith(target_dir):
                # any other directories may not be writable => avoid trying
                continue
            source_file = search_include_file(included_file, source_pos=self.pos, include=True)
            if not source_file:
                continue
            if target_file_dir != target_dir and not os.path.exists(target_file_dir):
                try:
                    os.makedirs(target_file_dir)
                except OSError as e:
                    import errno
                    if e.errno != errno.EEXIST:
                        raise
            rootwriter.save_annotation(source_file, target_file, coverage_xml=coverage_xml)

    def _serialize_lineno_map(self, env, ccodewriter):
        tb = env.context.gdb_debug_outputwriter
        markers = ccodewriter.buffer.allmarkers()

        d = defaultdict(list)
        for c_lineno, (src_desc, src_lineno) in enumerate(markers):
            if src_lineno > 0 and src_desc.filename is not None:
                d[src_desc, src_lineno].append(c_lineno + 1)

        tb.start('LineNumberMapping')
        for (src_desc, src_lineno), c_linenos in sorted(d.items()):
            assert src_desc.filename is not None
            tb.add_entry(
                'LineNumber',
                c_linenos=' '.join(map(str, c_linenos)),
                src_path=src_desc.filename,
                src_lineno=str(src_lineno),
            )
        tb.end('LineNumberMapping')
        tb.serialize()

    def find_referenced_modules(self, env, module_list, modules_seen):
        if env not in modules_seen:
            modules_seen[env] = 1
            for imported_module in env.cimported_modules:
                self.find_referenced_modules(imported_module, module_list, modules_seen)
            module_list.append(env)

    def sort_types_by_inheritance(self, type_dict, type_order, getkey):
        subclasses = defaultdict(list)  # maps type key to list of subclass keys
        for key in type_order:
            new_entry = type_dict[key]
            # collect all base classes to check for children
            base = new_entry.type.base_type
            while base:
                base_key = getkey(base)
                subclasses[base_key].append(key)
                base_entry = type_dict.get(base_key)
                if base_entry is None:
                    break
                base = base_entry.type.base_type

        # Simple topological sort using recursive DFS, based on
        # https://en.wikipedia.org/wiki/Topological_sorting#Depth-first_search
        seen = set()
        result = []
        def dfs(u):
            if u in seen:
                return
            seen.add(u)
            for v in subclasses[getkey(u.type)]:
                dfs(type_dict[v])
            result.append(u)

        for key in reversed(type_order):
            dfs(type_dict[key])

        result.reverse()
        return result

    def sort_type_hierarchy(self, module_list, env):
        # poor developer's OrderedDict
        vtab_dict, vtab_dict_order = {}, []
        vtabslot_dict, vtabslot_dict_order = {}, []

        for module in module_list:
            for entry in module.c_class_entries:
                if entry.used and not entry.in_cinclude:
                    type = entry.type
                    key = type.vtabstruct_cname
                    if not key:
                        continue
                    if key in vtab_dict:
                        # FIXME: this should *never* happen, but apparently it does
                        # for Cython generated utility code
                        from .UtilityCode import NonManglingModuleScope
                        assert isinstance(entry.scope, NonManglingModuleScope), str(entry.scope)
                        assert isinstance(vtab_dict[key].scope, NonManglingModuleScope), str(vtab_dict[key].scope)
                    else:
                        vtab_dict[key] = entry
                        vtab_dict_order.append(key)
            all_defined_here = module is env
            for entry in module.type_entries:
                if entry.used and (all_defined_here or entry.defined_in_pxd):
                    type = entry.type
                    if type.is_extension_type and not entry.in_cinclude:
                        type = entry.type
                        key = type.objstruct_cname
                        assert key not in vtabslot_dict, key
                        vtabslot_dict[key] = entry
                        vtabslot_dict_order.append(key)

        def vtabstruct_cname(entry_type):
            return entry_type.vtabstruct_cname
        vtab_list = self.sort_types_by_inheritance(
            vtab_dict, vtab_dict_order, vtabstruct_cname)

        def objstruct_cname(entry_type):
            return entry_type.objstruct_cname
        vtabslot_list = self.sort_types_by_inheritance(
            vtabslot_dict, vtabslot_dict_order, objstruct_cname)

        return (vtab_list, vtabslot_list)

    def sort_cdef_classes(self, env):
        key_func = operator.attrgetter('objstruct_cname')
        entry_dict, entry_order = {}, []
        for entry in env.c_class_entries:
            key = key_func(entry.type)
            assert key not in entry_dict, key
            entry_dict[key] = entry
            entry_order.append(key)
        env.c_class_entries[:] = self.sort_types_by_inheritance(
            entry_dict, entry_order, key_func)

    def generate_type_definitions(self, env, modules, vtab_list, vtabslot_list, code):
        # TODO: Why are these separated out?
        for entry in vtabslot_list:
            self.generate_objstruct_predeclaration(entry.type, code)
        vtabslot_entries = set(vtabslot_list)
        ctuple_names = set()
        for module in modules:
            definition = module is env
            type_entries = []
            for entry in module.type_entries:
                if entry.type.is_ctuple and entry.used:
                    if entry.name not in ctuple_names:
                        ctuple_names.add(entry.name)
                        type_entries.append(entry)
                elif definition or entry.defined_in_pxd:
                    type_entries.append(entry)
            type_entries = [t for t in type_entries if t not in vtabslot_entries]
            self.generate_type_header_code(type_entries, code)
        for entry in vtabslot_list:
            self.generate_objstruct_definition(entry.type, code)
            self.generate_typeobj_predeclaration(entry, code)
        for entry in vtab_list:
            self.generate_typeobj_predeclaration(entry, code)
            self.generate_exttype_vtable_struct(entry, code)
            self.generate_exttype_vtabptr_declaration(entry, code)
            self.generate_exttype_final_methods_declaration(entry, code)

    def generate_declarations_for_modules(self, env, modules, globalstate):
        typecode = globalstate['type_declarations']
        typecode.putln("")
        typecode.putln("/*--- Type declarations ---*/")
        # This is to work around the fact that array.h isn't part of the C-API,
        # but we need to declare it earlier than utility code.
        if 'cpython.array' in [m.qualified_name for m in modules]:
            typecode.putln('#ifndef _ARRAYARRAY_H')
            typecode.putln('struct arrayobject;')
            typecode.putln('typedef struct arrayobject arrayobject;')
            typecode.putln('#endif')
        vtab_list, vtabslot_list = self.sort_type_hierarchy(modules, env)
        self.generate_type_definitions(
            env, modules, vtab_list, vtabslot_list, typecode)
        modulecode = globalstate['module_declarations']
        for module in modules:
            defined_here = module is env
            modulecode.putln("")
            modulecode.putln("/* Module declarations from %s */" % module.qualified_name.as_c_string_literal())
            self.generate_c_class_declarations(module, modulecode, defined_here, globalstate)
            self.generate_cvariable_declarations(module, modulecode, defined_here)
            self.generate_cfunction_declarations(module, modulecode, defined_here)

    @staticmethod
    def _put_setup_code(code, name):
        code.put_code_here(UtilityCode.load(name, "ModuleSetupCode.c"))

    def generate_module_preamble(self, env, options, cimported_modules, metadata, code):
        code.put_generated_by()
        if metadata:
            code.putln("/* BEGIN: Cython Metadata")
            code.putln(json.dumps(metadata, indent=4, sort_keys=True))
            code.putln("END: Cython Metadata */")
            code.putln("")

        code.putln("#ifndef PY_SSIZE_T_CLEAN")
        code.putln("#define PY_SSIZE_T_CLEAN")
        code.putln("#endif /* PY_SSIZE_T_CLEAN */")
        self._put_setup_code(code, "InitLimitedAPI")

        for inc in sorted(env.c_includes.values(), key=IncludeCode.sortkey):
            if inc.location == inc.INITIAL:
                inc.write(code)
        code.putln("#ifndef Py_PYTHON_H")
        code.putln("    #error Python headers needed to compile C extensions, "
                   "please install development version of Python.")
        code.putln("#elif PY_VERSION_HEX < 0x03080000")
        code.putln("    #error Cython requires Python 3.8+.")
        code.putln("#else")
        code.globalstate["end"].putln("#endif /* Py_PYTHON_H */")

        from .. import __version__
        code.putln(f'#define __PYX_ABI_VERSION "{__version__.replace(".", "_")}"')
        code.putln('#define CYTHON_HEX_VERSION %s' % build_hex_version(__version__))
        code.putln("#define CYTHON_FUTURE_DIVISION %d" % (
            Future.division in env.context.future_directives))

        code.globalstate.use_utility_code(
            UtilityCode.load("CythonABIVersion", "ModuleSetupCode.c"))

        self._put_setup_code(code, "CModulePreamble")
        if env.context.options.cplus:
            self._put_setup_code(code, "CppInitCode")
        else:
            self._put_setup_code(code, "CInitCode")
        self._put_setup_code(code, "PythonCompatibility")
        self._put_setup_code(code, "MathInitCode")

        # Error handling and position macros.
        # Using "(void)cname" to prevent "unused" warnings.
        mark_errpos_code = (
            "#define __PYX_MARK_ERR_POS(f_index, lineno)  {"
            f" {Naming.filename_cname} = {Naming.filetable_cname}[f_index];"
            f" (void) {Naming.filename_cname};"
            f" {Naming.lineno_cname} = lineno;"
            f" (void) {Naming.lineno_cname};"
            "%s"  # for C line info
            f" (void) {Naming.clineno_cname}; "  # always suppress warnings
            "}"
        )
        cline_info = f" {Naming.clineno_cname} = {Naming.line_c_macro};"

        # Show the C code line in tracebacks or not? C macros take precedence over (deprecated) options.
        # 1) "CYTHON_CLINE_IN_TRACEBACK=0"  always disables C lines in tracebacks
        # 2) "CYTHON_CLINE_IN_TRACEBACK_RUNTIME=1" enables the feature + runtime configuration
        # 2a) "options.c_line_in_traceback=True"   changes the default to CYTHON_CLINE_IN_TRACEBACK_RUNTIME=1
        # 2b) "options.c_line_in_traceback=False"  changes the default to disable C lines
        # 4) "CYTHON_CLINE_IN_TRACEBACK=1"         enables C lines without runtime configuration
        # 5) if nothing is set, the default is to disable the feature

        default_cline_runtime = 0
        if options.c_line_in_traceback is not None:
            # explicitly set by user
            default_cline_runtime = int(options.c_line_in_traceback)

        code.putln("#ifndef CYTHON_CLINE_IN_TRACEBACK_RUNTIME")
        code.putln(f"#define CYTHON_CLINE_IN_TRACEBACK_RUNTIME {default_cline_runtime}")
        code.putln("#endif")

        code.putln("#ifndef CYTHON_CLINE_IN_TRACEBACK")
        code.putln("#define CYTHON_CLINE_IN_TRACEBACK CYTHON_CLINE_IN_TRACEBACK_RUNTIME")
        code.putln("#endif")

        code.putln("#if CYTHON_CLINE_IN_TRACEBACK")
        code.putln(mark_errpos_code % cline_info)
        code.putln("#else")
        code.putln(mark_errpos_code % "")
        code.putln("#endif")

        code.putln("#define __PYX_ERR(f_index, lineno, Ln_error) \\")
        code.putln("    { __PYX_MARK_ERR_POS(f_index, lineno) goto Ln_error; }")

        code.putln("")
        self.generate_extern_c_macro_definition(code, env.is_cpp())
        code.putln("")

        code.putln("#define %s" % self.api_name(Naming.h_guard_prefix, env))
        code.putln("#define %s" % self.api_name(Naming.api_guard_prefix, env))
        code.putln("/* Early includes */")
        self.generate_includes(env, cimported_modules, code, late=False)
        code.putln("")
        code.putln("#if defined(PYREX_WITHOUT_ASSERTIONS) && !defined(CYTHON_WITHOUT_ASSERTIONS)")
        code.putln("#define CYTHON_WITHOUT_ASSERTIONS")
        code.putln("#endif")
        code.putln("")

        if env.directives['ccomplex']:
            code.putln("")
            code.putln("#if !defined(CYTHON_CCOMPLEX)")
            code.putln("#define CYTHON_CCOMPLEX 1")
            code.putln("#endif")
            code.putln("")

        code.putln("#ifdef CYTHON_FREETHREADING_COMPATIBLE")
        code.putln("#if CYTHON_FREETHREADING_COMPATIBLE")
        code.putln("#define __Pyx_FREETHREADING_COMPATIBLE Py_MOD_GIL_NOT_USED")
        code.putln("#else")
        code.putln("#define __Pyx_FREETHREADING_COMPATIBLE Py_MOD_GIL_USED")
        code.putln("#endif")
        code.putln("#else")
        ft_compatible = "Py_MOD_GIL_NOT_USED" if env.directives["freethreading_compatible"] else "Py_MOD_GIL_USED"
        code.putln(f"#define __Pyx_FREETHREADING_COMPATIBLE {ft_compatible}")
        code.putln("#endif")

        c_string_type = env.directives['c_string_type']
        c_string_encoding = env.directives['c_string_encoding']
        if c_string_type not in ('bytes', 'bytearray') and not c_string_encoding:
            error(self.pos, "a default encoding must be provided if c_string_type is not a byte type")
        code.putln(f"#define __PYX_DEFAULT_STRING_ENCODING_IS_ASCII {int(c_string_encoding == 'ascii')}")
        code.putln(f"#define __PYX_DEFAULT_STRING_ENCODING_IS_UTF8 {int(c_string_encoding == 'utf8')}")
        if c_string_encoding not in ('ascii', 'utf8'):
            code.putln(f'#define __PYX_DEFAULT_STRING_ENCODING "{c_string_encoding}"')
        if c_string_type == 'bytearray':
            c_string_func_name = 'ByteArray'
        elif c_string_type == 'str':
            c_string_func_name = 'Unicode'
        else:
            c_string_func_name = c_string_type.title()
        code.putln(f'#define __Pyx_PyObject_FromString __Pyx_Py{c_string_func_name}_FromString')
        code.putln(f'#define __Pyx_PyObject_FromStringAndSize __Pyx_Py{c_string_func_name}_FromStringAndSize')
        code.put(UtilityCode.load_as_string("TypeConversions", "TypeConversion.c")[0])
        env.use_utility_code(UtilityCode.load_cached("FormatTypeName", "ObjectHandling.c"))

        # These utility functions are assumed to exist and used elsewhere.
        PyrexTypes.c_long_type.create_to_py_utility_code(env)
        PyrexTypes.c_long_type.create_from_py_utility_code(env)
        PyrexTypes.c_int_type.create_from_py_utility_code(env)

        code.put(Nodes.branch_prediction_macros)

        self._put_setup_code(code, "PretendToInitialize")
        code.putln('')
        code.putln('#if !CYTHON_USE_MODULE_STATE')
        code.putln('static PyObject *%s = NULL;' % env.module_cname)
        if Options.pre_import is not None:
            code.putln('static PyObject *%s;' % Naming.preimport_cname)
        code.putln('#endif')

        code.putln('static int %s;' % Naming.lineno_cname)
        code.putln('static int %s = 0;' % Naming.clineno_cname)
        code.putln('static const char * const %s = %s;' % (Naming.cfilenm_cname, Naming.file_c_macro))
        code.putln('static const char *%s;' % Naming.filename_cname)

        env.use_utility_code(UtilityCode.load_cached("FastTypeChecks", "ModuleSetupCode.c"))
        env.use_utility_code(UtilityCode.load("GetRuntimeVersion", "ModuleSetupCode.c"))
        env.use_utility_code(UtilityCode.load_cached("AddModuleRef", "ModuleSetupCode.c"))
        if has_np_pythran(env):
            env.use_utility_code(UtilityCode.load_cached("PythranConversion", "CppSupport.cpp"))

    def generate_extern_c_macro_definition(self, code, is_cpp):
        name = Naming.extern_c_macro
        code.putln("#ifdef CYTHON_EXTERN_C")
        # make sure that user overrides always take precedence
        code.putln('    #undef %s' % name)
        code.putln('    #define %s CYTHON_EXTERN_C' % name)
        code.putln("#elif defined(%s)" % name)
        code.putln("    #ifdef _MSC_VER")
        code.putln("    #pragma message (\"Please do not define the '%s' macro externally. Use 'CYTHON_EXTERN_C' instead.\")" % name)
        code.putln("    #else")
        code.putln("    #warning Please do not define the '%s' macro externally. Use 'CYTHON_EXTERN_C' instead." % name)
        code.putln("    #endif")
        code.putln("#else")
        if is_cpp:
            code.putln('    #define %s extern "C++"' % name)
        else:
            code.putln("  #ifdef __cplusplus")
            code.putln('    #define %s extern "C"' % name)
            code.putln("  #else")
            code.putln("    #define %s extern" % name)
            code.putln("  #endif")
        code.putln("#endif")

    def generate_dl_import_macro(self, code):
        code.putln("#ifndef DL_IMPORT")
        code.putln("  #define DL_IMPORT(_T) _T")
        code.putln("#endif")

    def generate_includes(self, env, cimported_modules, code, early=True, late=True):
        for inc in sorted(env.c_includes.values(), key=IncludeCode.sortkey):
            if inc.location == inc.EARLY:
                if early:
                    inc.write(code)
            elif inc.location == inc.LATE:
                if late:
                    inc.write(code)
        if early:
            code.putln_openmp("#include <omp.h>")

    def generate_filename_table(self, code):
        from os.path import isabs, basename
        code.putln("")
        code.putln("static const char* const %s[] = {" % Naming.filetable_cname)
        if code.globalstate.filename_list:
            for source_desc in code.globalstate.filename_list:
                file_path = source_desc.get_filenametable_entry()
                if isabs(file_path):
                    # never include absolute paths
                    file_path = source_desc.get_description()
                # Always use / as separator
                file_path = pathlib.Path(file_path).as_posix()
                escaped_filename = as_encoded_filename(file_path)
                code.putln('%s,' % escaped_filename.as_c_string_literal())
        else:
            # Some C compilers don't like an empty array
            code.putln("0")
        code.putln("};")

    def generate_type_predeclarations(self, env, code):
        pass

    def generate_type_header_code(self, type_entries, code):
        # Generate definitions of structs/unions/enums/typedefs/objstructs.
        #self.generate_gcc33_hack(env, code) # Is this still needed?
        # Forward declarations
        for entry in type_entries:
            if not entry.in_cinclude:
                #print "generate_type_header_code:", entry.name, repr(entry.type) ###
                type = entry.type
                if type.is_typedef:  # Must test this first!
                    pass
                elif type.is_struct_or_union or type.is_cpp_class:
                    self.generate_struct_union_predeclaration(entry, code)
                elif type.is_ctuple and not type.is_fused and entry.used:
                    self.generate_struct_union_predeclaration(entry.type.struct_entry, code)
                elif type.is_extension_type:
                    self.generate_objstruct_predeclaration(type, code)
        # Actual declarations
        for entry in type_entries:
            if not entry.in_cinclude:
                #print "generate_type_header_code:", entry.name, repr(entry.type) ###
                type = entry.type
                if type.is_typedef:  # Must test this first!
                    self.generate_typedef(entry, code)
                elif type.is_enum or type.is_cpp_enum:
                    self.generate_enum_definition(entry, code)
                elif type.is_struct_or_union:
                    self.generate_struct_union_definition(entry, code)
                elif type.is_ctuple and not type.is_fused and entry.used:
                    self.generate_struct_union_definition(entry.type.struct_entry, code)
                elif type.is_cpp_class:
                    self.generate_cpp_class_definition(entry, code)
                elif type.is_extension_type:
                    self.generate_objstruct_definition(type, code)
                if getattr(type, "scope", None):
                    for var_entry in type.scope.var_entries:
                        code.globalstate.use_entry_utility_code(var_entry)

    def generate_gcc33_hack(self, env, code):
        # Workaround for spurious warning generation in gcc 3.3
        code.putln("")
        for entry in env.c_class_entries:
            type = entry.type
            if not type.typedef_flag:
                name = type.objstruct_cname
                if name.startswith("__pyx_"):
                    tail = name[6:]
                else:
                    tail = name
                code.putln("typedef struct %s __pyx_gcc33_%s;" % (
                    name, tail))

    def generate_typedef(self, entry, code):
        base_type = entry.type.typedef_base_type
        enclosing_scope = entry.scope
        if base_type.is_numeric and not enclosing_scope.is_cpp_class_scope:
            try:
                writer = code.globalstate['numeric_typedefs']
            except KeyError:
                writer = code
        else:
            writer = code
        writer.mark_pos(entry.pos)
        writer.putln("typedef %s;" % base_type.declaration_code(entry.cname))

    def sue_predeclaration(self, type, kind, name):
        if type.typedef_flag:
            return "%s %s;\ntypedef %s %s %s;" % (
                kind, name,
                kind, name, name)
        else:
            return "%s %s;" % (kind, name)

    def generate_struct_union_predeclaration(self, entry, code):
        type = entry.type
        if type.is_cpp_class and type.templates:
            code.putln("template <typename %s>" % ", typename ".join(
                [T.empty_declaration_code() for T in type.templates]))
        code.putln(self.sue_predeclaration(type, type.kind, type.cname))

    def sue_header_footer(self, type, kind, name):
        header = "%s %s {" % (kind, name)
        footer = "};"
        return header, footer

    def generate_struct_union_definition(self, entry, code):
        code.mark_pos(entry.pos)
        type = entry.type
        scope = type.scope
        if scope:
            kind = type.kind
            packed = type.is_struct and type.packed
            if packed:
                kind = "%s %s" % (type.kind, "__Pyx_PACKED")
                code.globalstate.use_utility_code(packed_struct_utility_code)
            header, footer = \
                self.sue_header_footer(type, kind, type.cname)
            if packed:
                code.putln("#if defined(__SUNPRO_C)")
                code.putln("  #pragma pack(1)")
                code.putln("#elif !defined(__GNUC__)")
                code.putln("  #pragma pack(push, 1)")
                code.putln("#endif")
            code.putln(header)
            var_entries = scope.var_entries
            for attr in var_entries:
                code.putln(
                    "%s;" % attr.type.declaration_code(attr.cname))
            code.putln(footer)
            if packed:
                code.putln("#if defined(__SUNPRO_C)")
                code.putln("  #pragma pack()")
                code.putln("#elif !defined(__GNUC__)")
                code.putln("  #pragma pack(pop)")
                code.putln("#endif")

    def generate_cpp_constructor_code(self, arg_decls, arg_names, is_implementing, py_attrs, constructor, type, code):
        if is_implementing:
            code.putln("%s(%s) {" % (type.cname, ", ".join(arg_decls)))
            needs_gil = py_attrs or (constructor and not constructor.type.nogil)
            if needs_gil:
                code.put_ensure_gil()
            if py_attrs:
                for attr in py_attrs:
                    code.put_init_var_to_py_none(attr, nanny=False)
            if constructor:
                code.putln("%s(%s);" % (constructor.cname, ", ".join(arg_names)))
            if needs_gil:
                code.put_release_ensured_gil()
            code.putln("}")
        else:
            code.putln("%s(%s);" % (type.cname, ", ".join(arg_decls)))

    def generate_cpp_class_definition(self, entry, code):
        code.mark_pos(entry.pos)
        type = entry.type
        scope = type.scope
        if scope:
            if type.templates:
                code.putln("template <class %s>" % ", class ".join(
                    [T.empty_declaration_code() for T in type.templates]))
            # Just let everything be public.
            code.put("struct %s" % type.cname)
            if type.base_classes:
                base_class_decl = ", public ".join(
                    [base_class.empty_declaration_code() for base_class in type.base_classes])
                code.put(" : public %s" % base_class_decl)
            code.putln(" {")
            self.generate_type_header_code(scope.type_entries, code)
            py_attrs = [e for e in scope.entries.values()
                        if e.type.is_pyobject and not e.is_inherited]
            has_virtual_methods = False
            constructor = None
            destructor = None
            for attr in scope.var_entries:
                if attr.type.is_cfunction and attr.type.is_static_method:
                    code.put("static ")
                elif attr.name == "<init>":
                    constructor = scope.lookup_here("<init>")
                elif attr.name == "<del>":
                    destructor = attr
                elif attr.type.is_cfunction:
                    code.put("virtual ")
                    has_virtual_methods = True
                code.putln("%s;" % attr.type.declaration_code(attr.cname))
            is_implementing = 'init_module' in code.globalstate.parts

            if constructor or py_attrs:
                if constructor:
                    for constructor_alternative in constructor.all_alternatives():
                        arg_decls = []
                        arg_names = []
                        for arg in constructor_alternative.type.original_args[
                                :len(constructor_alternative.type.args)-constructor_alternative.type.optional_arg_count]:
                            arg_decls.append(arg.declaration_code())
                            arg_names.append(arg.cname)
                        if constructor_alternative.type.optional_arg_count:
                            arg_decls.append(constructor_alternative.type.op_arg_struct.declaration_code(Naming.optional_args_cname))
                            arg_names.append(Naming.optional_args_cname)
                        if not arg_decls:
                            default_constructor = True
                            arg_decls = []
                        self.generate_cpp_constructor_code(arg_decls, arg_names, is_implementing, py_attrs, constructor_alternative, type, code)
                else:
                    arg_decls = []
                    arg_names = []
                    self.generate_cpp_constructor_code(arg_decls, arg_names, is_implementing, py_attrs, constructor, type, code)

            if destructor or py_attrs or has_virtual_methods:
                if has_virtual_methods:
                    code.put("virtual ")
                if is_implementing:
                    code.putln("~%s() {" % type.cname)
                    if py_attrs:
                        code.put_ensure_gil()
                    if destructor:
                        code.putln("%s();" % destructor.cname)
                    if py_attrs:
                        for attr in py_attrs:
                            code.put_var_xdecref(attr, nanny=False)
                        code.put_release_ensured_gil()
                    code.putln("}")
                else:
                    code.putln("~%s();" % type.cname)
            if py_attrs:
                # Also need copy constructor and assignment operators.
                if is_implementing:
                    code.putln("%s(const %s& __Pyx_other) {" % (type.cname, type.cname))
                    code.put_ensure_gil()
                    for attr in scope.var_entries:
                        if not attr.type.is_cfunction:
                            code.putln("%s = __Pyx_other.%s;" % (attr.cname, attr.cname))
                            code.put_var_incref(attr, nanny=False)
                    code.put_release_ensured_gil()
                    code.putln("}")
                    code.putln("%s& operator=(const %s& __Pyx_other) {" % (type.cname, type.cname))
                    code.putln("if (this != &__Pyx_other) {")
                    code.put_ensure_gil()
                    for attr in scope.var_entries:
                        if not attr.type.is_cfunction:
                            code.put_var_xdecref(attr, nanny=False)
                            code.putln("%s = __Pyx_other.%s;" % (attr.cname, attr.cname))
                            code.put_var_incref(attr, nanny=False)
                    code.put_release_ensured_gil()
                    code.putln("}")
                    code.putln("return *this;")
                    code.putln("}")
                else:
                    code.putln("%s(const %s& __Pyx_other);" % (type.cname, type.cname))
                    code.putln("%s& operator=(const %s& __Pyx_other);" % (type.cname, type.cname))
            code.putln("};")

    def generate_enum_definition(self, entry, code):
        code.mark_pos(entry.pos)
        type = entry.type
        name = entry.cname or entry.name or ""

        kind = "enum class" if entry.type.is_cpp_enum else "enum"
        header, footer = self.sue_header_footer(type, kind, name)
        code.putln(header)
        enum_values = entry.enum_values
        if not enum_values:
            error(entry.pos, "Empty enum definition not allowed outside a 'cdef extern from' block")
        else:
            last_entry = enum_values[-1]
            # this does not really generate code, just builds the result value
            for value_entry in enum_values:
                if value_entry.value_node is not None:
                    value_entry.value_node.generate_evaluation_code(code)

            for value_entry in enum_values:
                if value_entry.value_node is None:
                    value_code = value_entry.cname.split("::")[-1]
                else:
                    value_code = ("%s = %s" % (
                        value_entry.cname.split("::")[-1],
                        value_entry.value_node.result()))
                if value_entry is not last_entry:
                    value_code += ","
                code.putln(value_code)
        code.putln(footer)

        if entry.type.is_enum:
            if entry.type.typedef_flag:
                # Not pre-declared.
                code.putln("typedef enum %s %s;" % (name, name))

    def generate_typeobj_predeclaration(self, entry, code):
        code.putln("")
        name = entry.type.typeobj_cname
        if name:
            if entry.visibility == 'extern' and not entry.in_cinclude:
                code.putln("%s %s %s;" % (
                    Naming.extern_c_macro,
                    PyrexTypes.public_decl("PyTypeObject", "DL_IMPORT"),
                    name))
            elif entry.visibility == 'public':
                code.putln("%s %s %s;" % (
                    Naming.extern_c_macro,
                    PyrexTypes.public_decl("PyTypeObject", "DL_EXPORT"),
                    name))
            # ??? Do we really need the rest of this? ???
            #else:
            #    code.putln("static PyTypeObject %s;" % name)

    def generate_exttype_vtable_struct(self, entry, code):
        if not entry.used:
            return

        code.mark_pos(entry.pos)
        # Generate struct declaration for an extension type's vtable.
        type = entry.type
        scope = type.scope

        self.specialize_fused_types(scope)

        if type.vtabstruct_cname:
            code.putln("")
            code.putln("struct %s {" % type.vtabstruct_cname)
            if type.base_type and type.base_type.vtabstruct_cname:
                code.putln("struct %s %s;" % (
                    type.base_type.vtabstruct_cname,
                    Naming.obj_base_cname))
            for method_entry in scope.cfunc_entries:
                if not method_entry.is_inherited:
                    code.putln("%s;" % method_entry.type.declaration_code("(*%s)" % method_entry.cname))
            code.putln("};")

    def generate_exttype_vtabptr_declaration(self, entry, code):
        if not entry.used:
            return

        code.mark_pos(entry.pos)
        # Generate declaration of pointer to an extension type's vtable.
        type = entry.type
        if type.vtabptr_cname:
            code.putln("static struct %s *%s;" % (
                type.vtabstruct_cname,
                type.vtabptr_cname))

    def generate_exttype_final_methods_declaration(self, entry, code):
        if not entry.used:
            return

        code.mark_pos(entry.pos)
        # Generate final methods prototypes
        for method_entry in entry.type.scope.cfunc_entries:
            if not method_entry.is_inherited and method_entry.final_func_cname:
                declaration = method_entry.type.declaration_code(
                    method_entry.final_func_cname)
                modifiers = code.build_function_modifiers(method_entry.func_modifiers)
                code.putln("static %s%s;" % (modifiers, declaration))

    def generate_objstruct_predeclaration(self, type, code):
        if not type.scope:
            return
        code.putln(self.sue_predeclaration(type, "struct", type.objstruct_cname))

    def generate_objstruct_definition(self, type, code):
        code.mark_pos(type.pos)
        # Generate object struct definition for an
        # extension type.
        if not type.scope:
            return  # Forward declared but never defined
        header, footer = \
            self.sue_header_footer(type, "struct", type.objstruct_cname)
        code.putln(header)
        base_type = type.base_type
        if base_type:
            basestruct_cname = base_type.objstruct_cname
            if basestruct_cname == "PyTypeObject":
                # User-defined subclasses of type are heap allocated.
                basestruct_cname = "PyHeapTypeObject"
            code.putln(
                "%s%s %s;" % (
                    ("struct ", "")[base_type.typedef_flag],
                    basestruct_cname,
                    Naming.obj_base_cname))
        else:
            code.putln(
                "PyObject_HEAD")
        if type.vtabslot_cname and not (type.base_type and type.base_type.vtabslot_cname):
            code.putln(
                "struct %s *%s;" % (
                    type.vtabstruct_cname,
                    type.vtabslot_cname))
        for attr in type.scope.var_entries:
            if attr.is_declared_generic:
                attr_type = py_object_type
            else:
                attr_type = attr.type
            if attr.is_cpp_optional:
                decl = attr_type.cpp_optional_declaration_code(attr.cname)
            else:
                decl = attr_type.declaration_code(attr.cname)
            code.globalstate.use_entry_utility_code(attr)
            code.putln("%s;" % decl)
        code.putln(footer)
        if type.objtypedef_cname is not None:
            # Only for exposing public typedef name.
            code.putln("typedef struct %s %s;" % (type.objstruct_cname, type.objtypedef_cname))

    def generate_c_class_declarations(self, env, code, definition, globalstate):
        module_state = globalstate['module_state']
        module_state_clear = globalstate['module_state_clear']
        module_state_traverse = globalstate['module_state_traverse']
        module_state_typeobj = module_state.insertion_point()
        for entry in env.c_class_entries:
            if definition or entry.defined_in_pxd:
                module_state.putln("PyTypeObject *%s;" % entry.type.typeptr_cname)
                module_state_clear.putln(
                    "Py_CLEAR(clear_module_state->%s);" %
                    entry.type.typeptr_cname)
                module_state_traverse.putln(
                    "Py_VISIT(traverse_module_state->%s);" %
                    entry.type.typeptr_cname)
                if entry.type.typeobj_cname is not None:
                    module_state_typeobj.putln("PyObject *%s;" % entry.type.typeobj_cname)
                    module_state_clear.putln(
                        "Py_CLEAR(clear_module_state->%s);" % (
                        entry.type.typeobj_cname))
                    module_state_traverse.putln(
                        "Py_VISIT(traverse_module_state->%s);" % (
                        entry.type.typeobj_cname))

    def generate_cvariable_declarations(self, env, code, definition):
        if env.is_cython_builtin:
            return
        for entry in env.var_entries:
            if (entry.in_cinclude or entry.in_closure or
                    (entry.visibility == 'private' and not (entry.defined_in_pxd or entry.used))):
                continue

            storage_class = None
            dll_linkage = None
            init = None

            if entry.visibility == 'extern':
                storage_class = Naming.extern_c_macro
                dll_linkage = "DL_IMPORT"
            elif entry.visibility == 'public':
                storage_class = Naming.extern_c_macro
                if definition:
                    dll_linkage = "DL_EXPORT"
                else:
                    dll_linkage = "DL_IMPORT"
            elif entry.visibility == 'private':
                storage_class = "static"
                dll_linkage = None
                if entry.init is not None:
                    init = entry.type.literal_code(entry.init)
            type = entry.type
            cname = entry.cname

            if entry.defined_in_pxd and not definition:
                storage_class = "static"
                dll_linkage = None
                type = CPtrType(type)
                cname = env.mangle(Naming.varptr_prefix, entry.name)
                init = 0

            if storage_class:
                code.put("%s " % storage_class)
            if entry.is_cpp_optional:
                code.put(type.cpp_optional_declaration_code(
                    cname, dll_linkage=dll_linkage))
            else:
                code.put(type.declaration_code(
                    cname, dll_linkage=dll_linkage))
            if init is not None:
                code.put_safe(" = %s" % init)
            code.putln(";")
            if entry.cname != cname:
                code.putln("#define %s (*%s)" % (entry.cname, cname))
            code.globalstate.use_entry_utility_code(entry)

    def generate_cfunction_declarations(self, env, code, definition):
        for entry in env.cfunc_entries:
            from_pyx = Options.cimport_from_pyx and not entry.visibility == 'extern'
            if (entry.used
                    or entry.visibility == 'public'
                    or entry.api
                    or from_pyx):
                generate_cfunction_declaration(entry, env, code, definition)

    def generate_variable_definitions(self, env, code):
        for entry in env.var_entries:
            if not entry.in_cinclude and entry.visibility == "public":
                code.put(entry.type.declaration_code(entry.cname))
                if entry.init is not None:
                    init = entry.type.literal_code(entry.init)
                    code.put_safe(" = %s" % init)
                code.putln(";")

    def generate_typeobj_definitions(self, env, code):
        full_module_name = env.qualified_name
        for entry in env.c_class_entries:
            #print "generate_typeobj_definitions:", entry.name
            #print "...visibility =", entry.visibility
            if entry.visibility != 'extern':
                type = entry.type
                scope = type.scope
                if scope:  # could be None if there was an error
                    self.generate_exttype_vtable(scope, code)
                    self.generate_new_function(scope, code, entry)
                    self.generate_del_function(scope, code)
                    self.generate_dealloc_function(scope, code)

                    if scope.needs_gc():
                        self.generate_traverse_function(scope, code, entry)
                        if scope.needs_tp_clear():
                            self.generate_clear_function(scope, code, entry)
                    if scope.defines_any_special(["__getitem__"]):
                        self.generate_getitem_int_function(scope, code)
                    if scope.defines_any_special(["__setitem__", "__delitem__"]):
                        self.generate_ass_subscript_function(scope, code)
                    if scope.defines_any_special(["__getslice__", "__setslice__", "__delslice__"]):
                        warning(self.pos,
                                "__getslice__, __setslice__, and __delslice__ are not supported by Python 3, "
                                "use __getitem__, __setitem__, and __delitem__ instead", 1)
                        code.putln("#error __getslice__, __setslice__, and __delslice__ not supported in Python 3.")
                    if scope.defines_any_special(["__setslice__", "__delslice__"]):
                        self.generate_ass_slice_function(scope, code)
                    if scope.defines_any_special(["__getattr__", "__getattribute__"]):
                        self.generate_getattro_function(scope, code)
                    if scope.defines_any_special(["__setattr__", "__delattr__"]):
                        self.generate_setattro_function(scope, code)
                    if scope.defines_any_special(["__get__"]):
                        self.generate_descr_get_function(scope, code)
                    if scope.defines_any_special(["__set__", "__delete__"]):
                        self.generate_descr_set_function(scope, code)
                    if not (scope.is_closure_class_scope or scope.is_defaults_class_scope) and scope.defines_any(["__dict__"]):
                        self.generate_dict_getter_function(scope, code)

                    if scope.defines_any_special(TypeSlots.richcmp_special_methods):
                        self.generate_richcmp_function(scope, code)
                    elif 'total_ordering' in scope.directives:
                        # Warn if this is used when it can't have any effect.
                        warning(scope.parent_type.pos,
                                "total_ordering directive used, but no comparison and equality methods defined")

                    for slot in TypeSlots.get_slot_table(code.globalstate.directives).PyNumberMethods:
                        if slot.is_binop and scope.defines_any_special(slot.user_methods):
                            self.generate_binop_function(scope, slot, code, entry.pos)

                    self.generate_property_accessors(scope, code)
                    self.generate_method_table(scope, code)
                    self.generate_getset_table(scope, code)
                    code.putln("#if CYTHON_USE_TYPE_SPECS")
                    self.generate_typeobj_spec(entry, code)
                    code.putln("#else")
                    self.generate_typeobj_definition(full_module_name, entry, code)
                    code.putln("#endif")

    def generate_exttype_vtable(self, scope, code):
        # Generate the definition of an extension type's vtable.
        type = scope.parent_type
        if type.vtable_cname:
            code.putln("static struct %s %s;" % (
                type.vtabstruct_cname,
                type.vtable_cname))

    def generate_self_cast(self, scope, code):
        type = scope.parent_type
        code.putln(
            "%s = (%s)o;" % (
                type.declaration_code("p"),
                type.empty_declaration_code()))

    @staticmethod
    def generate_freelist_condition(code, size_check, type_cname, type):
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("CheckTypeForFreelists", "ExtensionTypes.c"))
        if type.is_final_type:
            freelist_check = '__PYX_CHECK_FINAL_TYPE_FOR_FREELISTS'
        else:
            freelist_check = '__PYX_CHECK_TYPE_FOR_FREELISTS'
        obj_struct = type.declaration_code("", deref=True)
        typeptr_cname = code.name_in_slot_module_state(type.typeptr_cname)
        code.putln(
            f"if (likely((int)({size_check}) & {freelist_check}({type_cname}, {typeptr_cname}, sizeof({obj_struct}))))")

    def generate_new_function(self, scope, code, cclass_entry):
        tp_slot = TypeSlots.ConstructorSlot("tp_new", "__cinit__")
        slot_func = scope.mangle_internal("tp_new")
        if tp_slot.slot_code(scope) != slot_func:
            return  # never used

        type = scope.parent_type
        base_type = type.base_type

        have_entries, (py_attrs, py_buffers, memoryview_slices) = \
                        scope.get_refcounted_entries()
        is_final_type = scope.parent_type.is_final_type
        if scope.is_internal:
            # internal classes (should) never need None inits, normal zeroing will do
            py_attrs = []
        explicitly_constructable_attrs = [
            entry for entry in scope.var_entries
            if entry.type.needs_explicit_construction(scope)
        ]

        cinit_func_entry = scope.lookup_here("__cinit__")
        if cinit_func_entry and not cinit_func_entry.is_special:
            cinit_func_entry = None

        if base_type or (cinit_func_entry and not cinit_func_entry.trivial_signature):
            unused_marker = ''
        else:
            unused_marker = 'CYTHON_UNUSED '

        if base_type:
            freelist_size = 0  # not currently supported
        else:
            freelist_size = scope.directives.get('freelist', 0)
        freelist_name = scope.mangle_internal(Naming.freelist_name)
        freecount_name = scope.mangle_internal(Naming.freecount_name)

        if freelist_size:
            module_state = code.globalstate['module_state_contents']
            module_state.putln("")
            module_state.putln("#if CYTHON_USE_FREELISTS")
            module_state.putln("%s[%d];" % (
                scope.parent_type.declaration_code(freelist_name),
                freelist_size))
            module_state.putln("int %s;" % freecount_name)
            module_state.putln("#endif")

        code.start_slotfunc(
            scope, PyrexTypes.py_objptr_type, "tp_new",
            f"PyTypeObject *t, {unused_marker}PyObject *a, {unused_marker}PyObject *k", needs_prototype=True)

        need_self_cast = (type.vtabslot_cname or
                          (py_buffers or memoryview_slices or py_attrs) or
                          explicitly_constructable_attrs)
        if need_self_cast:
            code.putln("%s;" % scope.parent_type.declaration_code("p"))
        if base_type:
            tp_new = TypeSlots.get_base_slot_function(scope, tp_slot)
            base_type_typeptr_cname = base_type.typeptr_cname
            if not base_type.is_builtin_type:
                base_type_typeptr_cname = code.name_in_slot_module_state(base_type_typeptr_cname)
            if tp_new is None:
                tp_new = f"__Pyx_PyType_GetSlot({base_type_typeptr_cname}, tp_new, newfunc)"
            code.putln("PyObject *o = %s(t, a, k);" % tp_new)
        else:
            code.putln("PyObject *o;")
            if freelist_size:
                code.globalstate.use_utility_code(
                    UtilityCode.load_cached("IncludeStringH", "StringTools.c"))
                code.putln("#if CYTHON_USE_FREELISTS")
                freecount_name = code.name_in_slot_module_state(freecount_name)
                freelist_name = code.name_in_slot_module_state(freelist_name)
                self.generate_freelist_condition(code, f"{freecount_name} > 0", "t", type)
                code.putln("{")
                code.putln("o = (PyObject*)%s[--%s];" % (
                    freelist_name,
                    freecount_name))
                obj_struct = type.declaration_code("", deref=True)
                code.putln("#if CYTHON_USE_TYPE_SPECS")
                # We still hold a reference to the type object held by the previous
                # user of the freelist object - release it.
                code.putln("Py_DECREF(Py_TYPE(o));")
                code.putln("#endif")
                code.putln("memset(o, 0, sizeof(%s));" % obj_struct)
                code.putln("#if CYTHON_COMPILING_IN_LIMITED_API")
                # Although PyObject_INIT should be part of the Limited API, it causes
                # link errors on some combinations of Python versions and OSs.
                code.putln("(void) PyObject_Init(o, t);")
                code.putln("#else")
                code.putln("(void) PyObject_INIT(o, t);")
                code.putln("#endif")
                if scope.needs_gc():
                    code.putln("PyObject_GC_Track(o);")
                code.putln("} else")
                code.putln("#endif")
                code.putln("{")
            code.globalstate.use_utility_code(
                UtilityCode.load_cached("AllocateExtensionType", "ExtensionTypes.c")
            )
            code.putln(f"o = __Pyx_AllocateExtensionType(t, {is_final_type:d});")
        code.putln("if (unlikely(!o)) return 0;")
        if freelist_size and not base_type:
            code.putln('}')
        if need_self_cast:
            code.putln("p = %s;" % type.cast_code("o"))
        #if need_self_cast:
        #    self.generate_self_cast(scope, code)

        # from this point on, ensure DECREF(o) on failure
        needs_error_cleanup = False

        if type.vtabslot_cname:
            vtab_base_type = type
            while vtab_base_type.base_type and vtab_base_type.base_type.vtabstruct_cname:
                vtab_base_type = vtab_base_type.base_type
            if vtab_base_type is not type:
                struct_type_cast = "(struct %s*)" % vtab_base_type.vtabstruct_cname
            else:
                struct_type_cast = ""
            code.putln("p->%s = %s%s;" % (
                type.vtabslot_cname,
                struct_type_cast, type.vtabptr_cname))

        for entry in explicitly_constructable_attrs:
            entry.type.generate_explicit_construction(
                code, entry, extra_access_code="p->")

        for entry in py_attrs:
            if entry.name == "__dict__":
                needs_error_cleanup = True
                code.put("p->%s = PyDict_New(); if (unlikely(!p->%s)) goto bad;" % (
                    entry.cname, entry.cname))
            else:
                code.put_init_var_to_py_none(entry, "p->%s", nanny=False)

        for entry in memoryview_slices:
            code.putln("p->%s.data = NULL;" % entry.cname)
            code.putln("p->%s.memview = NULL;" % entry.cname)

        for entry in py_buffers:
            code.putln("p->%s.obj = NULL;" % entry.cname)

        if cclass_entry.cname == '__pyx_memoryviewslice':
            code.putln("p->from_slice.memview = NULL;")

        if cinit_func_entry:
            if cinit_func_entry.trivial_signature:
                cinit_args = f"o, {Naming.modulestateglobal_cname}->{Naming.empty_tuple}, NULL"
            else:
                cinit_args = "o, a, k"
            needs_error_cleanup = True
            code.putln("if (unlikely(%s(%s) < 0)) goto bad;" % (
                cinit_func_entry.func_cname, cinit_args))

        code.putln(
            "return o;")
        if needs_error_cleanup:
            code.putln("bad:")
            code.put_decref_clear("o", py_object_type, nanny=False)
            code.putln("return NULL;")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_del_function(self, scope, code):
        tp_slot = TypeSlots.get_slot_by_name("tp_finalize", scope.directives)
        slot_func_cname = scope.mangle_internal("tp_finalize")
        if tp_slot.slot_code(scope) != slot_func_cname:
            return  # never used

        entry = scope.lookup_here("__del__")
        if entry is None or not entry.is_special:
            return  # nothing to wrap
        code.putln("")

        if tp_slot.used_ifdef:
            code.putln("#if %s" % tp_slot.used_ifdef)

        code.start_slotfunc(scope, PyrexTypes.c_void_type, "tp_finalize", "PyObject *o", needs_funcstate=False)
        code.putln("PyObject *etype, *eval, *etb;")
        code.putln("PyErr_Fetch(&etype, &eval, &etb);")
        code.putln("%s(o);" % entry.func_cname)
        code.putln("PyErr_Restore(etype, eval, etb);")
        code.putln("}")
        code.exit_cfunc_scope()

        if tp_slot.used_ifdef:
            code.putln("#endif")

    def generate_dealloc_function(self, scope, code):
        tp_slot = TypeSlots.ConstructorSlot("tp_dealloc", '__dealloc__')
        slot_func = scope.mangle_internal("tp_dealloc")
        base_type = scope.parent_type.base_type
        if tp_slot.slot_code(scope) != slot_func:
            return  # never used

        slot_func_cname = scope.mangle_internal("tp_dealloc")
        code.start_slotfunc(scope, PyrexTypes.c_void_type, "tp_dealloc", "PyObject *o")

        is_final_type = scope.parent_type.is_final_type
        needs_gc = scope.needs_gc()
        needs_trashcan = scope.needs_trashcan()

        weakref_slot = scope.lookup_here("__weakref__") if not (scope.is_closure_class_scope or scope.is_defaults_class_scope) else None
        if weakref_slot not in scope.var_entries:
            weakref_slot = None

        dict_slot = scope.lookup_here("__dict__") if not (scope.is_closure_class_scope or scope.is_defaults_class_scope) else None
        if dict_slot not in scope.var_entries:
            dict_slot = None

        _, (py_attrs, _, memoryview_slices) = scope.get_refcounted_entries()
        explicitly_destructable_attrs = [
            entry for entry in scope.var_entries
            if entry.type.needs_explicit_destruction(scope)
        ]

        if py_attrs or explicitly_destructable_attrs or memoryview_slices or weakref_slot or dict_slot:
            self.generate_self_cast(scope, code)

        if not is_final_type or scope.may_have_finalize():
            # in Py3.4+, call tp_finalize() as early as possible
            code.putln("#if CYTHON_USE_TP_FINALIZE")
            if needs_gc:
                finalised_check = '!__Pyx_PyObject_GC_IsFinalized(o)'
            else:
                finalised_check = (
                    '(!PyType_IS_GC(Py_TYPE(o)) || !__Pyx_PyObject_GC_IsFinalized(o))')
            code.putln(
                "if (unlikely(__Pyx_PyObject_GetSlot(o, tp_finalize, destructor)) && %s) {" % finalised_check)

            code.putln("if (__Pyx_PyObject_GetSlot(o, tp_dealloc, destructor) == %s) {" % slot_func_cname)
            # if instance was resurrected by finaliser, return
            code.putln("if (PyObject_CallFinalizerFromDealloc(o)) return;")
            code.putln("}")
            code.putln("}")
            code.putln("#endif")

        if needs_gc:
            # We must mark this object as (gc) untracked while tearing
            # it down, lest the garbage collection is invoked while
            # running this destructor.
            code.putln("PyObject_GC_UnTrack(o);")

        if needs_trashcan:
            code.globalstate.use_utility_code(
                UtilityCode.load_cached("PyTrashcan", "ExtensionTypes.c"))
            code.putln("__Pyx_TRASHCAN_BEGIN(o, %s)" % slot_func_cname)

        if weakref_slot:
            # We must clean the weakreferences before calling the user's __dealloc__
            # because if the __dealloc__ releases the GIL, a weakref can be
            # dereferenced accessing the object in an inconsistent state or
            # resurrecting it.
            code.putln("if (p->__weakref__) PyObject_ClearWeakRefs(o);")

        # call the user's __dealloc__
        self.generate_usr_dealloc_call(scope, code)

        if dict_slot:
            code.putln("if (p->__dict__) PyDict_Clear(p->__dict__);")

        for entry in explicitly_destructable_attrs:
            entry.type.generate_explicit_destruction(code, entry, extra_access_code="p->")

        for entry in (py_attrs + memoryview_slices):
            code.put_xdecref_clear("p->%s" % entry.cname, entry.type, nanny=False,
                                   clear_before_decref=True, have_gil=True)

        if base_type:
            base_cname = base_type.typeptr_cname
            if not base_type.is_builtin_type:
                base_cname = code.name_in_slot_module_state(base_cname)
            tp_dealloc = TypeSlots.get_base_slot_function(scope, tp_slot)
            if tp_dealloc is not None:
                if needs_gc and base_type.scope and base_type.scope.needs_gc():
                    # We know that the base class uses GC, so probably expects it to be tracked.
                    # Undo the untracking above.
                    code.putln("PyObject_GC_Track(o);")
                code.putln("%s(o);" % tp_dealloc)
            elif base_type.is_builtin_type:
                if needs_gc and base_type.scope and base_type.scope.needs_gc():
                    # We know that the base class uses GC, so probably expects it to be tracked.
                    # Undo the untracking above.
                    code.putln("PyObject_GC_Track(o);")
                code.putln("__Pyx_PyType_GetSlot(%s, tp_dealloc, destructor)(o);" % base_cname)
            else:
                if needs_gc:
                    # We don't know if the base class uses GC or not, so must find out at runtime
                    # whether we should undo the untracking above or not.
                    code.putln("if (PyType_IS_GC(%s)) PyObject_GC_Track(o);" % base_cname)
                # This is an externally defined type.  Calling through the
                # cimported base type pointer directly interacts badly with
                # the module cleanup, which may already have cleared it.
                # In that case, fall back to traversing the type hierarchy.
                # If we're using the module state then always go through the
                # type hierarchy, because our access to the module state may
                # have been lost (at least for the limited API version of
                # using module state).
                code.putln("#if !CYTHON_USE_MODULE_STATE")
                code.putln("if (likely(%s)) __Pyx_PyType_GetSlot(%s, tp_dealloc, destructor)(o); else" % (
                    base_cname, base_cname))
                code.putln("#endif")
                code.putln("__Pyx_call_next_tp_dealloc(o, %s);" % slot_func_cname)
                code.globalstate.use_utility_code(
                    UtilityCode.load_cached("CallNextTpDealloc", "ExtensionTypes.c"))
        else:
            freelist_size = scope.directives.get('freelist', 0)
            if freelist_size:
                freelist_name = code.name_in_slot_module_state(
                    scope.mangle_internal(Naming.freelist_name))
                freecount_name = code.name_in_slot_module_state(
                    scope.mangle_internal(Naming.freecount_name))

                type = scope.parent_type
                code.putln("#if CYTHON_USE_FREELISTS")
                self.generate_freelist_condition(
                    code, f"{freecount_name} < {freelist_size}",
                    "Py_TYPE(o)", type)
                code.putln("{")
                code.putln("%s[%s++] = %s;" % (
                    freelist_name,
                    freecount_name,
                    type.cast_code("o")))
                # Deliberately don't DECREF the type object for objects returned to the freelist:
                # we hold a reference to the type to allow them to be cleaned up properly.
                code.putln("} else")
                code.putln("#endif")
                code.putln("{")
            code.putln("PyTypeObject *tp = Py_TYPE(o);")
            code.putln("#if CYTHON_USE_TYPE_SLOTS")
            # Asking for PyType_GetSlot(..., Py_tp_free) seems to cause an error in pypy
            code.putln("(*tp->tp_free)(o);")
            code.putln("#else")
            code.putln("{")
            code.putln("freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);")
            code.putln("if (tp_free) tp_free(o);")
            code.putln("}")
            code.putln("#endif")
            code.putln("#if CYTHON_USE_TYPE_SPECS")
            # Undo the INCREF of the type object in tp_new
            code.putln("Py_DECREF(tp);")
            code.putln("#endif")
            if freelist_size:
                code.putln("}")

        if needs_trashcan:
            code.putln("__Pyx_TRASHCAN_END")

        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_usr_dealloc_call(self, scope, code):
        entry = scope.lookup_here("__dealloc__")
        if not entry or not entry.is_special:
            return

        code.putln("{")
        code.putln("PyObject *etype, *eval, *etb;")
        code.putln("PyErr_Fetch(&etype, &eval, &etb);")
        # increase the refcount while we are calling into user code
        # to prevent recursive deallocation
        code.putln("__Pyx_SET_REFCNT(o, Py_REFCNT(o) + 1);")
        code.putln("%s(o);" % entry.func_cname)
        code.putln("__Pyx_SET_REFCNT(o, Py_REFCNT(o) - 1);")
        code.putln("PyErr_Restore(etype, eval, etb);")
        code.putln("}")

    def generate_traverse_function(self, scope, code, cclass_entry):
        tp_slot = TypeSlots.GCDependentSlot("tp_traverse")
        slot_func = scope.mangle_internal("tp_traverse")
        base_type = scope.parent_type.base_type
        if tp_slot.slot_code(scope) != slot_func:
            return  # never used

        code.start_slotfunc(scope, PyrexTypes.c_returncode_type, "tp_traverse", "PyObject *o, visitproc v, void *a")

        have_entries, (py_attrs, py_buffers, memoryview_slices) = (
            scope.get_refcounted_entries(include_gc_simple=False))

        needs_type_traverse = not base_type
        # we don't know statically if we need to traverse the type
        maybe_needs_type_traverse = False

        code.putln("int e;")

        if py_attrs or py_buffers:
            self.generate_self_cast(scope, code)

        if base_type:
            # want to call it explicitly if possible so inlining can be performed
            static_call = TypeSlots.get_base_slot_function(scope, tp_slot)
            if static_call:
                code.putln("e = %s(o, v, a); if (e) return e;" % static_call)
                # No need to call type traverse - base class will do it
            elif base_type.is_builtin_type:
                base_cname = base_type.typeptr_cname
                code.putln("{")
                code.putln(
                    f"traverseproc traverse = __Pyx_PyType_GetSlot({base_cname}, tp_traverse, traverseproc);")
                code.putln("if (!traverse); else { e = traverse(o,v,a); if (e) return e; }")
                code.putln("}")
                maybe_needs_type_traverse = True
            else:
                # This is an externally defined type.  Calling through the
                # cimported base type pointer directly interacts badly with
                # the module cleanup, which may already have cleared it.
                # In that case, fall back to traversing the type hierarchy.
                # If we're using the module state then always go through the
                # type hierarchy, because our access to the module state may
                # have been lost (at least for the limited API version of
                # using module state).
                base_cname = code.name_in_slot_module_state(base_type.typeptr_cname)
                code.putln("#if !CYTHON_USE_MODULE_STATE")
                code.putln("e = 0;")
                code.putln("if (likely(%s)) {" % base_cname)
                code.putln(
                    f"traverseproc traverse = __Pyx_PyType_GetSlot({base_cname}, tp_traverse, traverseproc);")
                code.putln("if (traverse) { e = traverse(o, v, a); }")
                code.putln("} else")
                code.putln("#endif")
                code.putln("{ e = __Pyx_call_next_tp_traverse(o, v, a, %s); }" % slot_func)
                code.putln("if (e) return e;")
                code.globalstate.use_utility_code(
                    UtilityCode.load_cached("CallNextTpTraverse", "ExtensionTypes.c"))
                maybe_needs_type_traverse = True
        if needs_type_traverse or maybe_needs_type_traverse:
            code.putln("{")
            code.putln(f"e = __Pyx_call_type_traverse(o, {int(not maybe_needs_type_traverse)}, v, a);")
            code.putln("if (e) return e;")
            code.putln("}")
            code.globalstate.use_utility_code(
                UtilityCode.load_cached("CallTypeTraverse", "ExtensionTypes.c"))

        for entry in py_attrs:
            var_code = "p->%s" % entry.cname
            var_as_pyobject = PyrexTypes.typecast(py_object_type, entry.type, var_code)
            code.putln("if (%s) {" % var_code)
            code.putln("e = (*v)(%s, a); if (e) return e;" % var_as_pyobject)
            code.putln("}")

        # Traverse buffer exporting objects.
        # Note: not traversing memoryview attributes of memoryview slices!
        # When triggered by the GC, it would cause multiple visits (gc_refs
        # subtractions which is not matched by its reference count!)
        for entry in py_buffers:
            cname = entry.cname + ".obj"
            code.putln("if (p->%s) {" % cname)
            code.putln("e = (*v)(p->%s, a); if (e) return e;" % cname)
            code.putln("}")

        code.putln("return 0;")
        code.putln("}")
        code.exit_cfunc_scope()

    def generate_clear_function(self, scope, code, cclass_entry):
        tp_slot = TypeSlots.get_slot_by_name("tp_clear", scope.directives)
        slot_func = scope.mangle_internal("tp_clear")
        base_type = scope.parent_type.base_type
        if tp_slot.slot_code(scope) != slot_func:
            return  # never used

        have_entries, (py_attrs, py_buffers, memoryview_slices) = (
            scope.get_refcounted_entries(include_gc_simple=False))

        if py_attrs or py_buffers or base_type:
            unused = ''
        else:
            unused = 'CYTHON_UNUSED '

        code.start_slotfunc(scope, PyrexTypes.c_returncode_type, "tp_clear", f"{unused}PyObject *o")

        if py_attrs and Options.clear_to_none:
            code.putln("PyObject* tmp;")

        if py_attrs or py_buffers:
            self.generate_self_cast(scope, code)

        if base_type:
            # want to call it explicitly if possible so inlining can be performed
            static_call = TypeSlots.get_base_slot_function(scope, tp_slot)
            if static_call:
                code.putln("%s(o);" % static_call)
            elif base_type.is_builtin_type:
                base_cname = base_type.typeptr_cname
                code.putln("{")
                code.putln(f"inquiry clear = __Pyx_PyType_GetSlot({base_cname}, tp_clear, inquiry);")
                code.putln("if (clear) clear(o);")
                code.putln("}")
            else:
                # This is an externally defined type.  Calling through the
                # cimported base type pointer directly interacts badly with
                # the module cleanup, which may already have cleared it.
                # In that case, fall back to traversing the type hierarchy.
                # If we're using the module state then always go through the
                # type hierarchy, because our access to the module state may
                # have been lost (at least for the limited API version of
                # using module state).
                base_cname = code.name_in_slot_module_state(base_type.typeptr_cname)
                code.putln("#if !CYTHON_USE_MODULE_STATE")
                code.putln("if (likely(%s)) {" % base_cname)
                code.putln(f"inquiry clear = __Pyx_PyType_GetSlot({base_cname}, tp_clear, inquiry);")
                code.putln("if (clear) clear(o);")
                code.putln("} else")
                code.putln("#endif")
                code.putln("{ __Pyx_call_next_tp_clear(o, %s); }" % slot_func)
                code.globalstate.use_utility_code(
                    UtilityCode.load_cached("CallNextTpClear", "ExtensionTypes.c"))

        if Options.clear_to_none:
            for entry in py_attrs:
                name = "p->%s" % entry.cname
                code.putln("tmp = ((PyObject*)%s);" % name)
                if entry.is_declared_generic:
                    code.put_init_to_py_none(name, py_object_type, nanny=False)
                else:
                    code.put_init_to_py_none(name, entry.type, nanny=False)
                code.putln("Py_XDECREF(tmp);")
        else:
            for entry in py_attrs:
                code.putln("Py_CLEAR(p->%s);" % entry.cname)

        for entry in py_buffers:
            # Note: shouldn't this call PyBuffer_Release ??
            code.putln("Py_CLEAR(p->%s.obj);" % entry.cname)

        if cclass_entry.cname == '__pyx_memoryviewslice':
            code.putln("__PYX_XCLEAR_MEMVIEW(&p->from_slice, 1);")

        code.putln("return 0;")
        code.putln("}")
        code.exit_cfunc_scope()

    def generate_getitem_int_function(self, scope, code):
        # This function is put into the sq_item slot when
        # a __getitem__ method is present. It converts its
        # argument to a Python integer and calls mp_subscript.
        code.start_slotfunc(scope, PyrexTypes.py_objptr_type, "sq_item", "PyObject *o, Py_ssize_t i", needs_funcstate=False)
        code.putln(
            "PyObject *r;")
        code.putln(
            "PyObject *x = PyLong_FromSsize_t(i); if(!x) return 0;")
        # Note that PyType_GetSlot only works on heap-types before 3.10, so not using type slots
        # and defining cdef classes as non-heap types is probably impossible
        code.putln("#if CYTHON_USE_TYPE_SLOTS || (!CYTHON_USE_TYPE_SPECS && __PYX_LIMITED_VERSION_HEX < 0x030A0000)")
        code.putln(
            "r = Py_TYPE(o)->tp_as_mapping->mp_subscript(o, x);")
        code.putln("#else")
        code.putln("r = ((binaryfunc)PyType_GetSlot(Py_TYPE(o), Py_mp_subscript))(o, x);")
        code.putln("#endif")
        code.putln(
            "Py_DECREF(x);")
        code.putln(
            "return r;")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_ass_subscript_function(self, scope, code):
        # Setting and deleting an item are both done through
        # the ass_subscript method, so we dispatch to user's __setitem__
        # or __delitem__, or raise an exception.
        base_type = scope.parent_type.base_type
        set_entry = scope.lookup_here("__setitem__")
        del_entry = scope.lookup_here("__delitem__")
        code.start_slotfunc(scope, PyrexTypes.c_returncode_type, "mp_ass_subscript", "PyObject *o, PyObject *i, PyObject *v")
        code.putln(
            "if (v) {")
        if set_entry:
            code.putln("return %s(o, i, v);" % set_entry.func_cname)
        else:
            code.putln(
                "__Pyx_TypeName o_type_name;")
            self.generate_guarded_basetype_call(
                base_type, "tp_as_mapping", "mp_ass_subscript", "objobjargproc", "o, i, v", code)
            code.putln(
                "o_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(o));")
            code.putln(
                "PyErr_Format(PyExc_NotImplementedError,")
            code.putln(
                '  "Subscript assignment not supported by " __Pyx_FMT_TYPENAME, o_type_name);')
            code.putln(
                "__Pyx_DECREF_TypeName(o_type_name);")
            code.putln(
                "return -1;")
        code.putln(
            "}")
        code.putln(
            "else {")
        if del_entry:
            code.putln(
                "return %s(o, i);" % (
                    del_entry.func_cname))
        else:
            code.putln(
                "__Pyx_TypeName o_type_name;")
            self.generate_guarded_basetype_call(
                base_type, "tp_as_mapping", "mp_ass_subscript", "objobjargproc", "o, i, v", code)
            code.putln(
                "o_type_name = __Pyx_PyType_GetFullyQualifiedName(Py_TYPE(o));")
            code.putln(
                "PyErr_Format(PyExc_NotImplementedError,")
            code.putln(
                '  "Subscript deletion not supported by " __Pyx_FMT_TYPENAME, o_type_name);')
            code.putln(
                "__Pyx_DECREF_TypeName(o_type_name);")
            code.putln(
                "return -1;")
        code.putln(
            "}")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_guarded_basetype_call(
            self, base_type, substructure, slot, functype, args, code):
        if base_type:
            base_tpname = code.typeptr_cname_in_module_state(base_type)
            # Note that the limited API versions will only work for non-heaptypes on Python3.10+.
            # I think that's unavoidable and the best we can do.
            if substructure:
                code.putln(
                    f"{functype} f = __Pyx_PyType_TryGetSubSlot({base_tpname}, {substructure}, {slot}, {functype});")
            else:
                code.putln(
                    f"{functype} f = __Pyx_PyType_TryGetSlot({base_tpname}, {slot}, {functype});")
            code.putln("if (f)")
            code.putln(f"return f({args});")

    def generate_richcmp_function(self, scope, code):
        if scope.lookup_here("__richcmp__"):
            # user implemented, nothing to do
            return
        # otherwise, we have to generate it from the Python special methods
        code.start_slotfunc(scope, PyrexTypes.py_objptr_type, "tp_richcompare", "PyObject *o1, PyObject *o2, int op")
        code.putln("switch (op) {")

        class_scopes = []
        cls = scope.parent_type
        while cls is not None and not cls.entry.visibility == 'extern':
            class_scopes.append(cls.scope)
            cls = cls.scope.parent_type.base_type
        assert scope in class_scopes

        extern_parent = None
        if cls and cls.entry.visibility == 'extern':
            # need to call up into base classes as we may not know all implemented comparison methods
            extern_parent = cls if cls.typeptr_cname else scope.parent_type.base_type

        total_ordering = 'total_ordering' in scope.directives

        comp_entry = {}

        for cmp_method in TypeSlots.richcmp_special_methods:
            for class_scope in class_scopes:
                entry = class_scope.lookup_here(cmp_method)
                if entry is not None:
                    comp_entry[cmp_method] = entry
                    break

        if total_ordering:
            # Check this is valid - we must have at least 1 operation defined.
            comp_names = [from_name for from_name, to_name in TOTAL_ORDERING if from_name in comp_entry]
            if not comp_names:
                if '__eq__' not in comp_entry and '__ne__'  not in comp_entry:
                    warning(scope.parent_type.pos,
                            "total_ordering directive used, but no comparison and equality methods defined")
                else:
                    warning(scope.parent_type.pos,
                          "total_ordering directive used, but no comparison methods defined")
                total_ordering = False
            else:
                if '__eq__' not in comp_entry and '__ne__' not in comp_entry:
                    warning(scope.parent_type.pos, "total_ordering directive used, but no equality method defined")
                    total_ordering = False

                # Same priority as functools, prefers
                # __lt__ to __le__ to __gt__ to __ge__
                ordering_source = max(comp_names)

        for cmp_method in TypeSlots.richcmp_special_methods:
            cmp_type = cmp_method.strip('_').upper()  # e.g. "__eq__" -> EQ
            entry = comp_entry.get(cmp_method)
            if entry is None and (not total_ordering or cmp_type in ('NE', 'EQ')):
                # No definition, fall back to superclasses.
                # eq/ne methods shouldn't use the total_ordering code.
                continue

            code.putln("case Py_%s: {" % cmp_type)
            if entry is None:
                assert total_ordering
                # We need to generate this from the other methods.
                invert_comp, comp_op, invert_equals = TOTAL_ORDERING[ordering_source, cmp_method]

                # First we always do the comparison.
                code.putln("PyObject *ret;")
                code.putln("ret = %s(o1, o2);" % comp_entry[ordering_source].func_cname)
                code.putln("if (likely(ret && ret != Py_NotImplemented)) {")
                code.putln("int order_res = __Pyx_PyObject_IsTrue(ret);")
                code.putln("Py_DECREF(ret);")
                code.putln("if (unlikely(order_res < 0)) return NULL;")
                # We may need to check equality too. For some combos it's never required.
                if invert_equals is not None:
                    # Implement the and/or check with an if.
                    if comp_op == '&&':
                        code.putln("if (%s order_res) {" % ('!!' if invert_comp else '!'))
                        code.putln("ret = __Pyx_NewRef(Py_False);")
                        code.putln("} else {")
                    elif comp_op == '||':
                        code.putln("if (%s order_res) {" % ('!' if invert_comp else ''))
                        code.putln("ret = __Pyx_NewRef(Py_True);")
                        code.putln("} else {")
                    else:
                        raise AssertionError('Unknown op %s' % (comp_op, ))
                    if '__eq__' in comp_entry:
                        eq_func = '__eq__'
                    else:
                        # Fall back to NE, which is defined here.
                        eq_func = '__ne__'
                        invert_equals = not invert_equals

                    code.putln("ret = %s(o1, o2);" % comp_entry[eq_func].func_cname)
                    code.putln("if (likely(ret && ret != Py_NotImplemented)) {")
                    code.putln("int eq_res = __Pyx_PyObject_IsTrue(ret);")
                    code.putln("Py_DECREF(ret);")
                    code.putln("if (unlikely(eq_res < 0)) return NULL;")
                    if invert_equals:
                        code.putln("ret = eq_res ? Py_False : Py_True;")
                    else:
                        code.putln("ret = eq_res ? Py_True : Py_False;")
                    code.putln("Py_INCREF(ret);")
                    code.putln("}")  # equals success
                    code.putln("}")  # Needs to try equals
                else:
                    # Convert direct to a boolean.
                    if invert_comp:
                        code.putln("ret = order_res ? Py_False : Py_True;")
                    else:
                        code.putln("ret = order_res ? Py_True : Py_False;")
                    code.putln("Py_INCREF(ret);")
                code.putln("}")  # comp_op
                code.putln("return ret;")
            else:
                code.putln("return %s(o1, o2);" % entry.func_cname)
            code.putln("}")  # Case

        if '__eq__' in comp_entry and '__ne__' not in comp_entry and not extern_parent:
            code.putln("case Py_NE: {")
            code.putln("PyObject *ret;")
            # Python itself does not do this optimisation, it seems...
            #code.putln("if (o1 == o2) return __Pyx_NewRef(Py_False);")
            code.putln("ret = %s(o1, o2);" % comp_entry['__eq__'].func_cname)
            code.putln("if (likely(ret && ret != Py_NotImplemented)) {")
            code.putln("int b = __Pyx_PyObject_IsTrue(ret);")
            code.putln("Py_DECREF(ret);")
            code.putln("if (unlikely(b < 0)) return NULL;")
            code.putln("ret = (b) ? Py_False : Py_True;")
            code.putln("Py_INCREF(ret);")
            code.putln("}")
            code.putln("return ret;")
            code.putln("}")

        code.putln("default: {")
        if extern_parent and extern_parent.typeptr_cname:
            code.putln("if (likely(%s->tp_richcompare)) return %s->tp_richcompare(o1, o2, op);" % (
                extern_parent.typeptr_cname, extern_parent.typeptr_cname))
        code.putln("return __Pyx_NewRef(Py_NotImplemented);")
        code.putln("}")

        code.putln("}")  # switch
        code.putln("}")
        code.exit_cfunc_scope()

    def generate_binop_function(self, scope, slot, code, pos):
        func_name = scope.mangle_internal(slot.slot_name)
        if scope.directives['c_api_binop_methods']:
            code.putln('#define %s %s' % (func_name, slot.left_slot.slot_code(scope)))
            return

        if slot.left_slot.signature in (TypeSlots.binaryfunc, TypeSlots.ibinaryfunc):
            slot_type = 'binaryfunc'
            extra_arg = extra_arg_decl = ''
        elif slot.left_slot.signature in (TypeSlots.powternaryfunc, TypeSlots.ipowternaryfunc):
            slot_type = 'ternaryfunc'
            extra_arg = ', extra_arg'
            extra_arg_decl = ', PyObject* extra_arg'
        else:
            error(pos, "Unexpected type slot signature: %s" % slot)
            return

        def get_slot_method_cname(method_name):
            entry = scope.lookup(method_name)
            return entry.func_cname if entry and entry.is_special else None

        def call_slot_method(method_name, reverse):
            func_cname = get_slot_method_cname(method_name)
            if func_cname:
                return "%s(%s%s)" % (
                    func_cname,
                    "right, left" if reverse else "left, right",
                    extra_arg)
            else:
                return '%s_maybe_call_slot(__Pyx_PyType_GetSlot(%s, tp_base, PyTypeObject*), left, right %s)' % (
                    func_name,
                    code.name_in_module_state(scope.parent_type.typeptr_cname),
                    extra_arg)

        if get_slot_method_cname(slot.left_slot.method_name) and not get_slot_method_cname(slot.right_slot.method_name):
            warning(pos, "Extension type implements %s() but not %s(). "
                         "The behaviour has changed from previous Cython versions to match Python semantics. "
                         "You can implement both special methods in a backwards compatible way." % (
                slot.left_slot.method_name,
                slot.right_slot.method_name,
            ))

        code.putln()
        preprocessor_guard = slot.preprocessor_guard_code()
        if preprocessor_guard:
            code.putln(preprocessor_guard)
        code.enter_cfunc_scope(scope)  # C class scope, not function scope

        overloads_left = int(bool(get_slot_method_cname(slot.left_slot.method_name)))
        overloads_right = int(bool(get_slot_method_cname(slot.right_slot.method_name)))
        parent_type_cname = scope.parent_type.typeptr_cname
        if scope.parent_type.is_extension_type:
            parent_type_cname = code.name_in_module_state(parent_type_cname)
        code.putln(
            TempitaUtilityCode.load_as_string(
                "BinopSlot", "ExtensionTypes.c",
                context={
                    "func_name": func_name,
                    "slot_name": slot.slot_name,
                    "overloads_left": overloads_left,
                    "overloads_right": overloads_right,
                    "call_left": call_slot_method(slot.left_slot.method_name, reverse=False),
                    "call_right": call_slot_method(slot.right_slot.method_name, reverse=True),
                    "type_cname": parent_type_cname,
                    "slot_type": slot_type,
                    "extra_arg": extra_arg,
                    "extra_arg_decl": extra_arg_decl,
                    })[1])

        code.exit_cfunc_scope()
        if preprocessor_guard:
            code.putln("#endif")

    def generate_getattro_function(self, scope, code):
        # First try to get the attribute using __getattribute__, if defined, or
        # PyObject_GenericGetAttr.
        #
        # If that raises an AttributeError, call the __getattr__ if defined.
        #
        # In both cases, defined can be in this class, or any base class.
        def lookup_here_or_base(n, tp=None, extern_return=None):
            # Recursive lookup
            if tp is None:
                tp = scope.parent_type
            r = tp.scope.lookup_here(n)
            if r is None:
                if tp.is_external and extern_return is not None:
                    return extern_return
                if tp.base_type is not None:
                    return lookup_here_or_base(n, tp.base_type)
            return r

        getattr_entry = lookup_here_or_base("__getattr__")
        getattribute_entry = lookup_here_or_base("__getattribute__")

        code.start_slotfunc(scope, PyrexTypes.py_objptr_type, "tp_getattro", "PyObject *o, PyObject *n", needs_funcstate=False)
        if getattribute_entry is not None:
            code.putln(
                "PyObject *v = %s(o, n);" % (
                    getattribute_entry.func_cname))
        else:
            code.putln(
                "PyObject *v = PyObject_GenericGetAttr(o, n);")
        if getattr_entry is not None:
            code.putln(
                "if (!v && PyErr_ExceptionMatches(PyExc_AttributeError)) {")
            code.putln(
                "PyErr_Clear();")
            code.putln(
                "v = %s(o, n);" % (
                    getattr_entry.func_cname))
            code.putln(
                "}")
        code.putln(
            "return v;")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_setattro_function(self, scope, code):
        # Setting and deleting an attribute are both done through
        # the setattro method, so we dispatch to user's __setattr__
        # or __delattr__ or fall back on PyObject_GenericSetAttr.
        base_type = scope.parent_type.base_type
        set_entry = scope.lookup_here("__setattr__")
        del_entry = scope.lookup_here("__delattr__")

        code.start_slotfunc(scope, PyrexTypes.c_returncode_type, "tp_setattro", "PyObject *o, PyObject *n, PyObject *v")
        code.putln(
            "if (v) {")
        if set_entry:
            code.putln(
                "return %s(o, n, v);" % (
                    set_entry.func_cname))
        else:
            self.generate_guarded_basetype_call(
                base_type, None, "tp_setattro", "setattrofunc", "o, n, v", code)
            code.putln(
                "return PyObject_GenericSetAttr(o, n, v);")
        code.putln(
            "}")
        code.putln(
            "else {")
        if del_entry:
            code.putln(
                "return %s(o, n);" % (
                    del_entry.func_cname))
        else:
            self.generate_guarded_basetype_call(
                base_type, None, "tp_setattro", "setattrofunc", "o, n, v", code)
            code.putln(
                "return PyObject_GenericSetAttr(o, n, 0);")
        code.putln(
            "}")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_descr_get_function(self, scope, code):
        # The __get__ function of a descriptor object can be
        # called with NULL for the second or third arguments
        # under some circumstances, so we replace them with
        # None in that case.
        user_get_entry = scope.lookup_here("__get__")

        code.start_slotfunc(scope, PyrexTypes.py_objptr_type, "tp_descr_get", "PyObject *o, PyObject *i, PyObject *c", needs_funcstate=False)
        code.putln(
            "PyObject *r = 0;")
        code.putln(
            "if (!i) i = Py_None;")
        code.putln(
            "if (!c) c = Py_None;")
        #code.put_incref("i", py_object_type)
        #code.put_incref("c", py_object_type)
        code.putln(
            "r = %s(o, i, c);" % (
                user_get_entry.func_cname))
        #code.put_decref("i", py_object_type)
        #code.put_decref("c", py_object_type)
        code.putln(
            "return r;")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_descr_set_function(self, scope, code):
        # Setting and deleting are both done through the __set__
        # method of a descriptor, so we dispatch to user's __set__
        # or __delete__ or raise an exception.
        base_type = scope.parent_type.base_type
        user_set_entry = scope.lookup_here("__set__")
        user_del_entry = scope.lookup_here("__delete__")

        code.start_slotfunc(scope, PyrexTypes.c_returncode_type, "tp_descr_set", "PyObject *o, PyObject *i, PyObject *v")
        code.putln(
            "if (v) {")
        if user_set_entry:
            code.putln(
                "return %s(o, i, v);" % (
                    user_set_entry.func_cname))
        else:
            self.generate_guarded_basetype_call(
                base_type, None, "tp_descr_set", "descrsetfunc", "o, i, v", code)
            code.putln(
                'PyErr_SetString(PyExc_NotImplementedError, "__set__");')
            code.putln(
                "return -1;")
        code.putln(
            "}")
        code.putln(
            "else {")
        if user_del_entry:
            code.putln(
                "return %s(o, i);" % (
                    user_del_entry.func_cname))
        else:
            self.generate_guarded_basetype_call(
                base_type, None, "tp_descr_set", "descrsetfunc", "o, i, v", code)
            code.putln(
                'PyErr_SetString(PyExc_NotImplementedError, "__delete__");')
            code.putln(
                "return -1;")
        code.putln(
            "}")
        code.putln(
            "}")
        code.exit_cfunc_scope()

    def generate_property_accessors(self, cclass_scope, code):
        for entry in cclass_scope.property_entries:
            property_scope = entry.scope
            if property_scope.defines_any(["__get__"]):
                self.generate_property_get_function(entry, code)
            if property_scope.defines_any(["__set__", "__del__"]):
                self.generate_property_set_function(entry, code)

    def generate_property_get_function(self, property_entry, code):
        property_scope = property_entry.scope
        property_entry.getter_cname = property_scope.parent_scope.mangle(
            Naming.prop_get_prefix, property_entry.name)
        get_entry = property_scope.lookup_here("__get__")

        code.putln("")
        code.putln(
            "static PyObject *%s(PyObject *o, CYTHON_UNUSED void *x) {" % (
                property_entry.getter_cname))
        code.putln(
            "return %s(o);" % (
                get_entry.func_cname))
        code.putln(
            "}")

    def generate_property_set_function(self, property_entry, code):
        property_scope = property_entry.scope
        property_entry.setter_cname = property_scope.parent_scope.mangle(
            Naming.prop_set_prefix, property_entry.name)
        set_entry = property_scope.lookup_here("__set__")
        del_entry = property_scope.lookup_here("__del__")

        code.putln("")
        code.putln(
            "static int %s(PyObject *o, PyObject *v, CYTHON_UNUSED void *x) {" % (
                property_entry.setter_cname))
        code.putln(
            "if (v) {")
        if set_entry:
            code.putln(
                "return %s(o, v);" % (
                    set_entry.func_cname))
        else:
            code.putln(
                'PyErr_SetString(PyExc_NotImplementedError, "__set__");')
            code.putln(
                "return -1;")
        code.putln(
            "}")
        code.putln(
            "else {")
        if del_entry:
            code.putln(
                "return %s(o);" % (
                    del_entry.func_cname))
        else:
            code.putln(
                'PyErr_SetString(PyExc_NotImplementedError, "__del__");')
            code.putln(
                "return -1;")
        code.putln(
            "}")
        code.putln(
            "}")

    def generate_typeobj_spec(self, entry, code):
        ext_type = entry.type
        scope = ext_type.scope

        members_slot = TypeSlots.get_slot_by_name("tp_members", code.globalstate.directives)
        members_slot.generate_substructure_spec(scope, code)

        buffer_slot = TypeSlots.get_slot_by_name("tp_as_buffer", code.globalstate.directives)
        if not buffer_slot.is_empty(scope):
            code.putln("#if !CYTHON_COMPILING_IN_LIMITED_API")
            buffer_slot.generate_substructure(scope, code)
            code.putln("#endif")

        if ext_type.typedef_flag:
            objstruct = ext_type.objstruct_cname
        else:
            objstruct = "struct %s" % ext_type.objstruct_cname

        weakref_entry = scope.lookup_here("__weakref__") if not scope.is_closure_class_scope else None
        if weakref_entry and weakref_entry.is_inherited:
            weakref_entry = None  # only generate it for the defining class
        generate_members = bool(weakref_entry)
        if generate_members:
            code.globalstate.use_utility_code(
                UtilityCode.load_cached("IncludeStructmemberH", "ModuleSetupCode.c"))
            code.putln("static PyMemberDef %s_members[] = {" % ext_type.typeobj_cname)
            code.putln("#if !CYTHON_USE_TYPE_SLOTS")
            if weakref_entry:
                # Note that unlike the assignment of tp_weaklistoffset in the type-ready code
                # used in the non-limited API case, this doesn't preserve the weaklistoffset
                # from base classes.
                # Practically that doesn't matter, but it isn't exactly the identical.
                code.putln('{"__weaklistoffset__", T_PYSSIZET, offsetof(%s, %s), READONLY, 0},'
                           % (objstruct, weakref_entry.cname))
            code.putln("#endif")
            code.putln("{0, 0, 0, 0, 0}")
            code.putln("};")

            if weakref_entry:
                position = format_position(weakref_entry.pos)
                weakref_warn_mesage = (
                    f"{position}: __weakref__ is unsupported in the Limited API when "
                    "running on Python <3.9.")
                # Note: Limited API rather than USE_TYPE_SPECS - we work round the issue
                # with USE_TYPE_SPECS outside the limited API
                code.putln("#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x03090000")
                code.putln("#if defined(__GNUC__) || defined(__clang__)")
                code.putln(f'#warning "{weakref_warn_mesage}"')
                code.putln("#elif defined(_MSC_VER)")
                code.putln(f'#pragma message("{weakref_warn_mesage}")')
                code.putln("#endif")
                code.putln("#endif")


        code.putln("static PyType_Slot %s_slots[] = {" % ext_type.typeobj_cname)
        for slot in TypeSlots.get_slot_table(code.globalstate.directives):
            slot.generate_spec(scope, code)
        if generate_members:
            code.putln("{Py_tp_members, (void*)%s_members}," % ext_type.typeobj_cname)
        code.putln("{0, 0},")
        code.putln("};")

        classname = scope.class_name.as_c_string_literal()
        code.putln("static PyType_Spec %s_spec = {" % ext_type.typeobj_cname)
        code.putln('"%s.%s",' % (self.full_module_name, classname.replace('"', '')))
        code.putln("sizeof(%s)," % objstruct)
        code.putln("0,")
        code.putln("%s," % TypeSlots.get_slot_by_name("tp_flags", scope.directives).slot_code(scope))
        code.putln("%s_slots," % ext_type.typeobj_cname)
        code.putln("};")

    def generate_typeobj_definition(self, modname, entry, code):
        type = entry.type
        scope = type.scope
        for suite in TypeSlots.get_slot_table(code.globalstate.directives).substructures:
            suite.generate_substructure(scope, code)
        code.putln("")
        if entry.visibility == 'public':
            header = "DL_EXPORT(PyTypeObject) %s = {"
        else:
            header = "static PyTypeObject %s = {"
        #code.putln(header % scope.parent_type.typeobj_cname)
        code.putln(header % type.typeobj_cname)
        code.putln(
            "PyVarObject_HEAD_INIT(0, 0)")
        classname = scope.class_name.as_c_string_literal()
        code.putln(
            '"%s."%s, /*tp_name*/' % (
                self.full_module_name,
                classname))
        if type.typedef_flag:
            objstruct = type.objstruct_cname
        else:
            objstruct = "struct %s" % type.objstruct_cname
        code.putln(
            "sizeof(%s), /*tp_basicsize*/" % objstruct)
        code.putln(
            "0, /*tp_itemsize*/")
        for slot in TypeSlots.get_slot_table(code.globalstate.directives):
            slot.generate(scope, code)
        code.putln(
            "};")

    def generate_method_table(self, env, code):
        if env.is_c_class_scope and not env.pyfunc_entries:
            return
        binding = env.directives['binding']

        code.putln("")
        wrapper_code_writer = code.insertion_point()

        code.putln(
            "static PyMethodDef %s[] = {" % (
                env.method_table_cname))
        for entry in env.pyfunc_entries:
            if not entry.fused_cfunction and not (binding and entry.is_overridable):
                code.put_pymethoddef(entry, ",", wrapper_code_writer=wrapper_code_writer)
        code.putln(
            "{0, 0, 0, 0}")
        code.putln(
            "};")

        if wrapper_code_writer.getvalue():
            wrapper_code_writer.putln("")

    def generate_dict_getter_function(self, scope, code):
        dict_attr = scope.lookup_here("__dict__")
        if not dict_attr or not dict_attr.is_variable:
            return
        func_name = scope.mangle_internal("__dict__getter")
        dict_name = dict_attr.cname

        code.putln("")
        code.putln("#if CYTHON_COMPILING_IN_LIMITED_API && __PYX_LIMITED_VERSION_HEX < 0x030A0000")
        code.putln("static PyObject *%s(PyObject *o, CYTHON_UNUSED void *x) {" % func_name)
        self.generate_self_cast(scope, code)
        code.putln("if (unlikely(!p->%s)){" % dict_name)
        code.putln("p->%s = PyDict_New();" % dict_name)
        code.putln("}")
        code.putln("Py_XINCREF(p->%s);" % dict_name)
        code.putln("return p->%s;" % dict_name)
        code.putln("}")
        code.putln("#else")
        # PyObject_GenericGetDict has the advantage that it's freethreading thread-safe,
        # handles both managed and unmanaged dicts (in case we switch to managed in future),
        # and can potentially do optimizations with per-class shared keys.
        code.putln(f"#define {func_name} PyObject_GenericGetDict")
        code.putln("#endif")

    def generate_getset_table(self, env, code):
        if env.property_entries:
            code.putln("")
            code.putln(
                "static struct PyGetSetDef %s[] = {" %
                env.getset_table_cname)
            for entry in env.property_entries:
                doc = entry.doc
                if doc:
                    if doc.is_unicode:
                        doc = doc.as_utf8_string()
                    doc_code = "PyDoc_STR(%s)" % doc.as_c_string_literal()
                else:
                    doc_code = "0"
                code.putln(
                    '{%s, %s, %s, %s, 0},' % (
                        entry.name.as_c_string_literal(),
                        entry.getter_cname or "0",
                        entry.setter_cname or "0",
                        doc_code))
            code.putln(
                "{0, 0, 0, 0, 0}")
            code.putln(
                "};")

    def create_import_star_conversion_utility_code(self, env):
        # Create all conversion helpers that are needed for "import *" assignments.
        # Must be done before code generation to support CythonUtilityCode.
        for name, entry in sorted(env.entries.items()):
            if entry.is_cglobal and entry.used:
                if not entry.type.is_pyobject:
                    entry.type.create_from_py_utility_code(env)

    def generate_import_star(self, env, code):
        env.use_utility_code(UtilityCode.load_cached("CStringEquals", "StringTools.c"))
        code.start_initcfunc(
            f"int {Naming.import_star_set}("
            f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname},"
            "PyObject *o, PyObject* py_name, const char *name)")

        code.putln("static const char* internal_type_names[] = {")
        for name, entry in sorted(env.entries.items()):
            if entry.is_type:
                code.putln('"%s",' % name)
        code.putln("0")
        code.putln("};")

        code.putln("const char** type_name = internal_type_names;")
        code.putln("while (*type_name) {")
        code.putln("if (__Pyx_StrEq(name, *type_name)) {")
        code.putln('PyErr_Format(PyExc_TypeError, "Cannot overwrite C type %s", name);')
        code.putln('goto bad;')
        code.putln("}")
        code.putln("type_name++;")
        code.putln("}")

        old_error_label = code.new_error_label()
        code.putln("if (0);")  # so the first one can be "else if"
        msvc_count = 0
        for name, entry in sorted(env.entries.items()):
            if entry.is_cglobal and entry.used and not entry.type.is_const:
                msvc_count += 1
                if msvc_count % 100 == 0:
                    code.putln("#ifdef _MSC_VER")
                    code.putln("if (0);  /* Workaround for MSVC C1061. */")
                    code.putln("#endif")
                code.putln('else if (__Pyx_StrEq(name, "%s")) {' % name)
                if entry.type.is_pyobject:
                    if entry.type.is_extension_type or entry.type.is_builtin_type:
                        type_test = entry.type.type_test_code(
                            env, "o")
                        code.putln("if (!(%s)) %s;" % (
                            type_test,
                            code.error_goto(entry.pos)))
                    code.putln("Py_INCREF(o);")
                    code.put_decref(entry.cname, entry.type, nanny=False)
                    code.putln("%s = %s;" % (
                        entry.cname,
                        PyrexTypes.typecast(entry.type, py_object_type, "o")))
                elif entry.type.create_from_py_utility_code(env):
                    # if available, utility code was already created in self.prepare_utility_code()
                    code.putln(entry.type.from_py_call_code(
                        'o', entry.cname, entry.pos, code))
                else:
                    code.putln('PyErr_Format(PyExc_TypeError, "Cannot convert Python object %s to %s");' % (
                        name, entry.type))
                    code.putln(code.error_goto(entry.pos))
                code.putln("}")
        code.putln("else {")
        code.putln("if (PyObject_SetAttr(%s, py_name, o) < 0) goto bad;" % Naming.module_cname)
        code.putln("}")
        code.putln("return 0;")
        if code.label_used(code.error_label):
            code.put_label(code.error_label)
            # This helps locate the offending name.
            code.put_add_traceback(EncodedString(self.full_module_name))
        code.error_label = old_error_label
        code.putln("bad:")
        code.putln("return -1;")
        code.putln("}")
        code.putln("")
        code.put_code_here(UtilityCode.load("ImportStar", "ImportExport.c"))
        code.exit_cfunc_scope()  # done with labels

    def generate_module_state_start(self, env, code):
        # TODO: Refactor to move module state struct decl closer to the static decl
        code.putln('typedef struct {')
        code.putln('PyObject *%s;' % env.module_dict_cname)
        code.putln('PyObject *%s;' % Naming.builtins_cname)
        code.putln('PyObject *%s;' % Naming.cython_runtime_cname)
        code.putln('PyObject *%s;' % Naming.empty_tuple)
        code.putln('PyObject *%s;' % Naming.empty_bytes)
        code.putln('PyObject *%s;' % Naming.empty_unicode)
        if Options.pre_import is not None:
            code.putln('PyObject *%s;' % Naming.preimport_cname)

    def generate_module_state_end(self, env, modules, globalstate):
        module_state = globalstate['module_state_end']
        module_state_clear = globalstate['module_state_clear_end']
        module_state_traverse = globalstate['module_state_traverse_end']
        module_state.putln('} %s;' % Naming.modulestatetype_cname)
        module_state.putln('')
        globalstate.use_utility_code(
            UtilityCode.load("MultiPhaseInitModuleState", "ModuleSetupCode.c")
        )
        module_state.putln("#if CYTHON_USE_MODULE_STATE")
        module_state.putln('#ifdef __cplusplus')
        module_state.putln('namespace {')
        module_state.putln('extern struct PyModuleDef %s;' % Naming.pymoduledef_cname)
        module_state.putln('} /* anonymous namespace */')
        module_state.putln('#else')
        module_state.putln('static struct PyModuleDef %s;' % Naming.pymoduledef_cname)
        module_state.putln('#endif')
        module_state.putln('')
        module_state.putln('#define %s (__Pyx_PyModule_GetState(__Pyx_State_FindModule(&%s)))' % (
            Naming.modulestateglobal_cname,
            Naming.pymoduledef_cname))
        module_state.putln('')
        module_state.putln('#define %s (__Pyx_State_FindModule(&%s))' % (
            env.module_cname,
            Naming.pymoduledef_cname))
        module_state.putln("#else")
        module_state.putln('static %s %s_static =' % (
            Naming.modulestatetype_cname,
            Naming.modulestateglobal_cname
        ))
        module_state.putln('#ifdef __cplusplus')
        # C++ likes to be initialized with {} to avoid "missing initializer" warnings
        # but it isn't valid C
        module_state.putln('    {};')
        module_state.putln('#else')
        module_state.putln('    {0};')
        module_state.putln('#endif')
        module_state.putln('static %s * const %s = &%s_static;' % (
            Naming.modulestatetype_cname,
            Naming.modulestateglobal_cname,
            Naming.modulestateglobal_cname
        ))
        module_state.putln("#endif")
        module_state_clear.putln("return 0;")
        module_state_clear.putln("}")
        module_state_clear.putln("#endif")
        module_state_traverse.putln("return 0;")
        module_state_traverse.putln("}")
        module_state_traverse.putln("#endif")


    def generate_module_state_clear(self, env, code):
        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln("static CYTHON_SMALL_CODE int %s_clear(PyObject *m) {" % Naming.module_cname)
        code.putln(f"{Naming.modulestatetype_cname} *clear_module_state = __Pyx_PyModule_GetState(m);")
        code.putln("if (!clear_module_state) return 0;")
        code.putln('Py_CLEAR(clear_module_state->%s);' %
            env.module_dict_cname)
        code.putln('Py_CLEAR(clear_module_state->%s);' %
            Naming.builtins_cname)
        code.putln('Py_CLEAR(clear_module_state->%s);' %
            Naming.cython_runtime_cname)
        code.putln('Py_CLEAR(clear_module_state->%s);' %
            Naming.empty_tuple)
        code.putln('Py_CLEAR(clear_module_state->%s);' %
            Naming.empty_bytes)
        code.putln('Py_CLEAR(clear_module_state->%s);' %
            Naming.empty_unicode)
        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        # In this case we have to remove the module from our lookup table ourself
        # because Python isn't going to do it.
        code.putln("__Pyx_State_RemoveModule(NULL);")
        code.putln("#endif")

    def generate_module_state_traverse(self, env, code):
        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln("static CYTHON_SMALL_CODE int %s_traverse(PyObject *m, visitproc visit, void *arg) {" % Naming.module_cname)
        code.putln(f"{Naming.modulestatetype_cname} *traverse_module_state = __Pyx_PyModule_GetState(m);")
        code.putln("if (!traverse_module_state) return 0;")
        code.putln(f'Py_VISIT(traverse_module_state->{env.module_dict_cname});')
        code.putln(f'Py_VISIT(traverse_module_state->{Naming.builtins_cname});')
        code.putln(f'Py_VISIT(traverse_module_state->{Naming.cython_runtime_cname});')
        code.putln(f'__Pyx_VISIT_CONST(traverse_module_state->{Naming.empty_tuple});')
        code.putln(f'__Pyx_VISIT_CONST(traverse_module_state->{Naming.empty_bytes});')
        code.putln(f'__Pyx_VISIT_CONST(traverse_module_state->{Naming.empty_unicode});')

    def generate_module_init_func(self, imported_modules, shared_utility_exporter, env, code):
        subfunction = self.mod_init_subfunction(self.pos, self.scope, code)

        self.generate_pymoduledef_struct(env, code)

        code.enter_cfunc_scope(self.scope)
        code.putln("")
        code.put_code_here(UtilityCode.load("PyModInitFuncType", "ModuleSetupCode.c"))

        modinit_func_name = EncodedString(f"PyInit_{env.module_name}")
        header3 = "__Pyx_PyMODINIT_FUNC %s(void)" % self.mod_init_func_cname('PyInit', env)
        # Optimise for small code size as the module init function is only executed once.
        code.putln("%s CYTHON_SMALL_CODE; /*proto*/" % header3)
        if self.scope.is_package:
            code.putln("#if !defined(CYTHON_NO_PYINIT_EXPORT) && (defined(_WIN32) || defined(WIN32) || defined(MS_WINDOWS))")
            code.putln("__Pyx_PyMODINIT_FUNC PyInit___init__(void) { return %s(); }" % (
                self.mod_init_func_cname('PyInit', env)))
            code.putln("#endif")
        # Hack for a distutils bug - https://bugs.python.org/issue39432
        # distutils attempts to make visible a slightly wrong PyInitU module name. Just create a dummy
        # function to keep it quiet
        wrong_punycode_module_name = self.wrong_punycode_module_name(env.module_name)
        if wrong_punycode_module_name:
            code.putln("#if !defined(CYTHON_NO_PYINIT_EXPORT) && (defined(_WIN32) || defined(WIN32) || defined(MS_WINDOWS))")
            code.putln("void %s(void) {} /* workaround for https://bugs.python.org/issue39432 */" % wrong_punycode_module_name)
            code.putln("#endif")
        code.putln(header3)

        # CPython 3.5+ supports multi-phase module initialisation (gives access to __spec__, __file__, etc.)
        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        code.putln("{")
        code.putln("return PyModuleDef_Init(&%s);" % Naming.pymoduledef_cname)
        code.putln("}")

        mod_create_func = UtilityCode.load("ModuleCreationPEP489", "ModuleSetupCode.c")
        code.put_code_here(mod_create_func)

        code.putln("")
        # main module init code lives in Py_mod_exec function, not in PyInit function
        code.putln("static CYTHON_SMALL_CODE int %s(PyObject *%s)" % (
            self.module_init_func_cname(),
            Naming.pymodinit_module_arg))
        code.putln("#endif")  # PEP489

        # start of module init/exec function (pre/post PEP 489)
        code.putln("{")
        code.putln('int stringtab_initialized = 0;')
        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln('int pystate_addmodule_run = 0;')
        code.putln("#endif")
        code.putln(f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname} = NULL;")

        tempdecl_code = code.insertion_point()

        profile = code.globalstate.directives['profile']
        linetrace = code.globalstate.directives['linetrace']
        if profile or linetrace:
            if linetrace:
                code.use_fast_gil_utility_code()
            code.globalstate.use_utility_code(UtilityCode.load_cached("Profile", "Profile.c"))

        code.put_declare_refcount_context()
        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        # Most extension modules simply can't deal with it, and Cython isn't ready either.
        # See issues listed here: https://docs.python.org/3/c-api/init.html#sub-interpreter-support
        code.putln("if (%s) {" % Naming.module_cname)
        # Hack: enforce single initialisation.
        code.putln("if (%s == %s) return 0;" % (
            Naming.module_cname,
            Naming.pymodinit_module_arg,
        ))
        code.putln('PyErr_SetString(PyExc_RuntimeError,'
                   ' "Module \'%s\' has already been imported. Re-initialisation is not supported.");' %
                   env.module_name.as_c_string_literal()[1:-1])
        code.putln("return -1;")
        code.putln("}")
        code.putln("#else")
        # Hack: enforce single initialisation also on reimports under different names (with PEP 3121/489).
        code.putln("if (%s) return __Pyx_NewRef(%s);" % (
            Naming.module_cname,
            Naming.module_cname,
        ))
        code.putln("#endif")

        code.putln("/*--- Module creation code ---*/")
        self.generate_module_creation_code(env, code)

        if profile or linetrace:
            tempdecl_code.put_trace_declarations()
            code.put_trace_frame_init()

        refnanny_import_code = UtilityCode.load("ImportRefnannyAPI", "ModuleSetupCode.c")
        code.put_code_here(refnanny_import_code)
        code.put_setup_refcount_context(modinit_func_name)

        code.putln("__Pyx_init_runtime_version();")
        env.use_utility_code(UtilityCode.load("CheckBinaryVersion", "ModuleSetupCode.c"))
        code.put_error_if_neg(self.pos, "__Pyx_check_binary_version("
                                        "__PYX_LIMITED_VERSION_HEX, "
                                        "__Pyx_get_runtime_version(), "
                                        "CYTHON_COMPILING_IN_LIMITED_API)"
        )

        empty_tuple = code.name_in_main_c_code_module_state(Naming.empty_tuple)
        code.putln("%s = PyTuple_New(0); %s" % (
            empty_tuple, code.error_goto_if_null(empty_tuple, self.pos)))
        empty_bytes = code.name_in_main_c_code_module_state(Naming.empty_bytes)
        code.putln("%s = PyBytes_FromStringAndSize(\"\", 0); %s" % (
            empty_bytes, code.error_goto_if_null(empty_bytes, self.pos)))
        empty_unicode = code.name_in_main_c_code_module_state(Naming.empty_unicode)
        code.putln("%s = PyUnicode_FromStringAndSize(\"\", 0); %s" % (
            empty_unicode, code.error_goto_if_null(empty_unicode, self.pos)))


        code.putln("/*--- Library function declarations ---*/")
        if env.directives['np_pythran']:
            code.put_error_if_neg(self.pos, "_import_array()")

        code.putln("/*--- Initialize various global constants etc. ---*/")
        code.put_error_if_neg(self.pos, f"__Pyx_InitConstants({Naming.modulestatevalue_cname})")
        code.putln("stringtab_initialized = 1;")
        code.put_error_if_neg(self.pos, "__Pyx_InitGlobals()")  # calls any utility code

        code.putln("if (%s) {" % self.is_main_module_flag_cname())
        code.put_error_if_neg(self.pos, 'PyObject_SetAttr(%s, %s, %s)' % (
            env.module_cname,
            code.intern_identifier(EncodedString("__name__")),
            code.intern_identifier(EncodedString("__main__"))))
        code.putln("}")

        # set up __file__ and __path__, then add the module to sys.modules
        self.generate_module_import_setup(env, code)

        if Options.cache_builtins:
            code.putln("/*--- Builtin init code ---*/")
            code.put_error_if_neg(
                self.pos,
                f"__Pyx_InitCachedBuiltins({Naming.modulestatevalue_cname})")

        code.putln("/*--- Constants init code ---*/")
        code.put_error_if_neg(
            self.pos,
            f"__Pyx_InitCachedConstants({Naming.modulestatevalue_cname})")
        # code objects come after the other globals (since they use strings and tuples)
        code.put_error_if_neg(
            self.pos,
            f"__Pyx_CreateCodeObjects({Naming.modulestatevalue_cname})")

        code.putln("/*--- Global type/function init code ---*/")

        with subfunction("Global init code") as inner_code:
            self.generate_global_init_code(env, inner_code)

        with subfunction("Variable export code") as inner_code:
            self.generate_c_variable_export_code(env, inner_code)

        with subfunction("Function export code") as inner_code:
            self.generate_c_function_export_code(env, inner_code)

        shared_utility_exporter.call_export_code(code)

        with subfunction("Type init code") as inner_code:
            self.generate_type_init_code(env, inner_code)

        with subfunction("Type import code") as inner_code:
            for module in imported_modules:
                self.generate_type_import_code_for_module(module, env, inner_code)

        with subfunction("Variable import code") as inner_code:
            for module in imported_modules:
                self.generate_c_variable_import_code_for_module(module, env, inner_code)

        with subfunction("Function import code") as inner_code:
            for module in imported_modules:
                self.specialize_fused_types(module)
                self.generate_c_function_import_code_for_module(module, env, inner_code)

        shared_utility_exporter.call_import_code(code)

        code.putln("/*--- Execution code ---*/")
        code.mark_pos(None)

        if profile or linetrace:
            assert code.funcstate.gil_owned
            code.put_trace_start(modinit_func_name, self.pos)
            code.funcstate.can_trace = True

        code.mark_pos(None)
        self.body.generate_execution_code(code)
        code.mark_pos(None)

        if profile or linetrace:
            code.funcstate.can_trace = False
            assert code.funcstate.gil_owned
            code.put_trace_return("Py_None", pos=self.pos)
            code.put_trace_exit()

        code.putln()
        code.putln("/*--- Wrapped vars code ---*/")
        self.generate_wrapped_entries_code(env, code)
        code.putln()

        if Options.generate_cleanup_code:
            code.globalstate.use_utility_code(
                UtilityCode.load_cached("RegisterModuleCleanup", "ModuleSetupCode.c"))
            code.putln("if (__Pyx_RegisterCleanup()) %s" % code.error_goto(self.pos))

        code.put_goto(code.return_label)
        code.put_label(code.error_label)
        for cname, type in code.funcstate.all_managed_temps():
            code.put_xdecref(cname, type)

        if profile or linetrace:
            code.put_trace_exception_propagating()
            code.put_trace_unwind(self.pos)

        code.putln('if (%s) {' % env.module_cname)
        code.putln(
            f'if ({code.name_in_main_c_code_module_state(env.module_dict_cname)} && stringtab_initialized) {{')
        # We can run into errors before the module or stringtab are initialized.
        # In this case it is not safe to add a traceback (because it uses the stringtab)
        code.put_add_traceback(EncodedString("init %s" % env.qualified_name))
        code.globalstate.use_utility_code(Nodes.traceback_utility_code)
        # Module reference and module dict are in global variables which might still be needed
        # for cleanup, atexit code, etc., so leaking is better than crashing.
        # At least clearing the module dict here might be a good idea, but could still break
        # user code in atexit or other global registries.
        ##code.put_decref_clear(env.module_dict_cname, py_object_type, nanny=False)
        code.putln('}')
        code.putln("#if !CYTHON_USE_MODULE_STATE")
        code.put_decref_clear(env.module_cname, py_object_type, nanny=False, clear_before_decref=True)
        code.putln("#else")
        # This section is mainly for the limited API. env.module_cname still owns a reference so
        # decrement that
        code.put_decref(env.module_cname, py_object_type, nanny=False)
        # Also remove the failed module from the module state lookup
        # fetch/restore the error indicator because PyState_RemoveModule might fail itself
        code.putln("if (pystate_addmodule_run) {")
        code.putln("PyObject *tp, *value, *tb;")
        code.putln("PyErr_Fetch(&tp, &value, &tb);")
        code.putln("PyState_RemoveModule(&%s);" % Naming.pymoduledef_cname)
        code.putln("PyErr_Restore(tp, value, tb);")
        code.putln("}")
        code.putln("#endif")
        code.putln('} else if (!PyErr_Occurred()) {')
        code.putln('PyErr_SetString(PyExc_ImportError, "init %s");' %
                   env.qualified_name.as_c_string_literal()[1:-1])
        code.putln('}')
        code.put_label(code.return_label)

        code.put_finish_refcount_context()

        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        code.putln("return (%s != NULL) ? 0 : -1;" % env.module_cname)
        code.putln("#else")
        code.putln("return %s;" % env.module_cname)
        code.putln("#endif")
        code.putln('}')

        tempdecl_code.put_temp_declarations(code.funcstate)

        code.exit_cfunc_scope()

    def mod_init_subfunction(self, pos, scope, orig_code):
        """
        Return a context manager that allows deviating the module init code generation
        into a separate function and instead inserts a call to it.

        Can be reused sequentially to create multiple functions.
        The functions get inserted at the point where the context manager was created.
        The call gets inserted where the context manager is used (on entry).
        """
        function_code = orig_code.insertion_point()

        class ModInitSubfunction:
            def __init__(self, code_type):
                cname = '_'.join(code_type.lower().split())
                assert re.match("^[a-z0-9_]+$", cname)
                self.cfunc_name = "__Pyx_modinit_%s" % cname
                self.description = code_type
                self.tempdecl_code = None
                self.call_code = None

            def set_call_code(self, code):
                self.call_code = code.insertion_point()

            def __enter__(self):
                if self.call_code is None:
                    self.call_code = orig_code.insertion_point()
                code = function_code
                code.start_initcfunc(
                    f"int {self.cfunc_name}({Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname})",
                    scope, refnanny=True)
                code.putln(f"CYTHON_UNUSED_VAR({Naming.modulestatevalue_cname});")
                self.tempdecl_code = code.insertion_point()
                code.put_setup_refcount_context(EncodedString(self.cfunc_name))
                # Leave a grepable marker that makes it easy to find the generator source.
                code.putln("/*--- %s ---*/" % self.description)
                return code

            def __exit__(self, exc_type, exc_value, exc_tb):
                if exc_type is not None:
                    # Don't generate any code or do any validations on errors.
                    self.tempdecl_code = self.call_code = None
                    return

                code = function_code
                code.put_finish_refcount_context()
                code.putln("return 0;")

                self.tempdecl_code.put_temp_declarations(code.funcstate)
                self.tempdecl_code = None

                needs_error_handling = code.label_used(code.error_label)
                if needs_error_handling:
                    code.put_label(code.error_label)
                    for cname, type in code.funcstate.all_managed_temps():
                        code.put_xdecref(cname, type)
                    code.put_finish_refcount_context()
                    code.putln("return -1;")
                code.putln("}")
                code.exit_cfunc_scope()

                if needs_error_handling:
                    self.call_code.putln(
                        self.call_code.error_goto_if_neg("%s(%s)" % (
                            self.cfunc_name, Naming.modulestatevalue_cname), pos))
                else:
                    self.call_code.putln(
                        f"(void){self.cfunc_name}({Naming.modulestatevalue_cname});")
                self.call_code = None

        return ModInitSubfunction

    def generate_module_import_setup(self, env, code):
        module_path = env.directives['set_initial_path']
        if module_path == 'SOURCEFILE':
            module_path = self.pos[0].filename

        if module_path:
            code.putln('if (!CYTHON_PEP489_MULTI_PHASE_INIT) {')
            code.putln('if (PyObject_SetAttrString(%s, "__file__", %s) < 0) %s;' % (
                env.module_cname,
                code.get_py_string_const(
                    EncodedString(decode_filename(module_path))),
                code.error_goto(self.pos)))
            code.putln("}")

            if env.is_package:
                # set __path__ to mark the module as package
                code.putln('if (!CYTHON_PEP489_MULTI_PHASE_INIT) {')
                temp = code.funcstate.allocate_temp(py_object_type, True)
                code.putln('%s = Py_BuildValue("[O]", %s); %s' % (
                    temp,
                    code.get_py_string_const(
                        EncodedString(decode_filename(
                            os.path.dirname(module_path)))),
                    code.error_goto_if_null(temp, self.pos)))
                code.put_gotref(temp, py_object_type)
                code.putln(
                    'if (PyObject_SetAttrString(%s, "__path__", %s) < 0) %s;' % (
                        env.module_cname, temp, code.error_goto(self.pos)))
                code.put_decref_clear(temp, py_object_type)
                code.funcstate.release_temp(temp)
                code.putln("}")

        elif env.is_package:
            # packages require __path__, so all we can do is try to figure
            # out the module path at runtime by rerunning the import lookup
            code.putln("if (!CYTHON_PEP489_MULTI_PHASE_INIT) {")
            code.globalstate.use_utility_code(UtilityCode.load(
                "SetPackagePathFromImportLib", "ImportExport.c"))
            code.putln(code.error_goto_if_neg(
                '__Pyx_SetPackagePathFromImportLib(%s)' % (
                    code.get_py_string_const(
                        EncodedString(self.full_module_name))),
                self.pos))
            code.putln("}")

        # CPython may not have put us into sys.modules yet, but relative imports and reimports require it
        fq_module_name = self.full_module_name
        if fq_module_name.endswith('.__init__'):
            fq_module_name = EncodedString(fq_module_name[:-len('.__init__')])
        fq_module_name_cstring = fq_module_name.as_c_string_literal()
        code.putln("{")
        code.putln("PyObject *modules = PyImport_GetModuleDict(); %s" %
                   code.error_goto_if_null("modules", self.pos))
        code.putln('if (!PyDict_GetItemString(modules, %s)) {' % fq_module_name_cstring)
        code.putln(code.error_goto_if_neg('PyDict_SetItemString(modules, %s, %s)' % (
            fq_module_name_cstring, env.module_cname), self.pos))
        code.putln("}")
        code.putln("}")

    def generate_module_cleanup_func(self, env, code):
        if not Options.generate_cleanup_code:
            return

        code.putln('static void %s(CYTHON_UNUSED PyObject *self) {' %
                   Naming.cleanup_cname)
        code.enter_cfunc_scope(env)
        code.putln(f"{Naming.modulestatetype_cname} *{Naming.modulestatevalue_cname};")

        # TODO - this should go away when module-state has been refactored more and
        # we are able to access the module state through "self". Currently the
        # `#define` for each entry forces us to access it through PyState_FindModule
        # which is sometime unreliable during destruction
        # (e.g. during interpreter shutdown).
        # In that case the safest thing is to give up.
        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln(f"if (!__Pyx_State_FindModule(&{Naming.pymoduledef_cname})) return;")
        code.putln("#endif")
        code.putln(f"{Naming.modulestatevalue_cname} = __Pyx_PyModule_GetState(self);")

        if Options.generate_cleanup_code >= 2:
            code.putln("/*--- Global cleanup code ---*/")
            rev_entries = list(env.var_entries)
            rev_entries.reverse()
            for entry in rev_entries:
                if entry.visibility != 'extern':
                    if entry.type.is_pyobject and entry.used:
                        if entry.is_cglobal:
                            # TODO - eventually these should probably be in the module state too
                            entry_cname = entry.cname
                        else:
                            entry_cname = code.name_in_module_state(entry.cname)
                        code.put_xdecref_clear(
                            entry_cname, entry.type,
                            clear_before_decref=True,
                            nanny=False)
                    if entry.type.needs_explicit_destruction(env):
                        entry.type.generate_explicit_destruction(code, entry)
        code.putln(f"__Pyx_CleanupGlobals({Naming.modulestatevalue_cname});")
        if Options.generate_cleanup_code >= 3:
            code.putln("/*--- Type import cleanup code ---*/")
            for ext_type in sorted(env.types_imported, key=operator.attrgetter('typeptr_cname')):
                typeptr_cname = code.name_in_main_c_code_module_state(ext_type.typeptr_cname)
                code.put_xdecref_clear(
                    typeptr_cname, ext_type,
                    clear_before_decref=True,
                    nanny=False)
        if Options.cache_builtins:
            code.putln("/*--- Builtin cleanup code ---*/")
            for entry in env.cached_builtins:
                code.put_xdecref_clear(
                    entry.cname, PyrexTypes.py_object_type,
                    clear_before_decref=True,
                    nanny=False)
        code.putln("/*--- Intern cleanup code ---*/")
        code.put_decref_clear(f"{code.name_in_main_c_code_module_state(Naming.empty_tuple)}",
                              PyrexTypes.py_object_type,
                              clear_before_decref=True,
                              nanny=False)
        for entry in env.c_class_entries:
            cclass_type = entry.type
            if cclass_type.is_external or cclass_type.base_type:
                continue
            if cclass_type.scope.directives.get('freelist', 0):
                scope = cclass_type.scope
                freelist_name = code.name_in_main_c_code_module_state(
                    scope.mangle_internal(Naming.freelist_name))
                freecount_name = code.name_in_main_c_code_module_state(
                    scope.mangle_internal(Naming.freecount_name))
                code.putln('#if CYTHON_USE_FREELISTS')
                code.putln("while (%s > 0) {" % freecount_name)
                code.putln("PyObject* o = (PyObject*)%s[--%s];" % (
                    freelist_name, freecount_name))
                code.putln("PyTypeObject *tp = Py_TYPE(o);")
                code.putln("#if CYTHON_USE_TYPE_SLOTS")
                code.putln("(*tp->tp_free)(o);")
                code.putln("#else")
                # Asking for PyType_GetSlot(..., Py_tp_free) seems to cause an error in pypy
                code.putln("freefunc tp_free = (freefunc)PyType_GetSlot(tp, Py_tp_free);")
                code.putln("if (tp_free) tp_free(o);")
                code.putln("#endif")
                code.putln("#if CYTHON_USE_TYPE_SPECS")
                # Release the reference that "o" owned for its type.
                code.putln("Py_DECREF(tp);")
                code.putln("#endif")
                code.putln("}")
                code.putln('#endif')  # CYTHON_USE_FREELISTS
#        for entry in env.pynum_entries:
#            code.put_decref_clear(entry.cname,
#                                  PyrexTypes.py_object_type,
#                                  nanny=False)
#        for entry in env.all_pystring_entries:
#            if entry.is_interned:
#                code.put_decref_clear(entry.pystring_cname,
#                                      PyrexTypes.py_object_type,
#                                      nanny=False)
#        for entry in env.default_entries:
#            if entry.type.is_pyobject and entry.used:
#                code.putln("Py_DECREF(%s); %s = 0;" % (
#                    code.entry_as_pyobject(entry), entry.cname))
        if Options.pre_import is not None:
            code.put_decref_clear(Naming.preimport_cname, py_object_type,
                                  nanny=False, clear_before_decref=True)
        for cname in [Naming.cython_runtime_cname, Naming.builtins_cname]:
            cname = code.name_in_main_c_code_module_state(cname)
            code.put_decref_clear(cname, py_object_type, nanny=False, clear_before_decref=True)
        code.put_decref_clear(
            code.name_in_main_c_code_module_state(env.module_dict_cname),
            py_object_type, nanny=False, clear_before_decref=True)

    def generate_main_method(self, env, code):
        module_is_main = self.is_main_module_flag_cname()
        if Options.embed == "main":
            wmain = "wmain"
        else:
            wmain = Options.embed
        main_method = TempitaUtilityCode.load_cached(
                "MainFunction", "Embed.c",
                context={
                    'module_name': env.module_name,
                    'module_is_main': module_is_main,
                    'main_method': Options.embed,
                    'wmain_method': wmain,
                    'embed_modules': tuple(Options.embed_modules)})
        code.globalstate.use_utility_code(main_method)

    def punycode_module_name(self, prefix, name):
        # adapted from PEP483
        if name.isascii():
            name = '_' + name
        else:
            name = 'U_' + name.encode('punycode').replace(b'-', b'_').decode('ascii')
        return "%s%s" % (prefix, name)

    def wrong_punycode_module_name(self, name):
        # to work around a distutils bug by also generating an incorrect symbol...
        if name.isascii():
            return None  # workaround is not needed
        return "PyInitU" + ("_"+name).encode('punycode').replace(b'-', b'_').decode('ascii')

    def mod_init_func_cname(self, prefix, env):
        # from PEP483
        return self.punycode_module_name(prefix, env.module_name)

    # Returns the name of the C-function that corresponds to the module initialisation.
    # (module initialisation == the cython code outside of functions)
    # Note that this should never be the name of a wrapper and always the name of the
    # function containing the actual code. Otherwise, cygdb will experience problems.
    def module_init_func_cname(self):
        env = self.scope
        return self.mod_init_func_cname(Naming.pymodule_exec_func_cname, env)

    def generate_pymoduledef_struct(self, env, code):
        if env.doc:
            doc = "%s" % code.get_string_const(env.doc)
        else:
            doc = "0"
        if Options.generate_cleanup_code:
            cleanup_func = "(freefunc)%s" % Naming.cleanup_cname
        else:
            cleanup_func = 'NULL'

        code.putln("")
        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        exec_func_cname = self.module_init_func_cname()
        code.putln("static PyObject* %s(PyObject *spec, PyModuleDef *def); /*proto*/" %
                   Naming.pymodule_create_func_cname)
        code.putln("static int %s(PyObject* module); /*proto*/" % exec_func_cname)

        code.putln("static PyModuleDef_Slot %s[] = {" % Naming.pymoduledef_slots_cname)
        code.putln("{Py_mod_create, (void*)%s}," % Naming.pymodule_create_func_cname)
        code.putln("{Py_mod_exec, (void*)%s}," % exec_func_cname)
        code.putln("#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING")
        code.putln("{Py_mod_gil, __Pyx_FREETHREADING_COMPATIBLE},")
        code.putln("#endif")
        code.putln("#if PY_VERSION_HEX >= 0x030C0000 && CYTHON_USE_MODULE_STATE")
        subinterp_option = {
            'no': 'Py_MOD_MULTIPLE_INTERPRETERS_NOT_SUPPORTED',
            'shared_gil': 'Py_MOD_MULTIPLE_INTERPRETERS_SUPPORTED',
            'own_gil': 'Py_MOD_PER_INTERPRETER_GIL_SUPPORTED'
        }.get(env.directives["subinterpreters_compatible"])
        code.putln("{Py_mod_multiple_interpreters, %s}," % subinterp_option)
        code.putln("#endif")
        code.putln("{0, NULL}")
        code.putln("};")
        if not env.module_name.isascii():
            code.putln("#else /* CYTHON_PEP489_MULTI_PHASE_INIT */")
            code.putln('#error "Unicode module names are only supported with multi-phase init'
                       ' as per PEP489"')
        code.putln("#endif")

        code.putln("")
        code.putln('#ifdef __cplusplus')
        code.putln('namespace {')
        code.putln("struct PyModuleDef %s =" % Naming.pymoduledef_cname)
        code.putln('#else')
        code.putln("static struct PyModuleDef %s =" % Naming.pymoduledef_cname)
        code.putln('#endif')
        code.putln('{')
        code.putln("  PyModuleDef_HEAD_INIT,")
        code.putln('  %s,' % env.module_name.as_c_string_literal())
        code.putln("  %s, /* m_doc */" % doc)
        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln(f"  sizeof({Naming.modulestatetype_cname}), /* m_size */")
        code.putln("#else")
        code.putln("  (CYTHON_PEP489_MULTI_PHASE_INIT) ? 0 : -1, /* m_size */")
        code.putln("#endif")
        code.putln("  %s /* m_methods */," % env.method_table_cname)
        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        code.putln("  %s, /* m_slots */" % Naming.pymoduledef_slots_cname)
        code.putln("#else")
        code.putln("  NULL, /* m_reload */")
        code.putln("#endif")
        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln("  %s_traverse, /* m_traverse */" % Naming.module_cname)
        code.putln("  %s_clear, /* m_clear */" % Naming.module_cname)
        code.putln("  %s /* m_free */" % cleanup_func)
        code.putln("#else")
        code.putln("  NULL, /* m_traverse */")
        code.putln("  NULL, /* m_clear */")
        code.putln("  %s /* m_free */" % cleanup_func)
        code.putln("#endif")
        code.putln("};")
        code.putln('#ifdef __cplusplus')
        code.putln('} /* anonymous namespace */')
        code.putln('#endif')

    def generate_module_creation_code(self, env, code):
        # Generate code to create the module object and
        # install the builtins.
        if env.doc:
            doc = "%s" % code.get_string_const(env.doc)
        else:
            doc = "0"

        # manage_ref is False (and refnanny calls are omitted) because refnanny isn't yet initialized.
        module_temp = code.funcstate.allocate_temp(py_object_type, manage_ref=False)
        code.putln("#if CYTHON_PEP489_MULTI_PHASE_INIT")
        code.putln("%s = %s;" % (
            module_temp,
            Naming.pymodinit_module_arg))
        code.put_incref(module_temp, py_object_type, nanny=False)
        code.putln("#else")
        code.putln(
            "%s = PyModule_Create(&%s); %s" % (
                module_temp,
                Naming.pymoduledef_cname,
                code.error_goto_if_null(module_temp, self.pos)))
        code.putln("#endif")

        code.putln("#if CYTHON_USE_MODULE_STATE")
        code.putln("{")
        # So that PyState_FindModule works in the init function:
        code.putln("int add_module_result = __Pyx_State_AddModule(%s, &%s);" % (
            module_temp, Naming.pymoduledef_cname))
        code.putln("%s = 0; /* transfer ownership from %s to %s pseudovariable */" % (
            module_temp, module_temp, env.module_name.as_c_string_literal()
        ))
        # At this stage the module likely has a refcount of 2 - one owned by the list
        # inside PyState_AddModule and one owned by "__pyx_m" (and returned from this
        # function as a new reference).
        code.putln(code.error_goto_if_neg("add_module_result", self.pos))
        code.putln("pystate_addmodule_run = 1;")
        code.putln("}")
        code.putln('#else')  # !CYTHON_USE_MODULE_STATE
        code.putln(f"{env.module_cname} = {module_temp};")
        code.putln("#endif")
        code.funcstate.release_temp(module_temp)

        code.putln("#if CYTHON_COMPILING_IN_CPYTHON_FREETHREADING")
        gil_option = ("Py_MOD_GIL_NOT_USED"
                      if env.directives["freethreading_compatible"]
                      else "Py_MOD_GIL_USED")
        code.putln(f"PyUnstable_Module_SetGIL({env.module_cname}, {gil_option});")
        code.putln("#endif")

        code.putln(f"{Naming.modulestatevalue_cname} = {Naming.modulestateglobal_cname};")
        code.putln("CYTHON_UNUSED_VAR(%s);" % module_temp)  # only used in limited API

        dict_cname = code.name_in_main_c_code_module_state(env.module_dict_cname)
        code.putln(
            "%s = PyModule_GetDict(%s); %s" % (
                dict_cname, env.module_cname,
                code.error_goto_if_null(dict_cname, self.pos)))
        code.put_incref(dict_cname, py_object_type, nanny=False)

        builtins_cname = code.name_in_main_c_code_module_state(Naming.builtins_cname)
        code.putln(
            '%s = __Pyx_PyImport_AddModuleRef(__Pyx_BUILTIN_MODULE_NAME); %s' % (
                builtins_cname,
                code.error_goto_if_null(builtins_cname, self.pos)))
        runtime_cname = code.name_in_main_c_code_module_state(Naming.cython_runtime_cname)
        code.putln(
            '%s = __Pyx_PyImport_AddModuleRef("cython_runtime"); %s' % (
                runtime_cname,
                code.error_goto_if_null(runtime_cname, self.pos)))
        code.putln(
            'if (PyObject_SetAttrString(%s, "__builtins__", %s) < 0) %s' % (
                env.module_cname,
                builtins_cname,
                code.error_goto(self.pos)))
        if Options.pre_import is not None:
            code.putln(
                '%s = __Pyx_PyImport_AddModuleRef("%s"); %s' % (
                    Naming.preimport_cname,
                    Options.pre_import,
                    code.error_goto_if_null(Naming.preimport_cname, self.pos)))

    def generate_global_init_code(self, env, code):
        # Generate code to initialise global PyObject *
        # variables to None.
        for entry in env.var_entries:
            if entry.visibility != 'extern':
                if entry.used:
                    entry.type.global_init_code(entry, code)
                if entry.type.needs_explicit_construction(env):
                    # TODO - this is slightly redundant with global_init_code
                    entry.type.generate_explicit_construction(code, entry)

    def generate_wrapped_entries_code(self, env, code):
        for name, entry in sorted(env.entries.items()):
            if (entry.create_wrapper
                    and not entry.is_type
                    and entry.scope is env):
                if not entry.type.create_to_py_utility_code(env):
                    error(entry.pos, "Cannot convert '%s' to Python object" % entry.type)
                code.putln("{")
                code.putln("PyObject* wrapped = %s(%s);"  % (
                    entry.type.to_py_function,
                    entry.cname))
                code.putln(code.error_goto_if_null("wrapped", entry.pos))
                code.putln(
                    'if (PyObject_SetAttrString(%s, "%s", wrapped) < 0) %s;' % (
                        env.module_cname,
                        name,
                        code.error_goto(entry.pos)))
                code.putln("}")

    def _select_exported_entries(self, all_entries):
        return [
            entry for entry in all_entries
            if entry.api or entry.defined_in_pxd or (Options.cimport_from_pyx and entry.visibility != 'extern')
        ]

    def generate_c_variable_export_code(self, env, code):
        """Generate code to create PyCFunction wrappers for exported C functions.
        """
        entries = self._select_exported_entries(env.var_entries)
        if not entries:
            return

        exports = [
            # (signature, name, cname)
            (entry.type.empty_declaration_code(), entry.name, entry.cname)
            for entry in entries
        ]
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("VoidPtrExport", "ImportExport.c"))

        _generate_export_code(code, self.pos, exports, "__Pyx_ExportVoidPtr", "void *{name}")

    def generate_c_function_export_code(self, env, code):
        """Generate code to create PyCFunction wrappers for exported C functions.
        """
        entries = self._select_exported_entries(env.cfunc_entries)
        if not entries:
            return

        exports = [
            # (signature, name, cname)
            (entry.type.signature_string(), entry.name, entry.cname)
            for entry in entries
        ]
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("FunctionExport", "ImportExport.c"))

        _generate_export_code(code, self.pos, exports, "__Pyx_ExportFunction", "void (*{name})(void)")

    def generate_type_import_code_for_module(self, module, env, code):
        # Generate type import code for all exported extension types in
        # an imported module.
        #if module.c_class_entries:
        with ModuleImportGenerator(code) as import_generator:
            for entry in module.c_class_entries:
                if entry.defined_in_pxd:
                    self.generate_type_import_code(env, entry.type, entry.pos, code, import_generator)

    def specialize_fused_types(self, pxd_env):
        """
        If fused c(p)def functions are defined in an imported pxd, but not
        used in this implementation file, we still have fused entries and
        not specialized ones. This method replaces any fused entries with their
        specialized ones.
        """
        for entry in pxd_env.cfunc_entries[:]:
            if entry.type.is_fused:
                # This call modifies the cfunc_entries in-place
                entry.type.get_all_specialized_function_types()

    def _select_imported_entries(self, all_entries, used_only=False):
        return [
            entry for entry in all_entries
            if entry.defined_in_pxd and (not used_only or entry.used)
        ]

    def generate_c_variable_import_code_for_module(self, module, env, code):
        """Generate import code for all exported C functions in a cimported module.
        """
        entries = self._select_imported_entries(module.var_entries, used_only=False)
        if not entries:
            return

        is_module_scope = env is module
        imports = [
            # (signature, name, cname)
            (entry.type.empty_declaration_code(),
             entry.name,
             entry.cname if is_module_scope else module.mangle(Naming.varptr_prefix, entry.name))
            for entry in entries
        ]
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("VoidPtrImport", "ImportExport.c"))

        _generate_import_code(
            code, self.pos, imports, module.qualified_name, f"__Pyx_ImportVoidPtr_{Naming.cyversion}", "void **{name}")

    def generate_c_function_import_code_for_module(self, module, env, code):
        """Generate import code for all exported C functions in a cimported module.
        """
        entries = self._select_imported_entries(module.cfunc_entries, used_only=True)
        if not entries:
            return

        imports = [
            # (signature, name, cname)
            (entry.type.signature_string(), entry.name, entry.cname)
            for entry in entries
        ]
        code.globalstate.use_utility_code(
            UtilityCode.load_cached("FunctionImport", "ImportExport.c"))

        _generate_import_code(
            code, self.pos, imports, module.qualified_name, f"__Pyx_ImportFunction_{Naming.cyversion}", "void (**{name})(void)")

    def generate_type_init_code(self, env, code):
        # Generate type import code for extern extension types
        # and type ready code for non-extern ones.
        with ModuleImportGenerator(code) as import_generator:
            for entry in env.c_class_entries:
                if entry.visibility == 'extern' and not entry.utility_code_definition:
                    self.generate_type_import_code(env, entry.type, entry.pos, code, import_generator)
                else:
                    self.generate_base_type_import_code(env, entry, code, import_generator)
                    self.generate_exttype_vtable_init_code(entry, code)
                    if entry.type.early_init:
                        self.generate_type_ready_code(entry, code)

    def generate_base_type_import_code(self, env, entry, code, import_generator):
        base_type = entry.type.base_type
        if (base_type and base_type.module_name != env.qualified_name and not
                (base_type.is_builtin_type or base_type.is_cython_builtin_type)
                 and not entry.utility_code_definition):
            self.generate_type_import_code(env, base_type, self.pos, code, import_generator)

    def generate_type_import_code(self, env, type, pos, code, import_generator):
        # If not already done, generate code to import the typeobject of an
        # extension type defined in another module, and extract its C method
        # table pointer if any.
        if type in env.types_imported:
            return
        if type.name not in Code.ctypedef_builtins_map:
            # see corresponding condition in generate_type_import_call() below!
            code.globalstate.use_utility_code(
                UtilityCode.load_cached("TypeImport", "ImportExport.c"))
        self.generate_type_import_call(type, code, import_generator, error_pos=pos)
        if type.vtabptr_cname:
            code.globalstate.use_utility_code(
                UtilityCode.load_cached('GetVTable', 'ImportExport.c'))
            code.putln("%s = (struct %s*)__Pyx_GetVtable(%s); %s" % (
                type.vtabptr_cname,
                type.vtabstruct_cname,
                code.name_in_main_c_code_module_state(type.typeptr_cname),
                code.error_goto_if_null(type.vtabptr_cname, pos)))
        env.types_imported.add(type)

    def generate_type_import_call(self, type, code, import_generator, error_code=None, error_pos=None, is_api=False):
        sizeof_objstruct = objstruct = type.objstruct_cname if type.typedef_flag else f"struct {type.objstruct_cname}"
        module_name = type.module_name
        type_name = type.name
        is_builtin = module_name in ('__builtin__', 'builtins')
        if not is_builtin:
            module_name = f'"{module_name}"'
        elif type_name in Code.ctypedef_builtins_map:
            # Fast path for special builtins, don't actually import
            code.putln(
                f'{code.name_in_module_state(type.typeptr_cname)} = {Code.ctypedef_builtins_map[type_name]};')
            return
        else:
            module_name = '__Pyx_BUILTIN_MODULE_NAME'
            if type_name in Code.renamed_py2_builtins_map:
                type_name = Code.renamed_py2_builtins_map[type_name]
            if objstruct in Code.basicsize_builtins_map:
                # Some builtin types have a tp_basicsize which differs from sizeof(...):
                sizeof_objstruct = Code.basicsize_builtins_map[objstruct]

        if not error_code:
            assert error_pos is not None
            error_code = code.error_goto(error_pos)

        module = import_generator.imported_module(module_name, error_code)
        typeptr_cname = type.typeptr_cname
        if not is_api:
            typeptr_cname = code.name_in_main_c_code_module_state(typeptr_cname)

        code.putln(
            f"{typeptr_cname} = __Pyx_ImportType_{Naming.cyversion}("
            f"{module}, {module_name}, {type.name.as_c_string_literal()},"
        )

        alignment_func = f"__PYX_GET_STRUCT_ALIGNMENT_{Naming.cyversion}"
        code.putln("#if defined(PYPY_VERSION_NUM) && PYPY_VERSION_NUM < 0x050B0000")
        code.putln(f'sizeof({objstruct}), {alignment_func}({objstruct}),')
        code.putln("#elif CYTHON_COMPILING_IN_LIMITED_API")
        if is_builtin:
            # Builtin types are opaque in when the limited API is enabled
            # and subsequents attempt to access their fields will trigger
            # compile errors. Skip the struct size check here so things keep
            # working when a builtin type is imported but not actually used.
            code.putln('0, 0,')
        else:
            code.putln(f'sizeof({objstruct}), {alignment_func}({objstruct}),')
        code.putln('#else')
        code.putln(f'sizeof({sizeof_objstruct}), {alignment_func}({sizeof_objstruct}),')
        code.putln("#endif")

        # check_size
        if type.check_size and type.check_size in ('error', 'warn', 'ignore'):
            check_size = type.check_size
        elif not type.is_external or type.is_subclassed:
            check_size = 'error'
        else:
            raise RuntimeError(
                f"invalid value for check_size '{type.check_size}' when compiling {module_name}.{type.name}")
        code.put(f'__Pyx_ImportType_CheckSize_{check_size.title()}_{Naming.cyversion});')

        code.putln(f' if (!{typeptr_cname}) {error_code}')
    def generate_type_ready_code(self, entry, code):
        Nodes.CClassDefNode.generate_type_ready_code(entry, code)

    def is_main_module_flag_cname(self):
        full_module_name = self.full_module_name.replace('.', '__')
        return self.punycode_module_name(Naming.module_is_main, full_module_name)

    def generate_exttype_vtable_init_code(self, entry, code):
        # Generate code to initialise the C method table of an
        # extension type.
        type = entry.type
        if type.vtable_cname:
            code.putln(
                "%s = &%s;" % (
                    type.vtabptr_cname,
                    type.vtable_cname))
            if type.base_type and type.base_type.vtabptr_cname:
                code.putln(
                    "%s.%s = *%s;" % (
                        type.vtable_cname,
                        Naming.obj_base_cname,
                        type.base_type.vtabptr_cname))

            c_method_entries = [
                entry for entry in type.scope.cfunc_entries
                if entry.func_cname]
            if c_method_entries:
                for meth_entry in c_method_entries:
                    vtable_type = meth_entry.vtable_type or meth_entry.type
                    cast = vtable_type.signature_cast_string()
                    code.putln(
                        "%s.%s = %s%s;" % (
                            type.vtable_cname,
                            meth_entry.cname,
                            cast,
                            meth_entry.func_cname))


# cimport/export code for functions and pointers.

def _deduplicate_inout_signatures(item_tuples):
    # We can save runtime space for identical signatures by reusing the same C strings.
    # To deduplicate the signatures, we sort by them and store duplicates as empty C strings.
    signatures, names, items = zip(*sorted(item_tuples))
    signatures = list(signatures)  # tuple -> list, to allow reassignments again

    last_sig = None
    for i, signature in enumerate(signatures):
        if signature == last_sig:
            signatures[i] = ''
        else:
            last_sig = signature

    return signatures, names, items


def _generate_import_export_code(code: Code.CCodeWriter, pos, inout_item_tuples, per_item_func, target, pointer_decl, use_pybytes, is_import):
    signatures, names, inout_items = _deduplicate_inout_signatures(inout_item_tuples)

    pyx = f"{Naming.pyrex_prefix}{'import' if is_import else 'export'}_"

    pointer_cast = pointer_decl.format(name='')
    sig_bytes = '\0'.join(signatures).encode('utf-8')
    names_bytes = '\0'.join(names).encode('utf-8')
    pointers = [f"({pointer_cast})&{item_cname}" for item_cname in inout_items]

    sigs_and_names_bytes = bytes_literal(sig_bytes + b'\0' + names_bytes, 'utf-8')
    if use_pybytes:
        code.putln(f"const char * {pyx}signature = __Pyx_PyBytes_AsString({code.get_py_string_const(sigs_and_names_bytes)});")
        code.putln("#if !CYTHON_ASSUME_SAFE_MACROS")
        code.putln(code.error_goto_if_null(f'{pyx}signature', pos))
        code.putln("#endif")
    else:
        code.putln(f"const char * {pyx}signature = {sigs_and_names_bytes.as_c_string_literal()};")

    code.putln(f"const char * {pyx}name = {pyx}signature + {len(sig_bytes) + 1};")
    code.putln(f"{pointer_decl.format(name=f'const {pyx}pointers[]')} = {{{', '.join(pointers)}, ({pointer_cast}) NULL}};")

    code.globalstate.use_utility_code(
        UtilityCode.load_cached("IncludeStringH", "StringTools.c"))

    code.putln(f"{pointer_decl.format(name=f'const *{pyx}pointer')} = {pyx}pointers;")
    code.putln(f"const char *{pyx}current_signature = {pyx}signature;")
    code.putln(f"while (*{pyx}pointer) {{")

    code.put_error_if_neg(
        pos,
        f"{per_item_func}({target}, {pyx}name, *{pyx}pointer, {pyx}current_signature)"
    )
    code.putln(f"++{pyx}pointer;")
    code.putln(f"{pyx}name = strchr({pyx}name, '\\0') + 1;")
    code.putln(f"{pyx}signature = strchr({pyx}signature, '\\0') + 1;")
    # Keep reusing the current signature until we find a new non-empty one.
    code.putln(f"if (*{pyx}signature != '\\0') {pyx}current_signature = {pyx}signature;")

    code.putln("}")  # while


def _generate_export_code(code: Code.CCodeWriter, pos, exports, export_func, pointer_decl):
    """Generate function/pointer export code.

    'exports' is a list of (signature, name, exported_cname) tuples.
    """
    code.putln("{")

    api_dict = code.funcstate.allocate_temp(py_object_type, manage_ref=True)
    code.globalstate.use_utility_code(
        UtilityCode.load_cached("GetApiDict", "ImportExport.c"))
    code.putln(
        f"{api_dict} = __Pyx_ApiExport_GetApiDict(); "
        f"{code.error_goto_if_null(api_dict, pos)}"
    )
    code.put_gotref(api_dict, py_object_type)

    _generate_import_export_code(code, pos, exports, export_func, api_dict, pointer_decl, use_pybytes=True, is_import=False)

    code.put_decref_clear(api_dict, py_object_type)
    code.funcstate.release_temp(api_dict)

    code.putln("}")


def _generate_import_code(code, pos, imports, qualified_module_name, import_func, pointer_decl):
    """Generate function/pointer import code.

    'imports' is a list of (signature, name, imported_cname) tuples.
    """
    code.putln("{")

    module_ref = code.funcstate.allocate_temp(py_object_type, manage_ref=True)
    code.putln(
        f'{module_ref} = PyImport_ImportModule({qualified_module_name.as_c_string_literal()}); '
        f'{code.error_goto_if_null(module_ref, pos)}'
    )
    code.put_gotref(module_ref, py_object_type)

    _generate_import_export_code(code, pos, imports, import_func, module_ref, pointer_decl, use_pybytes=True, is_import=True)

    code.put_decref_clear(module_ref, py_object_type)
    code.funcstate.release_temp(module_ref)

    code.putln("}")


# Module import helper

class ModuleImportGenerator:
    """
    Helper to generate module import while importing external types.
    This is used to avoid excessive re-imports of external modules when multiple types are looked up.
    """
    def __init__(self, code, imported_modules=None):
        self.code = code
        self.imported = {}
        if imported_modules:
            for name, cname in imported_modules.items():
                self.imported['"%s"' % name] = cname
        self.temps = []  # remember original import order for freeing

    def imported_module(self, module_name_string, error_code):
        if module_name_string in self.imported:
            return self.imported[module_name_string]

        code = self.code
        temp = code.funcstate.allocate_temp(py_object_type, manage_ref=True)
        self.temps.append(temp)
        code.putln('%s = PyImport_ImportModule(%s); if (unlikely(!%s)) %s' % (
            temp, module_name_string, temp, error_code))
        code.put_gotref(temp, py_object_type)
        self.imported[module_name_string] = temp
        return temp

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        code = self.code
        for temp in self.temps:
            code.put_decref_clear(temp, py_object_type)
            code.funcstate.release_temp(temp)


def generate_cfunction_declaration(entry, env, code, definition):
    from_cy_utility = entry.used and entry.utility_code_definition
    if entry.used and entry.inline_func_in_pxd or (not entry.in_cinclude and (
            definition or entry.defined_in_pxd or entry.visibility == 'extern' or from_cy_utility)):
        if entry.visibility == 'extern':
            storage_class = Naming.extern_c_macro
            dll_linkage = "DL_IMPORT"
        elif entry.visibility == 'public':
            storage_class = Naming.extern_c_macro
            dll_linkage = None
        elif entry.visibility == 'private':
            storage_class = "static"
            dll_linkage = None
        else:
            storage_class = "static"
            dll_linkage = None
        type = entry.type

        if entry.defined_in_pxd and not definition:
            storage_class = "static"
            dll_linkage = None
            type = CPtrType(type)

        header = type.declaration_code(
            entry.cname, dll_linkage=dll_linkage)
        modifiers = code.build_function_modifiers(entry.func_modifiers)
        code.putln("%s %s%s; /*proto*/" % (
            storage_class,
            modifiers,
            header))
        code.globalstate.use_entry_utility_code(entry)

#------------------------------------------------------------------------------------
#
#  Runtime support code
#
#------------------------------------------------------------------------------------

refnanny_utility_code = UtilityCode.load("Refnanny", "ModuleSetupCode.c")

packed_struct_utility_code = UtilityCode(proto="""
#if defined(__GNUC__)
#define __Pyx_PACKED __attribute__((__packed__))
#else
#define __Pyx_PACKED
#endif
""", impl="", proto_block='utility_code_proto_before_types')
