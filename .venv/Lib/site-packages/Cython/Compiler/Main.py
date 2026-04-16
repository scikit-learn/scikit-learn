#
#   Cython Top Level
#


import os
import re
import sys
import io

if sys.version_info[:2] < (3, 8):
    sys.stderr.write("Sorry, Cython requires Python 3.8+, found %d.%d\n" % tuple(sys.version_info[:2]))
    sys.exit(1)

# Do not import Parsing here, import it when needed, because Parsing imports
# Nodes, which globally needs debug command line options initialized to set a
# conditional metaclass. These options are processed by CmdLine called from
# main() in this file.
# import Parsing
from . import Errors
from .StringEncoding import EncodedString
from .Scanning import PyrexScanner, FileSourceDescriptor
from .Errors import PyrexError, CompileError, error, warning
from .Symtab import ModuleScope
from .. import Utils
from . import Options
from .Options import CompilationOptions, default_options
from .CmdLine import parse_command_line
from .Lexicon import (unicode_start_ch_any, unicode_continuation_ch_any,
                      unicode_start_ch_range, unicode_continuation_ch_range)


def _make_range_re(chrs):
    out = []
    for i in range(0, len(chrs), 2):
        out.append("{}-{}".format(chrs[i], chrs[i+1]))
    return "".join(out)

# py2 version looked like r"[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)*$"
module_name_pattern = "[{0}{1}][{0}{2}{1}{3}]*".format(
    unicode_start_ch_any, _make_range_re(unicode_start_ch_range),
    unicode_continuation_ch_any,
    _make_range_re(unicode_continuation_ch_range))
module_name_pattern = re.compile("{0}(\\.{0})*$".format(module_name_pattern))


standard_include_path = os.path.abspath(
    os.path.join(os.path.dirname(os.path.dirname(__file__)), 'Includes'))


class Context:
    #  This class encapsulates the context needed for compiling
    #  one or more Cython implementation files along with their
    #  associated and imported declaration files. It includes
    #  the root of the module import namespace and the list
    #  of directories to search for include files.
    #
    #  modules               {string : ModuleScope}
    #  include_directories   [string]
    #  future_directives     [object]
    #  language_level        int     currently 2 or 3 for Python 2/3

    cython_scope = None
    language_level = None  # warn when not set but default to Py2

    def __init__(self, include_directories, compiler_directives, cpp=False,
                 language_level=None, options=None):
        # cython_scope is a hack, set to False by subclasses, in order to break
        # an infinite loop.
        # Better code organization would fix it.

        from . import Builtin, CythonScope
        self.modules = {"__builtin__" : Builtin.builtin_scope}
        self.cython_scope = CythonScope.create_cython_scope(self)
        self.modules["cython"] = self.cython_scope
        self.include_directories = include_directories
        self.future_directives = set()
        self.compiler_directives = compiler_directives
        self.cpp = cpp
        self.options = options

        self.pxds = {}  # full name -> node tree
        self.utility_pxds = {}  # pxd name -> node tree
        self._interned = {}  # (type(value), value, *key_args) -> interned_value

        if language_level is not None:
            self.set_language_level(language_level)

        self.legacy_implicit_noexcept = self.compiler_directives.get('legacy_implicit_noexcept', False)

        self.gdb_debug_outputwriter = None

    @classmethod
    def from_options(cls, options):
        return cls(options.include_path, options.compiler_directives,
                   options.cplus, options.language_level, options=options)

    @property
    def shared_c_file_path(self):
        return self.options.shared_c_file_path if self.options else None

    @property
    def shared_utility_qualified_name(self):
        return self.options.shared_utility_qualified_name if self.options else None

    def set_language_level(self, level):
        from .Future import print_function, unicode_literals, absolute_import, division, generator_stop
        future_directives = set()
        if level == '3str':
            level = 3
        else:
            level = int(level)
        if level >= 3:
            future_directives.update([unicode_literals, print_function, absolute_import, division, generator_stop])
        self.language_level = level
        self.future_directives = future_directives
        if level >= 3:
            self.modules['builtins'] = self.modules['__builtin__']

    def intern_ustring(self, value, encoding=None):
        key = (EncodedString, value, encoding)
        try:
            return self._interned[key]
        except KeyError:
            pass
        value = EncodedString(value)
        if encoding:
            value.encoding = encoding
        self._interned[key] = value
        return value

    # pipeline creation functions can now be found in Pipeline.py

    def process_pxd(self, source_desc, scope, module_name):
        from . import Pipeline
        if isinstance(source_desc, FileSourceDescriptor) and source_desc._file_type == 'pyx':
            source = CompilationSource(source_desc, module_name, os.getcwd())
            result_sink = create_default_resultobj(source, self.options)
            pipeline = Pipeline.create_pyx_as_pxd_pipeline(self, result_sink)
            result = Pipeline.run_pipeline(pipeline, source)
        elif source_desc.in_utility_code:
            from . import ParseTreeTransforms
            transform = ParseTreeTransforms.CnameDirectivesTransform(self)
            pipeline = Pipeline.create_pxd_pipeline(self, scope, module_name)
            pipeline = Pipeline.insert_into_pipeline(
                pipeline, transform,
                before=ParseTreeTransforms.InterpretCompilerDirectives)
            result = Pipeline.run_pipeline(pipeline, source_desc)
        else:
            pipeline = Pipeline.create_pxd_pipeline(self, scope, module_name)
            result = Pipeline.run_pipeline(pipeline, source_desc)
        return result

    def nonfatal_error(self, exc):
        return Errors.report_error(exc)

    def _split_qualified_name(self, qualified_name, relative_import=False):
        # Splits qualified_name into parts in form of 2-tuples: (PART_NAME, IS_PACKAGE).
        qualified_name_parts = qualified_name.split('.')
        last_part = qualified_name_parts.pop()
        qualified_name_parts = [(p, True) for p in qualified_name_parts]
        if last_part != '__init__':
            # If Last part is __init__, then it is omitted. Otherwise, we need to check whether we can find
            # __init__.pyx/__init__.py file to determine if last part is package or not.
            is_package = False
            for suffix in ('.py', '.pyx'):
                path = self.search_include_directories(
                    qualified_name, suffix=suffix, source_pos=None, source_file_path=None, sys_path=not relative_import)
                if path:
                    is_package = self._is_init_file(path)
                    break

            qualified_name_parts.append((last_part, is_package))
        return qualified_name_parts

    @staticmethod
    def _is_init_file(path):
        return os.path.basename(path) in ('__init__.pyx', '__init__.py', '__init__.pxd') if path else False

    @staticmethod
    def _check_pxd_filename(pos, pxd_pathname, qualified_name):
        if not pxd_pathname:
            return
        pxd_filename = os.path.basename(pxd_pathname)
        if '.' in qualified_name and qualified_name == os.path.splitext(pxd_filename)[0]:
            warning(pos, "Dotted filenames ('%s') are deprecated."
                    " Please use the normal Python package directory layout." % pxd_filename, level=1)

    def find_module(self, module_name, from_module=None, pos=None, need_pxd=1,
                    absolute_fallback=True, relative_import=False):
        # Finds and returns the module scope corresponding to
        # the given relative or absolute module name. If this
        # is the first time the module has been requested, finds
        # the corresponding .pxd file and process it.
        # If from_module is not None, it must be a module scope,
        # and the module will first be searched for relative to
        # that module, provided its name is not a dotted name.
        debug_find_module = 0
        if debug_find_module:
            print("Context.find_module: module_name = %s, from_module = %s, pos = %s, need_pxd = %s" % (
                module_name, from_module, pos, need_pxd))

        scope = None
        pxd_pathname = None
        if from_module:
            if module_name:
                # from .module import ...
                qualified_name = from_module.qualify_name(module_name)
            else:
                # from . import ...
                qualified_name = from_module.qualified_name
                scope = from_module
                from_module = None
        else:
            qualified_name = module_name

        if not module_name_pattern.match(qualified_name):
            raise CompileError(pos or (module_name, 0, 0),
                               "'%s' is not a valid module name" % module_name)

        if from_module:
            if debug_find_module:
                print("...trying relative import")
            scope = from_module.lookup_submodule(module_name)
            if not scope:
                pxd_pathname = self.find_pxd_file(qualified_name, pos, sys_path=not relative_import)
                self._check_pxd_filename(pos, pxd_pathname, qualified_name)
                if pxd_pathname:
                    is_package = self._is_init_file(pxd_pathname)
                    scope = from_module.find_submodule(module_name, as_package=is_package)
        if not scope:
            if debug_find_module:
                print("...trying absolute import")
            if absolute_fallback:
                qualified_name = module_name
            scope = self
            for name, is_package in self._split_qualified_name(qualified_name, relative_import=relative_import):
                scope = scope.find_submodule(name, as_package=is_package)
        if debug_find_module:
            print("...scope = %s" % scope)
        if not scope.pxd_file_loaded:
            if debug_find_module:
                print("...pxd not loaded")
            if not pxd_pathname:
                if debug_find_module:
                    print("...looking for pxd file")
                # Only look in sys.path if we are explicitly looking
                # for a .pxd file.
                pxd_pathname = self.find_pxd_file(qualified_name, pos, sys_path=need_pxd and not relative_import)
                self._check_pxd_filename(pos, pxd_pathname, qualified_name)
                if debug_find_module:
                    print("......found %s" % pxd_pathname)
                if not pxd_pathname and need_pxd:
                    # Set pxd_file_loaded such that we don't need to
                    # look for the non-existing pxd file next time.
                    scope.pxd_file_loaded = True
                    package_pathname = self.search_include_directories(
                        qualified_name, suffix=".py", source_pos=pos, sys_path=not relative_import)
                    if package_pathname and package_pathname.endswith(Utils.PACKAGE_FILES):
                        pass
                    else:
                        error(pos, "'%s.pxd' not found" % qualified_name.replace('.', os.sep))
            if pxd_pathname:
                scope.pxd_file_loaded = True
                try:
                    if debug_find_module:
                        print("Context.find_module: Parsing %s" % pxd_pathname)
                    rel_path = module_name.replace('.', os.sep) + os.path.splitext(pxd_pathname)[1]
                    if not pxd_pathname.endswith(rel_path):
                        rel_path = pxd_pathname  # safety measure to prevent printing incorrect paths
                    source_desc = FileSourceDescriptor(pxd_pathname, rel_path)
                    err, result = self.process_pxd(source_desc, scope, qualified_name)
                    if err:
                        raise err
                    (pxd_codenodes, pxd_scope) = result
                    self.pxds[module_name] = (pxd_codenodes, pxd_scope)
                except CompileError:
                    pass
        return scope

    def find_pxd_file(self, qualified_name, pos=None, sys_path=True, source_file_path=None):
        # Search include path (and sys.path if sys_path is True) for
        # the .pxd file corresponding to the given fully-qualified
        # module name.
        # Will find either a dotted filename or a file in a
        # package directory. If a source file position is given,
        # the directory containing the source file is searched first
        # for a dotted filename, and its containing package root
        # directory is searched first for a non-dotted filename.
        pxd = self.search_include_directories(
            qualified_name, suffix=".pxd", source_pos=pos, sys_path=sys_path, source_file_path=source_file_path)
        if pxd is None and Options.cimport_from_pyx:
            return self.find_pyx_file(qualified_name, pos, sys_path=sys_path)
        return pxd

    def find_pyx_file(self, qualified_name, pos=None, sys_path=True, source_file_path=None):
        # Search include path for the .pyx file corresponding to the
        # given fully-qualified module name, as for find_pxd_file().
        return self.search_include_directories(
            qualified_name, suffix=".pyx", source_pos=pos, sys_path=sys_path, source_file_path=source_file_path)

    def find_include_file(self, filename, pos=None, source_file_path=None):
        # Search list of include directories for filename.
        # Reports an error and returns None if not found.
        path = self.search_include_directories(
            filename, source_pos=pos, include=True, source_file_path=source_file_path)
        if not path:
            error(pos, "'%s' not found" % filename)
        return path

    def search_include_directories(self, qualified_name,
                                   suffix=None, source_pos=None, include=False, sys_path=False, source_file_path=None):
        include_dirs = self.include_directories
        if sys_path:
            include_dirs = include_dirs + sys.path
        # include_dirs must be hashable for caching in @cached_function
        include_dirs = tuple(include_dirs + [standard_include_path])
        return search_include_directories(
            include_dirs, qualified_name, suffix or "", source_pos, include, source_file_path)

    def find_root_package_dir(self, file_path):
        return Utils.find_root_package_dir(file_path)

    def check_package_dir(self, dir, package_names):
        return Utils.check_package_dir(dir, tuple(package_names))

    def c_file_out_of_date(self, source_path, output_path):
        if not os.path.exists(output_path):
            return 1
        c_time = Utils.modification_time(output_path)
        if Utils.file_newer_than(source_path, c_time):
            return 1
        pxd_path = Utils.replace_suffix(source_path, ".pxd")
        if os.path.exists(pxd_path) and Utils.file_newer_than(pxd_path, c_time):
            return 1
        for kind, name in self.read_dependency_file(source_path):
            if kind == "cimport":
                dep_path = self.find_pxd_file(name, source_file_path=source_path)
            elif kind == "include":
                dep_path = self.search_include_directories(name, source_file_path=source_path)
            else:
                continue
            if dep_path and Utils.file_newer_than(dep_path, c_time):
                return 1
        return 0

    def find_cimported_module_names(self, source_path):
        return [ name for kind, name in self.read_dependency_file(source_path)
                 if kind == "cimport" ]

    def is_package_dir(self, dir_path):
        return Utils.is_package_dir(dir_path)

    def read_dependency_file(self, source_path):
        dep_path = Utils.replace_suffix(source_path, ".dep")
        if os.path.exists(dep_path):
            with open(dep_path) as f:
                chunks = [ line.split(" ", 1)
                           for line in (l.strip() for l in f)
                           if " " in line ]
            return chunks
        else:
            return ()

    def lookup_submodule(self, name):
        # Look up a top-level module. Returns None if not found.
        return self.modules.get(name, None)

    def find_submodule(self, name, as_package=False):
        # Find a top-level module, creating a new one if needed.
        scope = self.lookup_submodule(name)
        if not scope:
            scope = ModuleScope(name,
                parent_module = None, context = self, is_package=as_package)
            self.modules[name] = scope
        return scope

    def parse(self, source_desc, scope, pxd, full_module_name):
        if not isinstance(source_desc, FileSourceDescriptor):
            raise RuntimeError("Only file sources for code supported")
        scope.cpp = self.cpp
        # Parse the given source file and return a parse tree.
        num_errors = Errors.get_errors_count()
        try:
            with source_desc.get_file_object() as f:
                from . import Parsing
                s = PyrexScanner(f, source_desc, source_encoding = f.encoding,
                                 scope = scope, context = self)
                tree = Parsing.p_module(s, pxd, full_module_name)
                if self.options.formal_grammar:
                    try:
                        from ..Parser import ConcreteSyntaxTree
                    except ImportError:
                        raise RuntimeError(
                            "Formal grammar can only be used with compiled Cython with an available pgen.")
                    ConcreteSyntaxTree.p_module(source_desc.filename)
        except UnicodeDecodeError as e:
            #import traceback
            #traceback.print_exc()
            raise self._report_decode_error(source_desc, e)

        if Errors.get_errors_count() > num_errors:
            raise CompileError()
        return tree

    def _report_decode_error(self, source_desc, exc):
        msg = exc.args[-1]
        position = exc.args[2]
        encoding = exc.args[0]

        line = 1
        column = idx = 0
        with open(source_desc.filename, encoding='iso8859-1', newline='') as f:
            for line, data in enumerate(f, 1):
                idx += len(data)
                if idx >= position:
                    column = position - (idx - len(data)) + 1
                    break

        return error((source_desc, line, column),
                     "Decoding error, missing or incorrect coding=<encoding-name> "
                     "at top of source (cannot decode with encoding %r: %s)" % (encoding, msg))

    def extract_module_name(self, path, options):
        # Find fully_qualified module name from the full pathname
        # of a source file.
        dir, filename = os.path.split(path)
        module_name, _ = os.path.splitext(filename)
        if "." in module_name:
            return module_name
        names = [module_name]
        while self.is_package_dir(dir):
            parent, package_name = os.path.split(dir)
            if parent == dir:
                break
            names.append(package_name)
            dir = parent
        names.reverse()
        return ".".join(names)

    def setup_errors(self, options, result):
        Errors.init_thread()
        if options.use_listing_file:
            path = result.listing_file = Utils.replace_suffix(result.main_source_file, ".lis")
        else:
            path = None
        Errors.open_listing_file(path=path, echo_to_stderr=options.errors_to_stderr)

    def teardown_errors(self, err, options, result):
        source_desc = result.compilation_source.source_desc
        if not isinstance(source_desc, FileSourceDescriptor):
            raise RuntimeError("Only file sources for code supported")
        Errors.close_listing_file()
        result.num_errors = Errors.get_errors_count()
        if result.num_errors > 0:
            err = True
        if err and result.c_file:
            try:
                Utils.castrate_file(result.c_file, os.stat(source_desc.filename))
            except OSError:
                pass
            result.c_file = None


def get_output_filename(source_filename, cwd, options):
    if options.cplus:
        c_suffix = ".cpp"
    else:
        c_suffix = ".c"
    suggested_file_name = Utils.replace_suffix(source_filename, c_suffix)
    if options.output_file:
        out_path = os.path.join(cwd, options.output_file)
        if os.path.isdir(out_path):
            return os.path.join(out_path, os.path.basename(suggested_file_name))
        else:
            return out_path
    else:
        return suggested_file_name


def create_default_resultobj(compilation_source, options):
    result = CompilationResult()
    result.main_source_file = compilation_source.source_desc.filename
    result.compilation_source = compilation_source
    source_desc = compilation_source.source_desc
    result.c_file = get_output_filename(source_desc.filename,
                        compilation_source.cwd, options)
    result.embedded_metadata = options.embedded_metadata
    return result


def setup_source_object(source, source_ext, full_module_name, options, context):
    cwd = os.getcwd()
    abs_path = os.path.abspath(source)

    full_module_name = full_module_name or context.extract_module_name(source, options)
    full_module_name = EncodedString(full_module_name)

    Utils.raise_error_if_module_name_forbidden(full_module_name)
    if options.relative_path_in_code_position_comments:
        rel_path = full_module_name.replace('.', os.sep) + source_ext
        if not abs_path.endswith(rel_path):
            rel_path = source  # safety measure to prevent printing incorrect paths
    else:
        rel_path = abs_path
    source_desc = FileSourceDescriptor(abs_path, rel_path)
    return CompilationSource(source_desc, full_module_name, cwd)


def run_cached_pipeline(source, options, full_module_name, context, cache, fingerprint):
    cwd = os.getcwd()
    output_filename = get_output_filename(source, cwd, options)
    cached = cache.lookup_cache(output_filename, fingerprint)
    if cached:
        cache.load_from_cache(output_filename, cached)

        source_ext = os.path.splitext(source)[1]
        options.configure_language_defaults(source_ext[1:])  # py/pyx

        source = setup_source_object(source, source_ext, full_module_name, options, context)
        # Set up result object
        return create_default_resultobj(source, options)

    result = run_pipeline(source, options, full_module_name, context)
    if fingerprint:
        cache.store_to_cache(output_filename, fingerprint, result)
    return result


def run_pipeline(source, options, full_module_name, context):
    from . import Pipeline
    if options.verbose:
        sys.stderr.write("Compiling %s\n" % source)
    source_ext = os.path.splitext(source)[1]
    abs_path = os.path.abspath(source)
    options.configure_language_defaults(source_ext[1:])  # py/pyx

    source = setup_source_object(source, source_ext, full_module_name, options, context)
    # Set up result object
    result = create_default_resultobj(source, options)

    if options.annotate is None:
        # By default, decide based on whether an html file already exists.
        html_filename = os.path.splitext(result.c_file)[0] + ".html"
        if os.path.exists(html_filename):
            with open(html_filename, encoding="UTF-8") as html_file:
                if '<!-- Generated by Cython' in html_file.read(100):
                    options.annotate = True

    # Get pipeline
    if source_ext.lower() == '.py' or not source_ext:
        pipeline = Pipeline.create_py_pipeline(context, options, result)
    else:
        pipeline = Pipeline.create_pyx_pipeline(context, options, result)

    context.setup_errors(options, result)

    if '.' in source.full_module_name and '.' in os.path.splitext(os.path.basename(abs_path))[0]:
        warning((source.source_desc, 1, 0),
                "Dotted filenames ('%s') are deprecated."
                " Please use the normal Python package directory layout." % os.path.basename(abs_path), level=1)
    if re.search("[.]c(pp|[+][+]|xx)$", result.c_file, re.RegexFlag.IGNORECASE) and not context.cpp:
        warning((source.source_desc, 1, 0),
                "Filename implies a c++ file but Cython is not in c++ mode.",
                level=1)

    err, enddata = Pipeline.run_pipeline(pipeline, source)
    context.teardown_errors(err, options, result)
    if err is None and options.depfile:
        from ..Build.Dependencies import create_dependency_tree
        dependencies = create_dependency_tree(context).all_dependencies(result.main_source_file)
        Utils.write_depfile(result.c_file, result.main_source_file, dependencies)
    return result


# ------------------------------------------------------------------------
#
#  Main Python entry points
#
# ------------------------------------------------------------------------

class CompilationSource:
    """
    Contains the data necessary to start up a compilation pipeline for
    a single compilation unit.
    """
    def __init__(self, source_desc, full_module_name, cwd):
        self.source_desc = source_desc
        self.full_module_name = full_module_name
        self.cwd = cwd


class CompilationResult:
    """
    Results from the Cython compiler:

    c_file           string or None   The generated C source file
    h_file           string or None   The generated C header file
    i_file           string or None   The generated .pxi file
    api_file         string or None   The generated C API .h file
    listing_file     string or None   File of error messages
    object_file      string or None   Result of compiling the C file
    extension_file   string or None   Result of linking the object file
    num_errors       integer          Number of compilation errors
    compilation_source CompilationSource
    """

    c_file = None
    h_file = None
    i_file = None
    api_file = None
    listing_file = None
    object_file = None
    extension_file = None
    main_source_file = None
    num_errors = 0

    def get_generated_source_files(self):
        return [
            source_file for source_file in [self.c_file, self.h_file, self.i_file, self.api_file]
            if source_file
        ]


class CompilationResultSet(dict):
    """
    Results from compiling multiple Pyrex source files. A mapping
    from source file paths to CompilationResult instances. Also
    has the following attributes:

    num_errors   integer   Total number of compilation errors
    """

    num_errors = 0

    def add(self, source, result):
        self[source] = result
        self.num_errors += result.num_errors


def get_fingerprint(cache, source, options):
    from ..Build.Dependencies import create_dependency_tree
    from ..Build.Cache import FingerprintFlags
    context = Context.from_options(options)
    dependencies = create_dependency_tree(context)
    return cache.transitive_fingerprint(
            source, dependencies.all_dependencies(source), options,
            FingerprintFlags(
                'c++' if options.cplus else 'c',
                np_pythran=options.np_pythran
            )
    )


def compile_single(source, options, full_module_name, cache=None, context=None, fingerprint=None):
    """
    compile_single(source, options, full_module_name, cache, context, fingerprint)

    Compile the given Pyrex implementation file and return a CompilationResult.
    Always compiles a single file; does not perform timestamp checking or
    recursion.
    """

    if context is None:
        context = Context.from_options(options)

    if cache:
        fingerprint = fingerprint or get_fingerprint(cache, source, options)
        return run_cached_pipeline(source, options, full_module_name, context, cache, fingerprint)
    else:
        return run_pipeline(source, options, full_module_name, context)


def compile_multiple(sources, options, cache=None):
    """
    compile_multiple(sources, options, cache)

    Compiles the given sequence of Pyrex implementation files and returns
    a CompilationResultSet. Performs timestamp checking, caching and/or recursion
    if these are specified in the options.
    """
    if len(sources) > 1 and options.module_name:
        raise RuntimeError('Full module name can only be set '
                           'for single source compilation')
    # run_pipeline creates the context
    # context = Context.from_options(options)
    sources = [os.path.abspath(source) for source in sources]
    processed = set()
    results = CompilationResultSet()
    timestamps = options.timestamps
    context = None
    cwd = os.getcwd()
    for source in sources:
        if source not in processed:
            output_filename = get_output_filename(source, cwd, options)
            if context is None:
                context = Context.from_options(options)
            out_of_date = context.c_file_out_of_date(source, output_filename)
            if (not timestamps) or out_of_date:
                result = compile_single(source, options, full_module_name=options.module_name, cache=cache, context=context)
                results.add(source, result)
                # Compiling multiple sources in one context doesn't quite
                # work properly yet.
                context = None
            processed.add(source)
    if cache:
        cache.cleanup_cache()
    return results


def compile(source, options = None, full_module_name = None, **kwds):
    """
    compile(source [, options], [, <option> = <value>]...)

    Compile one or more Pyrex implementation files, with optional timestamp
    checking and recursing on dependencies.  The source argument may be a string
    or a sequence of strings.  If it is a string and no recursion or timestamp
    checking is requested, a CompilationResult is returned, otherwise a
    CompilationResultSet is returned.
    """
    options = CompilationOptions(defaults = options, **kwds)

    # cache is enabled when:
    # * options.cache is True (the default path to the cache base dir is used)
    # * options.cache is the explicit path to the cache base dir
    # unless annotations are generated
    cache = None
    if options.cache:
        if options.annotate or Options.annotate:
            if options.verbose:
                sys.stderr.write('Cache is ignored when annotations are enabled.\n')
        else:
            from ..Build.Cache import Cache
            cache_path = None if options.cache is True else options.cache
            cache = Cache(cache_path)

    if isinstance(source, str):
        if not options.timestamps:
            return compile_single(source, options, full_module_name, cache)
        source = [source]
    return compile_multiple(source, options, cache)


@Utils.cached_function
def search_include_directories(dirs, qualified_name, suffix="", pos=None, include=False, source_file_path=None):
    """
    Search the list of include directories for the given file name.

    If a source file path or position is given, first searches the directory
    containing that file.  Returns None if not found, but does not report an error.

    The 'include' option will disable package dereferencing.
    """
    if pos and not source_file_path:
        file_desc = pos[0]
        if not isinstance(file_desc, FileSourceDescriptor):
            raise RuntimeError("Only file sources for code supported")
        source_file_path = file_desc.filename
    if source_file_path:
        if include:
            dirs = (os.path.dirname(source_file_path),) + dirs
        else:
            dirs = (Utils.find_root_package_dir(source_file_path),) + dirs

    # search for dotted filename e.g. <dir>/foo.bar.pxd
    dotted_filename = qualified_name
    if suffix:
        dotted_filename += suffix

    for dirname in dirs:
        path = os.path.join(dirname, dotted_filename)
        if os.path.exists(path):
            return path

    # search for filename in package structure e.g. <dir>/foo/bar.pxd or <dir>/foo/bar/__init__.pxd
    if not include:

        names = qualified_name.split('.')
        package_names = tuple(names[:-1])
        module_name = names[-1]

        # search for standard packages first - PEP420
        namespace_dirs = []
        for dirname in dirs:
            package_dir, is_namespace = Utils.check_package_dir(dirname, package_names)
            if package_dir is not None:
                if is_namespace:
                    namespace_dirs.append(package_dir)
                    continue
                path = search_module_in_dir(package_dir, module_name, suffix)
                if path:
                    return path

        # search for namespaces second - PEP420
        for package_dir in namespace_dirs:
            path = search_module_in_dir(package_dir, module_name, suffix)
            if path:
                return path

    return None


@Utils.cached_function
def search_module_in_dir(package_dir, module_name, suffix):
    # matches modules of the form: <dir>/foo/bar.pxd
    path = Utils.find_versioned_file(package_dir, module_name, suffix)

    # matches modules of the form: <dir>/foo/bar/__init__.pxd
    if not path and suffix:
        path = Utils.find_versioned_file(os.path.join(package_dir, module_name), "__init__", suffix)

    return path


# ------------------------------------------------------------------------
#
#  Main command-line entry point
#
# ------------------------------------------------------------------------

def setuptools_main():
    return main(command_line = 1)


def main(command_line = 0):
    args = sys.argv[1:]
    any_failures = 0
    if command_line:
        try:
            options, sources = parse_command_line(args)
        except FileNotFoundError as e:
            print("{}: No such file or directory: '{}'".format(sys.argv[0], e.filename), file=sys.stderr)
            sys.exit(1)
    else:
        options = CompilationOptions(default_options)
        sources = args

    if options.show_version:
        Utils.print_version()

    if options.working_path!="":
        os.chdir(options.working_path)

    try:
        if options.shared_c_file_path:
            from ..Build.SharedModule import generate_shared_module
            generate_shared_module(options)
            return

        result = compile(sources, options)
        if result.num_errors > 0:
            any_failures = 1
    except (OSError, PyrexError) as e:
        sys.stderr.write(str(e) + '\n')
        any_failures = 1
    if any_failures:
        sys.exit(1)
