import gc
import hashlib
import inspect
import os
import re
import sys
import time

from distutils.core import Distribution, Extension
from distutils.command.build_ext import build_ext

import Cython
from ..Compiler.Main import Context
from ..Compiler.Options import (default_options, CompilationOptions,
    get_directive_defaults)

from ..Compiler.Visitor import CythonTransform, EnvTransform
from ..Compiler.ParseTreeTransforms import SkipDeclarations
from ..Compiler.TreeFragment import parse_from_strings
from .Dependencies import strip_string_literals, cythonize, cached_function
from .Cache import get_cython_cache_dir
from ..Compiler import Pipeline
import cython as cython_module

import importlib.util
from importlib.machinery import ExtensionFileLoader

def load_dynamic(name, path):
    spec = importlib.util.spec_from_file_location(name, loader=ExtensionFileLoader(name, path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class UnboundSymbols(EnvTransform, SkipDeclarations):
    def __init__(self):
        super(EnvTransform, self).__init__(context=None)
        self.unbound = set()
    def visit_NameNode(self, node):
        if not self.current_env().lookup(node.name):
            self.unbound.add(node.name)
        return node
    def __call__(self, node):
        super().__call__(node)
        return self.unbound


@cached_function
def unbound_symbols(code, context=None):
    if context is None:
        context = Context([], get_directive_defaults(),
                          options=CompilationOptions(default_options))
    from ..Compiler.ParseTreeTransforms import AnalyseDeclarationsTransform
    tree = parse_from_strings('(tree fragment)', code)
    for phase in Pipeline.create_pipeline(context, 'pyx'):
        if phase is None:
            continue
        tree = phase(tree)
        if isinstance(phase, AnalyseDeclarationsTransform):
            break
    import builtins
    return tuple(UnboundSymbols()(tree) - set(dir(builtins)))


def unsafe_type(arg, context=None):
    py_type = type(arg)
    if py_type is int:
        return 'long'
    else:
        return safe_type(arg, context)


def safe_type(arg, context=None):
    py_type = type(arg)
    if py_type in (list, tuple, dict, str):
        return py_type.__name__
    elif py_type is complex:
        return 'double complex'
    elif py_type is float:
        return 'double'
    elif py_type is bool:
        return 'bint'
    elif 'numpy' in sys.modules and isinstance(arg, sys.modules['numpy'].ndarray):
        return 'numpy.ndarray[numpy.%s_t, ndim=%s]' % (arg.dtype.name, arg.ndim)
    else:
        for base_type in py_type.__mro__:
            if base_type.__module__ in ('__builtin__', 'builtins'):
                return 'object'
            module = context.find_module(base_type.__module__, need_pxd=False)
            if module:
                entry = module.lookup(base_type.__name__)
                if entry.is_type:
                    return '%s.%s' % (base_type.__module__, base_type.__name__)
        return 'object'


def _get_build_extension():
    dist = Distribution()
    # Ensure the build respects distutils configuration by parsing
    # the configuration files
    config_files = dist.find_config_files()
    dist.parse_config_files(config_files)
    build_extension = build_ext(dist)
    build_extension.finalize_options()
    return build_extension


@cached_function
def _create_context(cython_include_dirs):
    return Context(
        list(cython_include_dirs),
        get_directive_defaults(),
        options=CompilationOptions(default_options)
    )


_cython_inline_cache = {}
_cython_inline_default_context = _create_context(('.',))


def _populate_unbound(kwds, unbound_symbols, locals=None, globals=None):
    for symbol in unbound_symbols:
        if symbol not in kwds:
            if locals is None or globals is None:
                calling_frame = inspect.currentframe().f_back.f_back.f_back
                if locals is None:
                    locals = calling_frame.f_locals
                if globals is None:
                    globals = calling_frame.f_globals
            if not isinstance(locals, dict):
                # FrameLocalsProxy is stricter than dict on how it looks up keys
                # and this means our "EncodedStrings" don't match the keys in locals.
                # Therefore copy to a dict.
                locals = dict(locals)
            if symbol in locals:
                kwds[symbol] = locals[symbol]
            elif symbol in globals:
                kwds[symbol] = globals[symbol]
            else:
                print("Couldn't find %r" % symbol)


def _inline_key(orig_code, arg_sigs, language_level):
    key = orig_code, arg_sigs, sys.version_info, sys.executable, language_level, Cython.__version__
    return hashlib.sha256(str(key).encode('utf-8')).hexdigest()


def cython_inline(code, get_type=unsafe_type,
                  lib_dir=os.path.join(get_cython_cache_dir(), 'inline'),
                  cython_include_dirs=None, cython_compiler_directives=None,
                  force=False, quiet=False, locals=None, globals=None, language_level=None, **kwds):

    if get_type is None:
        get_type = lambda x: 'object'
    ctx = _create_context(tuple(cython_include_dirs)) if cython_include_dirs else _cython_inline_default_context

    cython_compiler_directives = dict(cython_compiler_directives) if cython_compiler_directives else {}
    if language_level is None and 'language_level' not in cython_compiler_directives:
        language_level = '3'
    if language_level is not None:
        cython_compiler_directives['language_level'] = language_level

    key_hash = None

    # Fast path if this has been called in this session.
    _unbound_symbols = _cython_inline_cache.get(code)
    if _unbound_symbols is not None:
        _populate_unbound(kwds, _unbound_symbols, locals, globals)
        args = sorted(kwds.items())
        arg_sigs = tuple([(get_type(value, ctx), arg) for arg, value in args])
        key_hash = _inline_key(code, arg_sigs, language_level)
        invoke = _cython_inline_cache.get((code, arg_sigs, key_hash))
        if invoke is not None:
            arg_list = [arg[1] for arg in args]
            return invoke(*arg_list)

    orig_code = code
    code, literals = strip_string_literals(code)
    code = strip_common_indent(code)
    if locals is None:
        locals = inspect.currentframe().f_back.f_back.f_locals
    if globals is None:
        globals = inspect.currentframe().f_back.f_back.f_globals
    try:
        _cython_inline_cache[orig_code] = _unbound_symbols = unbound_symbols(code)
        _populate_unbound(kwds, _unbound_symbols, locals, globals)
    except AssertionError:
        if not quiet:
            # Parsing from strings not fully supported (e.g. cimports).
            print("Could not parse code as a string (to extract unbound symbols).")

    cimports = []
    for name, arg in list(kwds.items()):
        if arg is cython_module:
            cimports.append('\ncimport cython as %s' % name)
            del kwds[name]
    arg_names = sorted(kwds)
    arg_sigs = tuple([(get_type(kwds[arg], ctx), arg) for arg in arg_names])
    if key_hash is None:
        key_hash = _inline_key(orig_code, arg_sigs, language_level)
    module_name = "_cython_inline_" + key_hash

    if module_name in sys.modules:
        module = sys.modules[module_name]

    else:
        build_extension = None
        if cython_inline.so_ext is None:
            # Figure out and cache current extension suffix
            build_extension = _get_build_extension()
            cython_inline.so_ext = build_extension.get_ext_filename('')

        lib_dir = os.path.abspath(lib_dir)
        module_path = os.path.join(lib_dir, module_name + cython_inline.so_ext)

        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        if force or not os.path.isfile(module_path):
            cflags = []
            define_macros = []
            c_include_dirs = []
            qualified = re.compile(r'([.\w]+)[.]')
            for type, _ in arg_sigs:
                m = qualified.match(type)
                if m:
                    cimports.append('\ncimport %s' % m.groups()[0])
                    # one special case
                    if m.groups()[0] == 'numpy':
                        import numpy
                        c_include_dirs.append(numpy.get_include())
                        define_macros.append(("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"))
                        # cflags.append('-Wno-unused')
            module_body, func_body = extract_func_code(code)
            params = ', '.join(['%s %s' % a for a in arg_sigs])
            module_code = """
%(module_body)s
%(cimports)s
def __invoke(%(params)s):
%(func_body)s
    return locals()
            """ % {'cimports': '\n'.join(cimports),
                   'module_body': module_body,
                   'params': params,
                   'func_body': func_body }
            for key, value in literals.items():
                module_code = module_code.replace(key, value)
            pyx_file = os.path.join(lib_dir, module_name + '.pyx')
            fh = open(pyx_file, 'w')
            try:
                fh.write(module_code)
            finally:
                fh.close()
            extension = Extension(
                name=module_name,
                sources=[pyx_file],
                include_dirs=c_include_dirs or None,
                extra_compile_args=cflags or None,
                define_macros=define_macros or None,
            )
            if build_extension is None:
                build_extension = _get_build_extension()
            build_extension.extensions = cythonize(
                [extension],
                include_path=cython_include_dirs or ['.'],
                compiler_directives=cython_compiler_directives,
                quiet=quiet)
            build_extension.build_temp = os.path.dirname(pyx_file)
            build_extension.build_lib  = lib_dir
            build_extension.run()

        if sys.platform == 'win32' and sys.version_info >= (3, 8):
            with os.add_dll_directory(os.path.abspath(lib_dir)):
                module = load_dynamic(module_name, module_path)
        else:
            module = load_dynamic(module_name, module_path)

    _cython_inline_cache[orig_code, arg_sigs, key_hash] = module.__invoke
    arg_list = [kwds[arg] for arg in arg_names]
    return module.__invoke(*arg_list)


# The code template used for cymeit benchmark runs.
# We keep the benchmark repetition separate from the benchmarked code
# to prevent the C compiler from doing unhelpful loop optimisations.
_CYMEIT_TEMPLATE = """
def __PYX_repeat_benchmark(benchmark, timer, size_t number):
    cdef size_t i

    t0 = timer()
    for i in range(number):
        benchmark()
    t1 = timer()
    return t1 - t0

def __PYX_make_benchmark():
    {setup_code}

    def __PYX_run_benchmark():
        {benchmark_code}

    return __PYX_run_benchmark
"""


def cymeit(code, setup_code=None, import_module=None, directives=None, timer=time.perf_counter, repeat=9):
    """Benchmark a Cython code string similar to 'timeit'.

    'setup_code': string of setup code that will be run before taking the timings.

    'import_module': a module namespace to run the benchmark in
                     (usually a compiled Cython module).

    'directives': Cython directives to use when compiling the benchmark code.

    'timer': The timer function. Defaults to 'time.perf_counter', returning float seconds.
             Nanosecond timers are detected (and can only be used) if they return integers.

    'repeat': The number of timings to take and return.

    Returns a tuple: (list of single-loop timings, number of loops run for each)
    """
    import textwrap

    # Compile the benchmark code as an inline closure function.

    setup_code = strip_common_indent(setup_code) if setup_code else ''
    code = strip_common_indent(code) if code.strip() else 'pass'

    module_namespace = __import__(import_module).__dict__ if import_module else None

    cymeit_code = _CYMEIT_TEMPLATE.format(
        setup_code=textwrap.indent(setup_code, ' '*4).strip(),
        benchmark_code=textwrap.indent(code, ' '*8).strip(),

    )

    namespace = cython_inline(
        cymeit_code,
        cython_compiler_directives=directives,
        locals=module_namespace,
    )

    make_benchmark = namespace['__PYX_make_benchmark']
    repeat_benchmark = namespace['__PYX_repeat_benchmark']

    # Based on 'timeit' in CPython 3.13.

    def timeit(number):
        benchmark = make_benchmark()

        gcold = gc.isenabled()
        gc.disable()
        try:
            timing = repeat_benchmark(benchmark, timer, number)
        finally:
            if gcold:
                gc.enable()
        return timing

    # Find a sufficiently large number of loops, warm up the system.

    timer_returns_nanoseconds = isinstance(timer(), int)
    one_second = 1_000_000_000 if timer_returns_nanoseconds else 1.0

    # Run for at least 0.2 seconds, either as integer nanoseconds or floating point seconds.
    min_runtime = one_second // 5 if timer_returns_nanoseconds else one_second / 5

    def autorange():
        i = 1
        while True:
            for j in 1, 2, 5:
                number = i * j
                time_taken = timeit(number)
                assert isinstance(time_taken, int if timer_returns_nanoseconds else float)
                if time_taken >= min_runtime:
                    return number
                elif timer_returns_nanoseconds and (time_taken < 10 and number >= 10):
                    # Arbitrary sanity check to prevent endless loops for non-ns timers.
                    raise RuntimeError(f"Timer seems to return non-ns timings: {timer}")
            i *= 10

    autorange()  # warmup
    number = autorange()

    # Run and repeat the benchmark.
    timings = [
        timeit(number)
        for _ in range(repeat)
    ]

    half = number // 2  # for integer rounding

    timings = [
        (timing + half) // number if timer_returns_nanoseconds else timing / number
        for timing in timings
    ]

    return (timings, number)


# Cached suffix used by cython_inline above.  None should get
# overridden with actual value upon the first cython_inline invocation
cython_inline.so_ext = None

_find_non_space = re.compile(r'\S').search


def strip_common_indent(code):
    min_indent = None
    lines = code.splitlines()
    for line in lines:
        match = _find_non_space(line)
        if not match:
            continue  # blank
        indent = match.start()
        if line[indent] == '#':
            continue  # comment
        if min_indent is None or min_indent > indent:
            min_indent = indent
    for ix, line in enumerate(lines):
        match = _find_non_space(line)
        if not match or not line or line[indent:indent+1] == '#':
            continue
        lines[ix] = line[min_indent:]
    return '\n'.join(lines)


module_statement = re.compile(r'^((cdef +(extern|class))|cimport|(from .+ cimport)|(from .+ import +[*]))')
def extract_func_code(code):
    module = []
    function = []
    current = function
    code = code.replace('\t', ' ')
    lines = code.split('\n')
    for line in lines:
        if not line.startswith(' '):
            if module_statement.match(line):
                current = module
            else:
                current = function
        current.append(line)
    return '\n'.join(module), '    ' + '\n    '.join(function)


def get_body(source):
    ix = source.index(':')
    if source[:5] == 'lambda':
        return "return %s" % source[ix+1:]
    else:
        return source[ix+1:]


# Lots to be done here... It would be especially cool if compiled functions
# could invoke each other quickly.
class RuntimeCompiledFunction:

    def __init__(self, f):
        self._f = f
        self._body = get_body(inspect.getsource(f))

    def __call__(self, *args, **kwds):
        all = inspect.getcallargs(self._f, *args, **kwds)
        return cython_inline(self._body, locals=self._f.__globals__, globals=self._f.__globals__, **all)
