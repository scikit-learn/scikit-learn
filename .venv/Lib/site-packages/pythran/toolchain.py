'''
This module contains all the stuff to make your way from python code to
a dynamic library, see __init__.py for exported interfaces.
'''

from pythran.backend import Cxx, Python
from pythran.config import cfg
from pythran.cxxgen import PythonModule, Include, Line, Statement
from pythran.cxxgen import FunctionBody, FunctionDeclaration, Value, Block
from pythran.cxxgen import ReturnStatement
from pythran.errors import PythranCompileError
from pythran.middlend import refine, mark_unexported_functions
from pythran.passmanager import PassManager
from pythran.tables import pythran_ward
from pythran.types.type_dependencies import pytype_to_deps
from pythran.types.conversion import pytype_to_ctype
from pythran.spec import load_specfile, Spec
from pythran.spec import spec_to_string
from pythran.syntax import check_specs, check_exports, PythranSyntaxError
from pythran.version import __version__
from pythran.utils import cxxid
import pythran.frontend as frontend

import sysconfig

from tempfile import mkdtemp, NamedTemporaryFile
import gast as ast
import importlib
import logging
import os.path
import shutil
import glob
import hashlib
from functools import reduce
import sys

logger = logging.getLogger('pythran')


def _extract_specs_dependencies(specs):
    """ Extract types dependencies from specs for each exported signature. """
    deps = set()
    # for each function
    for signatures in specs.functions.values():
        # for each signature
        for signature in signatures:
            # for each argument
            for t in signature:
                deps.update(pytype_to_deps(t))
    # and each capsule
    for signature in specs.capsules.values():
        # for each argument
        for t in signature:
            deps.update(pytype_to_deps(t))

    return sorted(deps)


def _parse_optimization(optimization):
    '''Turns an optimization of the form
        my_optim
        my_package.my_optim
        into the associated symbol'''
    splitted = optimization.split('.')
    if len(splitted) == 1:
        splitted = ['pythran', 'optimizations'] + splitted
    return reduce(getattr, splitted[1:], __import__(splitted[0]))


def _write_temp(content, suffix):
    '''write `content` to a temporary XXX`suffix` file and return the filename.
       It is user's responsibility to delete when done.'''
    with NamedTemporaryFile(mode='w', suffix=suffix, delete=False) as out:
        out.write(content)
        return out.name


def has_argument(module, fname):
    '''Checks if a given function has arguments'''
    for n in module.body:
        if isinstance(n, ast.FunctionDef) and n.name == fname:
            return [cxxid(arg.id) for arg in n.args.args]
    return []


def front_middle_end(module_name, code, optimizations=None, module_dir=None,
                     entry_points=None, report_times = False):
    """Front-end and middle-end compilation steps"""

    pm = PassManager(module_name, module_dir, code)

    # front end
    ir, docstrings = frontend.parse(pm, code)

    if entry_points is not None:
        ir = mark_unexported_functions(ir, entry_points)

    # middle-end
    if optimizations is None:
        optimizations = cfg.get('pythran', 'optimizations').split()
    optimizations = [_parse_optimization(opt) for opt in optimizations]
    refine(pm, ir, optimizations, report_times)

    return pm, ir, docstrings


# PUBLIC INTERFACE STARTS HERE


def generate_py(module_name, code, optimizations=None, module_dir=None, report_times=False):
    '''python + pythran spec -> py code

    Prints and returns the optimized python code.

    '''

    pm, ir, _ = front_middle_end(module_name, code, optimizations, module_dir, report_times=report_times)
    return pm.dump(Python, ir)


def generate_cxx(module_name, code, specs=None, optimizations=None,
                 module_dir=None, report_times=False):
    '''python + pythran spec -> c++ code
    returns a PythonModule object and an error checker

    the error checker can be used to print more detailed info on the origin of
    a compile error (e.g. due to bad typing)

    '''
    if specs:
        entry_points = set(specs.keys())
    else:
        entry_points = None

    pm, ir, docstrings = front_middle_end(module_name, code, optimizations,
                                          module_dir,
                                          report_times=report_times,
                                          entry_points=entry_points)

    # back-end
    content = pm.dump(Cxx, ir)

    # instantiate the meta program
    if specs is None:

        class Generable(object):
            def __init__(self, content):
                self.content = content

            def __str__(self):
                return str(self.content)

            generate = __str__

        mod = Generable(content)

        def error_checker():
            from pythran.types import tog
            tog.typecheck(ir)

    else:

        # uniform typing
        if isinstance(specs, dict):
            specs = Spec(specs, {})

        def error_checker():
            from pythran.types import tog
            types = tog.typecheck(ir)
            check_specs(specs, types)

        specs.to_docstrings(docstrings)
        check_exports(pm, ir, specs)

        if isinstance(code, bytes):
            code_bytes = code
        else:
            code_bytes = code.encode('ascii', 'ignore')
        metainfo = {'hash': hashlib.sha256(code_bytes).hexdigest(),
                    'version': __version__}

        mod = PythonModule(module_name, docstrings, metainfo)
        mod.add_to_includes(
            Include("pythonic/core.hpp"),
            Include("pythonic/python/core.hpp"),
            # FIXME: only include these when needed
            Include("pythonic/types/bool.hpp"),
            Include("pythonic/types/int.hpp"),
            Line("#ifdef _OPENMP\n#include <omp.h>\n#endif")
        )
        mod.add_to_includes(*[Include(inc) for inc in
                              _extract_specs_dependencies(specs)])
        mod.add_to_includes(*content.body)
        mod.add_to_includes(
            Include("pythonic/python/exception_handler.hpp"),
        )

        def warded(module_name, internal_name):
            return pythran_ward + '{0}::{1}'.format(module_name, internal_name)

        for function_name, signatures in specs.functions.items():
            internal_func_name = cxxid(function_name)
            # global variables are functions with no signatures :-)
            if not signatures:
                mod.add_global_var(function_name,
                                   "{}()()".format(warded(module_name,
                                                          internal_func_name)))

            for sigid, signature in enumerate(signatures):
                numbered_function_name = "{0}{1}".format(internal_func_name,
                                                         sigid)
                arguments_types = [pytype_to_ctype(t) for t in signature]
                arguments_names = has_argument(ir, function_name)
                arguments = [n for n, _ in
                             zip(arguments_names, arguments_types)]
                name_fmt = pythran_ward + "{0}::{1}::type{2}"
                args_list = ", ".join(arguments_types)
                specialized_fname = name_fmt.format(module_name,
                                                    internal_func_name,
                                                    "<{0}>".format(args_list)
                                                    if arguments_names else "")
                result_type = "typename %s::result_type" % specialized_fname
                mod.add_pyfunction(
                    FunctionBody(
                        FunctionDeclaration(
                            Value(
                                result_type,
                                numbered_function_name),
                            [Value(t + '&&', a)
                             for t, a in zip(arguments_types, arguments)]),
                        Block([Statement("""
                            PyThreadState *_save = PyEval_SaveThread();
                            try {{
                                auto res = {0}()({1});
                                PyEval_RestoreThread(_save);
                                return res;
                            }}
                            catch(...) {{
                                PyEval_RestoreThread(_save);
                                throw;
                            }}
                            """.format(warded(module_name,
                                              internal_func_name),
                                       ', '.join(arguments)))])
                    ),
                    function_name,
                    arguments_types,
                    signature
                )

        for function_name, signature in specs.capsules.items():
            internal_func_name = cxxid(function_name)

            arguments_types = [pytype_to_ctype(t) for t in signature]
            arguments_names = has_argument(ir, function_name)
            arguments = [n for n, _ in
                         zip(arguments_names, arguments_types)]
            name_fmt = pythran_ward + "{0}::{1}::type{2}"
            args_list = ", ".join(arguments_types)
            specialized_fname = name_fmt.format(module_name,
                                                internal_func_name,
                                                "<{0}>".format(args_list)
                                                if arguments_names else "")
            result_type = "typename %s::result_type" % specialized_fname
            docstring = spec_to_string(function_name, signature)
            mod.add_capsule(
                FunctionBody(
                    FunctionDeclaration(
                        Value(result_type, function_name),
                        [Value(t, a)
                         for t, a in zip(arguments_types, arguments)]),
                    Block([ReturnStatement("{0}()({1})".format(
                        warded(module_name, internal_func_name),
                        ', '.join(arguments)))])
                ),
                function_name,
                docstring
            )

        if specs.ufuncs:
            mod.add_to_includes(
                Include("pythonic/include/types/numpy_ufunc.hpp"),
            )

        for function_name, signatures in specs.ufuncs.items():
            internal_func_name = cxxid(function_name)

            for signature in signatures:
                arguments_types = [pytype_to_ctype(t) for t in signature]
                numpy_types = ['pythonic::c_type_to_numpy_type<{}>::value'.format(t) for t in
                               arguments_types]
                arguments_names = has_argument(ir, function_name)
                arguments = [n for n, _ in
                             zip(arguments_names, arguments_types)]
                name_fmt = pythran_ward + "{0}::{1}::type{2}"
                args_list = ", ".join(arguments_types)
                specialized_fname = name_fmt.format(module_name,
                                                    internal_func_name,
                                                    "<{0}>".format(args_list)
                                                    if arguments_names else "")
                result_type = "typename %s::result_type" % specialized_fname
                numpy_result_type = 'pythonic::c_type_to_numpy_type<{}>::value'.format(result_type)
                numpy_types.append(numpy_result_type)

                mod.add_ufunc(
                    FunctionBody(
                        FunctionDeclaration(
                            Value(result_type, function_name),
                            [Value(t, a)
                             for t, a in zip(arguments_types, arguments)]),
                        Block([ReturnStatement("{0}()({1})".format(
                            warded(module_name, internal_func_name),
                            ', '.join(arguments)))])
                    ),
                    function_name,
                    pythran_ward + module_name + "::" + internal_func_name,
                    arguments_types + [result_type],
                    numpy_types
                )

    return mod, error_checker
    return mod, error_checker


def compile_cxxfile(module_name, cxxfile, output_binary=None, **kwargs):
    '''c++ file -> native module
    Return the filename of the produced shared library
    Raises PythranCompileError on failure

    '''
    # local import so that we don't depend on setuptools for the code generation
    # part
    from pythran.dist import PythranExtension, PythranBuildExt

    from setuptools import setup

    builddir = mkdtemp()
    buildtmp = mkdtemp()

    extension = PythranExtension(module_name,
                                 [cxxfile],
                                 **kwargs)

    try:
        setup(name=module_name,
              ext_modules=[extension],
              cmdclass={"build_ext": PythranBuildExt},
              # fake CLI call
              script_name='setup.py',
              script_args=['--verbose'
                           if logger.isEnabledFor(logging.INFO)
                           else '--quiet',
                           'build_ext',
                           '--build-lib', builddir,
                           '--build-temp', buildtmp]
              )
    except SystemExit as e:
        raise PythranCompileError(str(e))

    def copy(src_file, dest_file):
        # not using shutil.copy because it fails to copy stat across devices
        with open(src_file, 'rb') as src:
            with open(dest_file, 'wb') as dest:
                dest.write(src.read())

    ext = sysconfig.get_config_var('EXT_SUFFIX')
    if getattr(extension, 'py_limited_api', False):
        _, ext = os.path.splitext(ext)
        ext = f".abi3{ext}"

    # Copy all generated files including the module name prefix (.pdb, ...)
    for f in glob.glob(os.path.join(builddir, module_name + "*")):
        if f.endswith(ext):
            if output_binary:
                output_binary = output_binary.replace('%{ext}', ext)
            else:
                output_binary = os.path.join(os.getcwd(), module_name + ext)
            copy(f, output_binary)
        else:
            if output_binary:
                output_binary = output_binary.replace('%{ext}', '')
                output_directory = os.path.dirname(output_binary)
            else:
                output_directory = os.getcwd()
            copy(f, os.path.join(output_directory, os.path.basename(f)))
    shutil.rmtree(builddir)
    shutil.rmtree(buildtmp)

    logger.info("Generated module: " + module_name)
    logger.info("Output: " + output_binary)

    return output_binary


def compile_cxxcode(module_name, cxxcode, output_binary=None, keep_temp=False,
                    **kwargs):
    '''c++ code (string) -> temporary file -> native module.
    Returns the generated .so.

    '''

    # Get a temporary C++ file to compile
    fdpath = _write_temp(cxxcode, '.cpp')
    output_binary = compile_cxxfile(module_name, fdpath,
                                    output_binary, **kwargs)
    if not keep_temp:
        # remove tempfile
        os.remove(fdpath)
    else:
        logger.warning("Keeping temporary generated file:" + fdpath)

    return output_binary


def compile_pythrancode(module_name, pythrancode, specs=None,
                        opts=None, cpponly=False, pyonly=False,
                        output_file=None, module_dir=None, report_times=False, **kwargs):
    '''Pythran code (string) -> c++ code -> native module

    if `cpponly` is set to true, return the generated C++ filename
    if `pyonly` is set to true, prints the generated Python filename,
       unless `output_file` is set
    otherwise, return the generated native library filename
    '''

    if pyonly:
        # Only generate the optimized python code
        content = generate_py(module_name, pythrancode, opts, module_dir, report_times)
        if output_file is None:
            print(content)
            return None
        else:
            tmp_file = _write_temp(content, '.py')
            output_file = output_file.format('.py')
            shutil.move(tmp_file, output_file)
            logger.info("Generated Python source file: " + output_file)

    # Autodetect the Pythran spec if not given as parameter
    from pythran.spec import spec_parser
    if specs is None:
        specs = spec_parser(pythrancode)

    # Generate C++, get a PythonModule object
    module, error_checker = generate_cxx(module_name, pythrancode, specs, opts,
                                         module_dir, report_times)

    if 'ENABLE_PYTHON_MODULE' in kwargs.get('undef_macros', []):
        module.preamble.insert(0, Line('#undef ENABLE_PYTHON_MODULE'))
        module.preamble.insert(0, Line('#define PY_MAJOR_VERSION {}'.
                                       format(sys.version_info.major)))

    if cpponly:
        # User wants only the C++ code
        tmp_file = _write_temp(str(module), '.cpp')
        if output_file:
            output_file = output_file.replace('%{ext}', '.cpp')
        else:
            output_file = module_name + ".cpp"
        shutil.move(tmp_file, output_file)
        logger.info("Generated C++ source file: " + output_file)
    else:
        if not specs:
            raise ValueError("Empty spec files while generating native module")

        # Compile to binary
        try:
            output_file = compile_cxxcode(module_name,
                                          str(module),
                                          output_binary=output_file,
                                          **kwargs)
        except PythranCompileError:
            logger.warning("Compilation error, "
                           "trying hard to find its origin...")
            error_checker()
            logger.warning("Nope, I'm going to flood you with C++ errors!")
            raise

    return output_file


def import_pythrancode(pythrancode, **kwargs):
    # It should be safe to delete the library once it's loaded in memory.
    digest = hashlib.sha256(pythrancode.encode()).hexdigest()
    module_name = "pythranized_{}".format(digest)
    tmpfile = None
    try:
        tmpfile = compile_pythrancode(module_name, pythrancode, **kwargs)
        spec = importlib.util.spec_from_file_location(module_name, tmpfile)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        if tmpfile is not None:
            os.remove(tmpfile)


def compile_pythranfile(file_path, output_file=None, module_name=None,
                        cpponly=False, pyonly=False, report_times=False, **kwargs):
    """
    Pythran file -> c++ file -> native module.

    Returns the generated .so (or .cpp if `cpponly` is set to true).

    Usage without an existing spec file

    >>> with open('pythran_test.py', 'w') as fd:
    ...    _ = fd.write('def foo(i): return i ** 2')
    >>> cpp_path = compile_pythranfile('pythran_test.py', cpponly=True)

    Usage with an existing spec file:

    >>> with open('pythran_test.pythran', 'w') as fd:
    ...    _ = fd.write('export foo(int)')
    >>> so_path = compile_pythranfile('pythran_test.py')

    Specify the output file:

    >>> import sysconfig
    >>> ext = sysconfig.get_config_vars()["EXT_SUFFIX"]
    >>> so_path = compile_pythranfile('pythran_test.py', output_file='foo'+ext)
    """
    if not output_file:
        # derive module name from input file name
        _, basename = os.path.split(file_path)
        module_name = module_name or os.path.splitext(basename)[0]

    else:
        # derive module name from destination output_file name
        _, basename = os.path.split(output_file.replace('%{ext}', ''))
        module_name = module_name or basename.split(".", 1)[0]

    module_dir = os.path.dirname(file_path)

    # Look for an extra spec file
    spec_file = os.path.splitext(file_path)[0] + '.pythran'
    if os.path.isfile(spec_file):
        specs = load_specfile(spec_file)
        kwargs.setdefault('specs', specs)

    try:
        with open(file_path) as fd:
            output_file = compile_pythrancode(module_name, fd.read(),
                                              output_file=output_file,
                                              cpponly=cpponly, pyonly=pyonly,
                                              module_dir=module_dir,
                                              report_times=report_times,
                                              **kwargs)
    except PythranSyntaxError as e:
        if e.filename is None:
            e.filename = file_path
        raise

    return output_file


def import_pythranfile(pythranpath, **kwargs):
    with open(pythranpath) as fd:
        return import_pythrancode(fd.read(), **kwargs)


def test_compile():
    '''Simple passthrough compile test.
    May raises PythranCompileError Exception.

    '''
    code = '''
        #include <pythonic/core.hpp>
    '''

    output_file = compile_cxxcode('test', code)
    output_file and os.remove(output_file)
