"""
=====================
Cython related magics
=====================

Magic command interface for interactive work with Cython

.. note::

  The ``Cython`` package needs to be installed separately. It
  can be obtained using ``easy_install`` or ``pip``.

Usage
=====

To enable the magics below, execute ``%load_ext cython``.

``%%cython``

{CYTHON_DOC}

``%%cython_inline``

{CYTHON_INLINE_DOC}

``%%cython_pyximport``

{CYTHON_PYXIMPORT_DOC}

Author:
* Brian Granger

Code moved from IPython and adapted by:
* Martín Gaitán

Parts of this code were taken from Cython.inline.
"""
#-----------------------------------------------------------------------------
# Copyright (C) 2010-2011, IPython Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file ipython-COPYING.rst, distributed with this software.
#-----------------------------------------------------------------------------


import io
import os
import re
import sys
import time
import copy
import distutils.log
import textwrap

IO_ENCODING = sys.getfilesystemencoding()

import hashlib
from distutils.core import Distribution, Extension
from distutils.command.build_ext import build_ext

from IPython.core import display
from IPython.core import magic_arguments
from IPython.core.magic import Magics, magics_class, cell_magic
try:
    from IPython.paths import get_ipython_cache_dir
except ImportError:
    # older IPython version
    from IPython.utils.path import get_ipython_cache_dir
from IPython.utils.text import dedent

from ..Shadow import __version__ as cython_version
from ..Compiler.Errors import CompileError
from .Inline import cython_inline, load_dynamic
from .Dependencies import cythonize
from ..Utils import captured_fd, print_captured


PGO_CONFIG = {
    'gcc': {
        'gen': ['-fprofile-generate', '-fprofile-dir={TEMPDIR}'],
        'use': ['-fprofile-use', '-fprofile-correction', '-fprofile-dir={TEMPDIR}'],
    },
    # blind copy from 'configure' script in CPython 3.7
    'icc': {
        'gen': ['-prof-gen'],
        'use': ['-prof-use'],
    }
}
PGO_CONFIG['mingw32'] = PGO_CONFIG['gcc']


@magics_class
class CythonMagics(Magics):

    def __init__(self, shell):
        super().__init__(shell)
        self._reloads = {}
        self._code_cache = {}
        self._pyximport_installed = False

    def _import_all(self, module):
        mdict = module.__dict__
        if '__all__' in mdict:
            keys = mdict['__all__']
        else:
            keys = [k for k in mdict if not k.startswith('_')]

        for k in keys:
            try:
                self.shell.push({k: mdict[k]})
            except KeyError:
                msg = "'module' object has no attribute '%s'" % k
                raise AttributeError(msg)

    @cell_magic
    def cython_inline(self, line, cell):
        """Compile and run a Cython code cell using Cython.inline.

        This magic simply passes the body of the cell to Cython.inline
        and returns the result. If the variables `a` and `b` are defined
        in the user's namespace, here is a simple example that returns
        their sum::

            %%cython_inline
            return a+b

        For most purposes, we recommend the usage of the `%%cython` magic.
        """
        locs = self.shell.user_global_ns
        globs = self.shell.user_ns
        return cython_inline(cell, locals=locs, globals=globs)

    @cell_magic
    def cython_pyximport(self, line, cell):
        """Compile and import a Cython code cell using pyximport.

        The contents of the cell are written to a `.pyx` file in the current
        working directory, which is then imported using `pyximport`. This
        magic requires a module name to be passed::

            %%cython_pyximport modulename
            def f(x):
                return 2.0*x

        The compiled module is then imported and all of its symbols are
        injected into the user's namespace. For most purposes, we recommend
        the usage of the `%%cython` magic.
        """
        module_name = line.strip()
        if not module_name:
            raise ValueError('module name must be given')
        fname = module_name + '.pyx'
        with open(fname, 'w', encoding='utf-8') as f:
            f.write(cell)
        if 'pyximport' not in sys.modules or not self._pyximport_installed:
            import pyximport
            pyximport.install()
            self._pyximport_installed = True
        if module_name in self._reloads:
            module = self._reloads[module_name]
            # Note: reloading extension modules is not actually supported
            # (requires PEP-489 reinitialisation support).
            # Don't know why this should ever have worked as it reads here.
            # All we really need to do is to update the globals below.
            #reload(module)
        else:
            __import__(module_name)
            module = sys.modules[module_name]
            self._reloads[module_name] = module
        self._import_all(module)

    @magic_arguments.magic_arguments()
    @magic_arguments.argument(
        '-a', '--annotate', action='store_const', const='default', dest='annotate',
        help="Produce a colorized HTML version of the source."
    )
    @magic_arguments.argument(
        '--annotate-fullc', action='store_const', const='fullc', dest='annotate',
        help="Produce a colorized HTML version of the source "
             "which includes entire generated C/C++-code."
    )
    @magic_arguments.argument(
        '-+', '--cplus', action='store_true', default=False,
        help="Output a C++ rather than C file."
    )
    @magic_arguments.argument(
        '-3', dest='language_level', action='store_const', const=3, default=None,
        help="Select Python 3 syntax."
    )
    @magic_arguments.argument(
        '-2', dest='language_level', action='store_const', const=2, default=None,
        help="Select Python 2 syntax."
    )
    @magic_arguments.argument(
        '-f', '--force', action='store_true', default=False,
        help="Force the compilation of a new module, even if the source has been "
             "previously compiled."
    )
    @magic_arguments.argument(
        '-c', '--compile-args', action='append', default=[],
        help="Extra flags to pass to compiler via the `extra_compile_args` "
             "Extension flag (can be specified  multiple times)."
    )
    @magic_arguments.argument(
        '--link-args', action='append', default=[],
        help="Extra flags to pass to linker via the `extra_link_args` "
             "Extension flag (can be specified  multiple times)."
    )
    @magic_arguments.argument(
        '-l', '--lib', action='append', default=[],
        help="Add a library to link the extension against (can be specified "
             "multiple times)."
    )
    @magic_arguments.argument(
        '-n', '--name',
        help="Specify a name for the Cython module."
    )
    @magic_arguments.argument(
        '-L', dest='library_dirs', metavar='dir', action='append', default=[],
        help="Add a path to the list of library directories (can be specified "
             "multiple times)."
    )
    @magic_arguments.argument(
        '-I', '--include', action='append', default=[],
        help="Add a path to the list of include directories (can be specified "
             "multiple times)."
    )
    @magic_arguments.argument(
        '-S', '--src', action='append', default=[],
        help="Add a path to the list of src files (can be specified "
             "multiple times)."
    )
    @magic_arguments.argument(
        '--pgo', dest='pgo', action='store_true', default=False,
        help=("Enable profile guided optimisation in the C compiler. "
              "Compiles the cell twice and executes it in between to generate a runtime profile.")
    )
    @magic_arguments.argument(
        '--verbose', dest='quiet', action='store_false', default=True,
        help=("Print debug information like generated .c/.cpp file location "
              "and exact gcc/g++ command invoked.")
    )
    @cell_magic
    def cython(self, line, cell):
        """Compile and import everything from a Cython code cell.

        The contents of the cell are written to a `.pyx` file in the
        directory returned by `get_ipython_cache_dir()/cython` using a filename
        with the hash of the code. This file is then cythonized and compiled.
        The resulting module is imported and all of its symbols are injected
        into the user's namespace. The usage is similar to that of
        `%%cython_pyximport` but you don't have to pass a module name::

            %%cython
            def f(x):
                return 2.0*x

        To compile OpenMP codes, pass the required  `--compile-args`
        and `--link-args`.  For example with gcc::

            %%cython --compile-args=-fopenmp --link-args=-fopenmp
            ...

        To enable profile guided optimisation, pass the ``--pgo`` option.
        Note that the cell itself needs to take care of establishing a suitable
        profile when executed. This can be done by implementing the functions to
        optimise, and then calling them directly in the same cell on some realistic
        training data like this::

            %%cython --pgo
            def critical_function(data):
                for item in data:
                    ...

            # execute function several times to build profile
            from somewhere import some_typical_data
            for _ in range(100):
                critical_function(some_typical_data)

        In Python 3.5 and later, you can distinguish between the profile and
        non-profile runs as follows::

            if "_pgo_" in __name__:
                ...  # execute critical code here
        """
        args = magic_arguments.parse_argstring(self.cython, line)
        code = cell if cell.endswith('\n') else cell + '\n'
        lib_dir = os.path.join(get_ipython_cache_dir(), 'cython')
        key = (code, line, sys.version_info, sys.executable, cython_version)

        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)

        if args.pgo:
            key += ('pgo',)
        if args.force:
            # Force a new module name by adding the current time to the
            # key which is hashed to determine the module name.
            key += (time.time(),)

        if args.name:
            module_name = str(args.name)  # no-op in Py3
        else:
            module_name = "_cython_magic_" + hashlib.sha256(str(key).encode('utf-8')).hexdigest()
        html_file = os.path.join(lib_dir, module_name + '.html')
        module_path = os.path.join(lib_dir, module_name + self.so_ext)

        have_module = os.path.isfile(module_path)
        need_cythonize = args.pgo or not have_module

        if args.annotate:
            if not os.path.isfile(html_file):
                need_cythonize = True

        extension = None
        if need_cythonize:
            extensions = self._cythonize(module_name, code, lib_dir, args, quiet=args.quiet)
            if extensions is None:
                # Compilation failed and printed error message
                return None
            assert len(extensions) == 1
            extension = extensions[0]
            self._code_cache[key] = module_name

            if args.pgo:
                self._profile_pgo_wrapper(extension, lib_dir)

        def print_compiler_output(stdout, stderr, where):
            # On windows, errors are printed to stdout, we redirect both to sys.stderr.
            print_captured(stdout, where, "Content of stdout:\n")
            print_captured(stderr, where, "Content of stderr:\n")

        get_stderr = get_stdout = None
        try:
            with captured_fd(1) as get_stdout:
                with captured_fd(2) as get_stderr:
                    self._build_extension(
                        extension, lib_dir, pgo_step_name='use' if args.pgo else None, quiet=args.quiet)
        except (distutils.errors.CompileError, distutils.errors.LinkError):
            # Build failed, print error message from compiler/linker
            print_compiler_output(get_stdout(), get_stderr(), sys.stderr)
            return None

        # Build seems ok, but we might still want to show any warnings that occurred
        print_compiler_output(get_stdout(), get_stderr(), sys.stdout)

        module = load_dynamic(module_name, module_path)
        self._import_all(module)

        if args.annotate:
            try:
                with open(html_file, encoding='utf-8') as f:
                    annotated_html = f.read()
            except OSError as e:
                # File could not be opened. Most likely the user has a version
                # of Cython before 0.15.1 (when `cythonize` learned the
                # `force` keyword argument) and has already compiled this
                # exact source without annotation.
                print('Cython completed successfully but the annotated '
                      'source could not be read.', file=sys.stderr)
                print(e, file=sys.stderr)
            else:
                return display.HTML(self.clean_annotated_html(annotated_html))

    def _profile_pgo_wrapper(self, extension, lib_dir):
        """
        Generate a .c file for a separate extension module that calls the
        module init function of the original module.  This makes sure that the
        PGO profiler sees the correct .o file of the final module, but it still
        allows us to import the module under a different name for profiling,
        before recompiling it into the PGO optimised module.  Overwriting and
        reimporting the same shared library is not portable.
        """
        extension = copy.copy(extension)  # shallow copy, do not modify sources in place!
        module_name = extension.name
        pgo_module_name = '_pgo_' + module_name
        pgo_wrapper_c_file = os.path.join(lib_dir, pgo_module_name + '.c')
        with open(pgo_wrapper_c_file, 'w', encoding='utf-8') as f:
            f.write(textwrap.dedent("""
            #include "Python.h"
            extern PyMODINIT_FUNC PyInit_%(module_name)s(void);
            PyMODINIT_FUNC PyInit_%(pgo_module_name)s(void); /*proto*/
            PyMODINIT_FUNC PyInit_%(pgo_module_name)s(void) {
                return PyInit_%(module_name)s();
            }
            """ % {'module_name': module_name, 'pgo_module_name': pgo_module_name}))

        extension.sources = extension.sources + [pgo_wrapper_c_file]  # do not modify in place!
        extension.name = pgo_module_name

        self._build_extension(extension, lib_dir, pgo_step_name='gen')

        # import and execute module code to generate profile
        so_module_path = os.path.join(lib_dir, pgo_module_name + self.so_ext)
        load_dynamic(pgo_module_name, so_module_path)

    def _cythonize(self, module_name, code, lib_dir, args, quiet=True):
        pyx_file = os.path.join(lib_dir, module_name + '.pyx')

        c_include_dirs = args.include
        c_src_files = list(map(str, args.src))
        if 'numpy' in code:
            import numpy
            c_include_dirs.append(numpy.get_include())
        with open(pyx_file, 'w', encoding='utf-8') as f:
            f.write(code)
        extension = Extension(
            name=module_name,
            sources=[pyx_file] + c_src_files,
            include_dirs=c_include_dirs,
            library_dirs=args.library_dirs,
            extra_compile_args=args.compile_args,
            extra_link_args=args.link_args,
            libraries=args.lib,
            language='c++' if args.cplus else 'c',
        )
        try:
            opts = dict(
                quiet=quiet,
                annotate=args.annotate,
                force=True,
                language_level=min(3, sys.version_info[0]),
            )
            if args.language_level is not None:
                assert args.language_level in (2, 3)
                opts['language_level'] = args.language_level
            return cythonize([extension], **opts)
        except CompileError:
            return None

    def _build_extension(self, extension, lib_dir, temp_dir=None, pgo_step_name=None, quiet=True):
        build_extension = self._get_build_extension(
            extension, lib_dir=lib_dir, temp_dir=temp_dir, pgo_step_name=pgo_step_name)
        old_threshold = None
        try:
            if not quiet:
                old_threshold = distutils.log.set_threshold(distutils.log.DEBUG)
            build_extension.run()
        finally:
            if not quiet and old_threshold is not None:
                distutils.log.set_threshold(old_threshold)

    def _add_pgo_flags(self, build_extension, step_name, temp_dir):
        compiler_type = build_extension.compiler.compiler_type
        if compiler_type == 'unix':
            compiler_cmd = build_extension.compiler.compiler_so
            # TODO: we could try to call "[cmd] --version" for better insights
            if not compiler_cmd:
                pass
            elif 'clang' in compiler_cmd or 'clang' in compiler_cmd[0]:
                compiler_type = 'clang'
            elif 'icc' in compiler_cmd or 'icc' in compiler_cmd[0]:
                compiler_type = 'icc'
            elif 'gcc' in compiler_cmd or 'gcc' in compiler_cmd[0]:
                compiler_type = 'gcc'
            elif 'g++' in compiler_cmd or 'g++' in compiler_cmd[0]:
                compiler_type = 'gcc'
        config = PGO_CONFIG.get(compiler_type)
        orig_flags = []
        if config and step_name in config:
            flags = [f.format(TEMPDIR=temp_dir) for f in config[step_name]]
            for extension in build_extension.extensions:
                orig_flags.append((extension.extra_compile_args, extension.extra_link_args))
                extension.extra_compile_args = extension.extra_compile_args + flags
                extension.extra_link_args = extension.extra_link_args + flags
        else:
            print("No PGO %s configuration known for C compiler type '%s'" % (step_name, compiler_type),
                  file=sys.stderr)
        return orig_flags

    @property
    def so_ext(self):
        """The extension suffix for compiled modules."""
        try:
            return self._so_ext
        except AttributeError:
            self._so_ext = self._get_build_extension().get_ext_filename('')
            return self._so_ext

    def _clear_distutils_mkpath_cache(self):
        """clear distutils mkpath cache

        prevents distutils from skipping re-creation of dirs that have been removed
        """
        try:
            from distutils.dir_util import _path_created
        except ImportError:
            pass
        else:
            _path_created.clear()

    def _get_build_extension(self, extension=None, lib_dir=None, temp_dir=None,
                             pgo_step_name=None, _build_ext=build_ext):
        self._clear_distutils_mkpath_cache()
        dist = Distribution()
        config_files = dist.find_config_files()
        try:
            config_files.remove('setup.cfg')
        except ValueError:
            pass
        dist.parse_config_files(config_files)

        if not temp_dir:
            temp_dir = lib_dir
        add_pgo_flags = self._add_pgo_flags

        if pgo_step_name:
            base_build_ext = _build_ext
            class _build_ext(_build_ext):
                def build_extensions(self):
                    add_pgo_flags(self, pgo_step_name, temp_dir)
                    base_build_ext.build_extensions(self)

        build_extension = _build_ext(dist)
        build_extension.finalize_options()
        if temp_dir:
            build_extension.build_temp = temp_dir
        if lib_dir:
            build_extension.build_lib = lib_dir
        if extension is not None:
            build_extension.extensions = [extension]
        return build_extension

    @staticmethod
    def clean_annotated_html(html, include_style=True):
        """Clean up the annotated HTML source.

        Strips the link to the generated C or C++ file, which we do not
        present to the user.

        Returns an HTML snippet (no <html>, <head>, or <body>),
        containing only the style tag(s) and _contents_ of the body,
        appropriate for embedding multiple times in cell output.
        """
        # extract CSS and body, rather than full HTML document
        chunks = []
        if include_style:
            styles = re.findall("<style.*</style>", html, re.MULTILINE | re.DOTALL)
            chunks.extend(styles)
        # extract body
        body = re.search(
            r"<body[^>]*>(.+)</body>", html, re.MULTILINE | re.DOTALL
        ).group(1)

        # exclude link to generated file
        r = re.compile('<p>Raw output: <a href="(.*)">(.*)</a>')
        for line in body.splitlines():
            if not r.match(line):
                chunks.append(line)
        return "\n".join(chunks)

__doc__ = __doc__.format(
    # rST doesn't see the -+ flag as part of an option list, so we
    # hide it from the module-level docstring.
    CYTHON_DOC=dedent(CythonMagics.cython.__doc__
                                  .replace('-+, --cplus', '--cplus    ')),
    CYTHON_INLINE_DOC=dedent(CythonMagics.cython_inline.__doc__),
    CYTHON_PYXIMPORT_DOC=dedent(CythonMagics.cython_pyximport.__doc__),
)
