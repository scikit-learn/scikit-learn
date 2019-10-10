"""Helpers for OpenMP support during the build."""

# This code is adapted for a large part from the astropy openmp helpers, which
# can be found at: https://github.com/astropy/astropy-helpers/blob/master/astropy_helpers/openmp_helpers.py  # noqa


import os
import sys
import glob
import tempfile
import textwrap
import warnings
import subprocess

from numpy.distutils.ccompiler import new_compiler
from distutils.sysconfig import customize_compiler
from distutils.errors import CompileError, LinkError


CCODE = textwrap.dedent(
    """\
    #include <stdio.h>
    int main(void) {
    printf("success\\n");
    return 0;
    }
    """)


CCODE_OPENMP = textwrap.dedent(
    """\
    #include <omp.h>
    #include <stdio.h>
    int main(void) {
    #pragma omp parallel
    printf("nthreads=%d\\n", omp_get_num_threads());
    return 0;
    }
    """)


def compile_test_program(code, extra_preargs=[], extra_postargs=[]):
    """Check that some C code can be compiled and run"""
    ccompiler = new_compiler()
    customize_compiler(ccompiler)

    # extra_(pre/post)args can be a callable to make it possible to get it's
    # value from the compiler
    if callable(extra_preargs):
        extra_preargs = extra_preargs(ccompiler)
    if callable(extra_postargs):
        extra_postargs = extra_postargs(ccompiler)

    start_dir = os.path.abspath('.')

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.chdir(tmp_dir)

        # Write test program
        with open('test_program.c', 'w') as f:
            f.write(code)

        os.mkdir('objects')

        # Compile, test program
        ccompiler.compile(['test_program.c'], output_dir='objects',
                          extra_postargs=extra_postargs)

        # Link test program
        objects = glob.glob(
            os.path.join('objects', '*' + ccompiler.obj_extension))
        ccompiler.link_executable(objects, 'test_program',
                                  extra_preargs=extra_preargs,
                                  extra_postargs=extra_postargs)

        # Run test program
        output = subprocess.check_output('./test_program')
        output = output.decode(sys.stdout.encoding or 'utf-8').splitlines()

        os.chdir(start_dir)

    return output


def basic_check_build():
    """Check basic compilation and linking of C code"""
    try:
        output = compile_test_program(CCODE)
        if "success" in output[0]:
            compiler_ok = True
        else:
            compiler_ok = False
    except (CompileError, LinkError, subprocess.CalledProcessError):
        compiler_ok = False

    return compiler_ok


def get_openmp_flag(compiler):
    if hasattr(compiler, 'compiler'):
        compiler = compiler.compiler[0]
    else:
        compiler = compiler.__class__.__name__

    if sys.platform == "win32" and ('icc' in compiler or 'icl' in compiler):
        return ['/Qopenmp']
    elif sys.platform == "win32":
        return ['/openmp']
    elif sys.platform == "darwin" and ('icc' in compiler or 'icl' in compiler):
        return ['-openmp']
    elif sys.platform == "darwin" and 'openmp' in os.getenv('CPPFLAGS', ''):
        # -fopenmp can't be passed as compile flag when using Apple-clang.
        # OpenMP support has to be enabled during preprocessing.
        #
        # For example, our macOS wheel build jobs use the following environment
        # variables to build with Apple-clang and the brew installed "libomp":
        #
        # export CPPFLAGS="$CPPFLAGS -Xpreprocessor -fopenmp"
        # export CFLAGS="$CFLAGS -I/usr/local/opt/libomp/include"
        # export CXXFLAGS="$CXXFLAGS -I/usr/local/opt/libomp/include"
        # export LDFLAGS="$LDFLAGS -L/usr/local/opt/libomp/lib -lomp"
        # export DYLD_LIBRARY_PATH=/usr/local/opt/libomp/lib
        return []
    # Default flag for GCC and clang:
    return ['-fopenmp']


def check_openmp_support():
    """Check whether OpenMP test code can be compiled and run"""
    if os.getenv('SKLEARN_NO_OPENMP'):
        # Build explicitly without OpenMP support
        return "explicitly disabled"

    if not basic_check_build():
        # simple build without OpenMP fails. No need to check OpenMP (avoids
        # possible wrong interpretation of error message)
        return "unrelated fail"

    extra_preargs = os.getenv('LDFLAGS', None)
    if extra_preargs is not None:
        extra_preargs = extra_preargs.split(" ")
    else:
        extra_preargs = []

    extra_postargs = get_openmp_flag

    try:
        output = compile_test_program(CCODE_OPENMP,
                                      extra_preargs=extra_preargs,
                                      extra_postargs=extra_postargs)

        if 'nthreads=' in output[0]:
            nthreads = int(output[0].strip().split('=')[1])
            openmp_supported = (len(output) == nthreads)
        else:
            openmp_supported = False

    except (CompileError, LinkError, subprocess.CalledProcessError):
        openmp_supported = False

    message = textwrap.dedent(
        """
                            ***

        It seems that scikit-learn cannot be built with OpenMP support.

        - Make sure you have followed the installation instructions:

            https://scikit-learn.org/dev/developers/advanced_installation.html

        - If your compiler supports OpenMP but the build still fails, please
          submit a bug report at:

            https://github.com/scikit-learn/scikit-learn/issues

        - If you want to build scikit-learn without OpenMP support, you can set
          the environment variable SKLEARN_NO_OPENMP and rerun the build
          command. Note however that some estimators will run in sequential
          mode and their `n_jobs` parameter will have no effect anymore.

                            ***
        """)

    if not openmp_supported:
        warnings.warns(message)
        return "unsupported"

    return "supported"
