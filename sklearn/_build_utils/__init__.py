"""
Utilities useful during the build.
"""
# author: Andy Mueller, Gael Varoquaux
# license: BSD


import os
import os.path as op
import sklearn
import contextlib
import shutil
from glob import glob
import textwrap

from distutils.version import LooseVersion

from .pre_build_helpers import basic_check_build
from .openmp_helpers import check_openmp_support
from .._min_dependencies import CYTHON_MIN_VERSION


DEFAULT_ROOT = 'sklearn'

VCOMP140_SRC_PATH = "C:\\Windows\System32\\vcomp140.dll"  # noqa
GLOB_PATTERN = "lib.*/sklearn"


def _check_cython_version():
    message = ('Please install Cython with a version >= {0} in order '
               'to build a scikit-learn from source.').format(
                    CYTHON_MIN_VERSION)
    try:
        import Cython
    except ModuleNotFoundError as e:
        # Re-raise with more informative error message instead:
        raise ModuleNotFoundError(message) from e

    if LooseVersion(Cython.__version__) < CYTHON_MIN_VERSION:
        message += (' The current version of Cython is {} installed in {}.'
                    .format(Cython.__version__, Cython.__path__))
        raise ValueError(message)


def cythonize_extensions(top_path, config):
    """Check that a recent Cython is available and cythonize extensions"""
    _check_cython_version()
    from Cython.Build import cythonize

    # Fast fail before cythonization if compiler fails compiling basic test
    # code even without OpenMP
    basic_check_build()

    # check simple compilation with OpenMP. If it fails scikit-learn will be
    # built without OpenMP and the test test_openmp_supported in the test suite
    # will fail.
    # `check_openmp_support` compiles a small test program to see if the
    # compilers are properly configured to build with OpenMP. This is expensive
    # and we only want to call this function once.
    # The result of this check is cached as a private attribute on the sklearn
    # module (only at build-time) to be used twice:
    # - First to set the value of SKLEARN_OPENMP_PARALLELISM_ENABLED, the
    #   cython build-time variable passed to the cythonize() call.
    # - Then in the build_ext subclass defined in the top-level setup.py file
    #   to actually build the compiled extensions with OpenMP flags if needed.
    sklearn._OPENMP_SUPPORTED = check_openmp_support()

    n_jobs = 1
    with contextlib.suppress(ImportError):
        import joblib
        if LooseVersion(joblib.__version__) > LooseVersion("0.13.0"):
            # earlier joblib versions don't account for CPU affinity
            # constraints, and may over-estimate the number of available
            # CPU particularly in CI (cf loky#114)
            n_jobs = joblib.cpu_count()

    config.ext_modules = cythonize(
        config.ext_modules,
        nthreads=n_jobs,
        compile_time_env={
            'SKLEARN_OPENMP_PARALLELISM_ENABLED': sklearn._OPENMP_SUPPORTED},
        compiler_directives={'language_level': 3})


def gen_from_templates(templates, top_path):
    """Generate cython files from a list of templates"""
    # Lazy import because cython is not a runtime dependency.
    from Cython import Tempita

    for template in templates:
        outfile = template.replace('.tp', '')

        # if the template is not updated, no need to output the cython file
        if not (os.path.exists(outfile) and
                os.stat(template).st_mtime < os.stat(outfile).st_mtime):

            with open(template, "r") as f:
                tmpl = f.read()

            tmpl_ = Tempita.sub(tmpl)

            with open(outfile, "w") as f:
                f.write(tmpl_)


def make_distributor_init(sklearn_dirname, dll_filename):
    """Create a _distributor_init.py file for the vcomp dll.

    This file is imported first when importing the sklearn package so as
    to pre-load the vendored vcomp dll.
    """
    distributor_init = op.join(sklearn_dirname, '_distributor_init.py')
    with open(distributor_init, 'wt') as f:
        f.write(textwrap.dedent("""
            '''
            Helper to preload the OpenMP dll to prevent "dll not found"
            errors.
            Once a DLL is preloaded, its namespace is made available to any
            subsequent DLL. This file originated in the scikit-learn-wheels
            github repo, and is created as part of the scripts that build the
            wheel.
            '''
            import os
            import os.path as op
            from ctypes import WinDLL


            if os.name == 'nt':
                # Pre-load the DLL stored in sklearn/.libs by convention.
                dll_path = op.join(op.dirname(__file__), '.libs', '{}')
                WinDLL(op.abspath(dll_path))

        """.format(dll_filename)))
    return op.abspath(distributor_init)


def embed_vcomp140(build_dirname):
    # TODO: use threadpoolctl to dynamically locate the right vcomp dll
    # instead? This would require first in-place building scikit-learn
    # to make it "importable".
    if not op.exists(VCOMP140_SRC_PATH):
        raise ValueError(f"Could not find {VCOMP140_SRC_PATH}.")

    if not op.isdir(build_dirname):
        raise RuntimeError(f"Could not find {build_dirname} folder. "
                           "Run 'python setup.py build' first.")

    target_folder_glob_pattern = op.join(build_dirname, GLOB_PATTERN)
    target_folders = glob(target_folder_glob_pattern)

    if len(target_folders) == 0:
        raise RuntimeError(f"Could not find folder matching "
                           f"{target_folder_glob_pattern}.")

    if len(target_folders) > 1:
        raise RuntimeError(f"Found too many target folders: "
                           f"{', '.join(target_folders)}.")

    target_folder = op.abspath(op.join(target_folders[0], ".libs"))

    # create the "sklearn/.libs" subfolder
    if not op.exists(target_folder):
        os.mkdir(target_folder)

    print(f"Copying {VCOMP140_SRC_PATH} to: \n{target_folder}.")
    shutil.copy2(VCOMP140_SRC_PATH, target_folder)

    sklearn_dirname = op.join(build_dirname, "..", "sklearn")
    dll_filename = op.basename(VCOMP140_SRC_PATH)

    # Generate the _distributor_init file in the source tree
    print("Generating the '_distributor_init.py' file in: ")
    print(make_distributor_init(sklearn_dirname, dll_filename))
