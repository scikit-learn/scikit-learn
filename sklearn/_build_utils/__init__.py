"""
Utilities useful during the build.
"""
# author: Andy Mueller, Gael Varoquaux
# license: BSD


import os
import sklearn
import contextlib

from distutils.version import LooseVersion

from .pre_build_helpers import basic_check_build
from .openmp_helpers import check_openmp_support


DEFAULT_ROOT = 'sklearn'
# on conda, this is the latest for python 3.5
CYTHON_MIN_VERSION = '0.28.5'


def build_from_c_and_cpp_files(extensions):
    """Modify the extensions to build from the .c and .cpp files.

    This is useful for releases, this way cython is not required to
    run python setup.py install.
    """
    for extension in extensions:
        sources = []
        for sfile in extension.sources:
            path, ext = os.path.splitext(sfile)
            if ext in ('.pyx', '.py'):
                if extension.language == 'c++':
                    ext = '.cpp'
                else:
                    ext = '.c'
                sfile = path + ext
            sources.append(sfile)
        extension.sources = sources


def maybe_cythonize_extensions(top_path, config):
    """Tweaks for building extensions between release and development mode."""
    # Fast fail before cythonization if compiler fails compiling basic test
    # code even without OpenMP
    basic_check_build()

    # check simple compilation with OpenMP. If it fails scikit-learn will be
    # built without OpenMP and scikit-learn will warn the user when imported.
    # `check_openmp_support` compiles a small test program to see if the
    # compilers are properly configured to build with OpenMP. This is expensive
    # and we only want to call this function once.
    # The result of this check is cached as a private attribute on the sklearn
    # module (only at build-time) to be used twice:
    # - First to set the value of SKLEARN_OPENMP_SUPPORTED, the cython
    #   build-time variable passed to the cythonize() call.
    # - Then in the build_ext subclass defined in the top-level setup.py file
    #   to actually build the compiled extensions with OpenMP flags if needed.
    sklearn._OPENMP_SUPPORTED = check_openmp_support()

    is_release = os.path.exists(os.path.join(top_path, 'PKG-INFO'))

    if is_release:
        build_from_c_and_cpp_files(config.ext_modules)
    else:
        message = ('Please install cython with a version >= {0} in order '
                   'to build a scikit-learn development version.').format(
                       CYTHON_MIN_VERSION)
        try:
            import Cython
            if LooseVersion(Cython.__version__) < CYTHON_MIN_VERSION:
                message += ' Your version of Cython was {0}.'.format(
                    Cython.__version__)
                raise ValueError(message)
            from Cython.Build import cythonize
        except ImportError as exc:
            exc.args += (message,)
            raise

        n_jobs = 1
        with contextlib.suppress(ImportError):
            import joblib
            if LooseVersion(joblib.__version__) > LooseVersion("0.13.0"):
                # earlier joblib versions don't account for CPU affinity
                # constraints, and may over-estimate the number of available
                # CPU particularly in CI (cf loky#114)
                n_jobs = joblib.effective_n_jobs()

        config.ext_modules = cythonize(
            config.ext_modules,
            nthreads=n_jobs,
            compile_time_env={
                'SKLEARN_OPENMP_SUPPORTED': sklearn._OPENMP_SUPPORTED},
            compiler_directives={'language_level': 3})


def gen_from_templates(templates, top_path):
    """Generate cython files from a list of templates"""
    is_release = os.path.exists(os.path.join(top_path, 'PKG-INFO'))
    # Files are already cythonized, nothing to do.
    if is_release:
        return

    # Lazy import because cython is not a dependency when building from
    # source distribution.
    from Cython import Tempita # noqa

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
