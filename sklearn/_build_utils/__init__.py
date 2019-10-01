"""
Utilities useful during the build.
"""
# author: Andy Mueller, Gael Varoquaux
# license: BSD


import os

from distutils.version import LooseVersion
import contextlib

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
    with_openmp = check_openmp_support()

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
            compile_time_env={'SKLEARN_OPENMP_SUPPORTED': with_openmp},
            compiler_directives={'language_level': 3})
