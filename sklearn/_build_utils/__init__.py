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

# The following places need to be in sync with regard to Cython version:
# - .circleci config file
# - sklearn/_build_utils/__init__.py
# - advanced installation guide
CYTHON_MIN_VERSION = '0.28.5'


def _check_cython_version():
    message = ('Please install Cython with a version >= {0} in order '
               'to build a scikit-learn from source.').format(
                    CYTHON_MIN_VERSION)
    try:
        import Cython
    except ModuleNotFoundError:
        # Re-raise with more informative error message instead:
        raise ModuleNotFoundError(message)

    if LooseVersion(Cython.__version__) < CYTHON_MIN_VERSION:
        message += (' The current version of Cython is {} installed in {}.'
                    .format(Cython.__version__, Cython.__path__))
        raise ValueError(message)


def cythonize_extensions(top_path, config):
    """Check that a recent Cython is available and cythonize extensions"""
    _check_cython_version()
    from Cython.Build import cythonize

    with_openmp = check_openmp_support()
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
        compile_time_env={'SKLEARN_OPENMP_SUPPORTED': with_openmp},
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
