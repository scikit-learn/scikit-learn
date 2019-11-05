"""
Utilities useful during the build.
"""
# author: Andy Mueller, Gael Varoquaux
# license: BSD


import os
from os.path import join, dirname
from distutils.version import LooseVersion
import contextlib

from .openmp_helpers import check_openmp_support


DEFAULT_ROOT = 'sklearn'


def get_cython_version():
    # Try to find pyproject.toml
    pyproject_toml = join(dirname(__file__), '..', '..', 'pyproject.toml')
    if not os.path.exists(pyproject_toml):
        raise ImportError("Unable to locate file pyproject.toml")

    # Try to find the minimum version from pyproject.toml
    with open(pyproject_toml) as pt:
        for line in pt:
            if "cython" not in line.lower():
                continue
            _, line = line.split('=')
            required_version, _ = line.split('"')
            return required_version
        else:
            raise KeyError("Can't find cython min version in pyproject.toml")


CYTHON_MIN_VERSION = get_cython_version()


def cythonize_extensions(top_path, config):
    """Tweaks for building extensions between release and development mode."""
    with_openmp = check_openmp_support()

    # Could be due to a too old pip version and build isolation, check that
    try:
        # Note, pip may not be installed or not have been used
        import pip
        if LooseVersion(pip.__version__) < LooseVersion('18.0.0'):
            raise RuntimeError("Cannot cythonize, possibly due "
                               "to `pip` being too old, found version {}, "
                               "needed is >= 18.0.0.".format(
                                   pip.__version__))
    except ImportError:
        raise RuntimeError("pip not installed!")

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
