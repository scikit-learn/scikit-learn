#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Cournapeau David <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
# License: 3-clause BSD

import sys
import os
import platform
import shutil

# We need to import setuptools before because it monkey-patches distutils
import setuptools  # noqa
from distutils.command.clean import clean as Clean
from distutils.command.sdist import sdist

import traceback
import importlib

try:
    import builtins
except ImportError:
    # Python 2 compat: just to be able to declare that Python >=3.8 is needed.
    import __builtin__ as builtins

# This is a bit (!) hackish: we are setting a global variable so that the
# main sklearn __init__ can detect if it is being loaded by the setup
# routine, to avoid attempting to load components that aren't built yet:
# the numpy distutils extensions that are used by scikit-learn to
# recursively build the compiled extensions in sub-packages is based on the
# Python import machinery.
builtins.__SKLEARN_SETUP__ = True


DISTNAME = "scikit-learn"
DESCRIPTION = "A set of python modules for machine learning and data mining"
with open("README.rst") as f:
    LONG_DESCRIPTION = f.read()
MAINTAINER = "Andreas Mueller"
MAINTAINER_EMAIL = "amueller@ais.uni-bonn.de"
URL = "http://scikit-learn.org"
DOWNLOAD_URL = "https://pypi.org/project/scikit-learn/#files"
LICENSE = "new BSD"
PROJECT_URLS = {
    "Bug Tracker": "https://github.com/scikit-learn/scikit-learn/issues",
    "Documentation": "https://scikit-learn.org/stable/documentation.html",
    "Source Code": "https://github.com/scikit-learn/scikit-learn",
}

# We can actually import a restricted version of sklearn that
# does not need the compiled code
import sklearn  # noqa
import sklearn._min_dependencies as min_deps  # noqa
from sklearn.externals._packaging.version import parse as parse_version  # noqa


VERSION = sklearn.__version__

# See: https://numpy.org/doc/stable/reference/c-api/deprecations.html
DEFINE_MACRO_NUMPY_C_API = (
    "NPY_NO_DEPRECATED_API",
    "NPY_1_7_API_VERSION",
)

# XXX: add new extensions to this list when they
# are not using the old NumPy C API (i.e. version 1.7)
# TODO: when Cython>=3.0 is used, make sure all Cython extensions
# use the newest NumPy C API by `#defining` `NPY_NO_DEPRECATED_API` to be
# `NPY_1_7_API_VERSION`, and remove this list.
# See: https://github.com/cython/cython/blob/1777f13461f971d064bd1644b02d92b350e6e7d1/docs/src/userguide/migrating_to_cy30.rst#numpy-c-api # noqa
USE_NEWEST_NUMPY_C_API = (
    "sklearn.__check_build._check_build",
    "sklearn._loss._loss",
    "sklearn.cluster._k_means_common",
    "sklearn.cluster._k_means_lloyd",
    "sklearn.cluster._k_means_elkan",
    "sklearn.cluster._k_means_minibatch",
    "sklearn.datasets._svmlight_format_fast",
    "sklearn.decomposition._cdnmf_fast",
    "sklearn.ensemble._hist_gradient_boosting._gradient_boosting",
    "sklearn.ensemble._hist_gradient_boosting.histogram",
    "sklearn.ensemble._hist_gradient_boosting.splitting",
    "sklearn.ensemble._hist_gradient_boosting._binning",
    "sklearn.ensemble._hist_gradient_boosting._predictor",
    "sklearn.ensemble._hist_gradient_boosting._bitset",
    "sklearn.ensemble._hist_gradient_boosting.common",
    "sklearn.ensemble._hist_gradient_boosting.utils",
    "sklearn.feature_extraction._hashing_fast",
    "sklearn.manifold._barnes_hut_tsne",
    "sklearn.metrics.cluster._expected_mutual_info_fast",
    "sklearn.metrics._pairwise_distances_reduction._datasets_pair",
    "sklearn.metrics._pairwise_distances_reduction._gemm_term_computer",
    "sklearn.metrics._pairwise_distances_reduction._base",
    "sklearn.metrics._pairwise_distances_reduction._argkmin",
    "sklearn.metrics._pairwise_distances_reduction._radius_neighborhood",
    "sklearn.metrics._pairwise_fast",
    "sklearn.neighbors._partition_nodes",
    "sklearn.tree._splitter",
    "sklearn.tree._utils",
    "sklearn.utils._cython_blas",
    "sklearn.utils._fast_dict",
    "sklearn.utils._openmp_helpers",
    "sklearn.utils._weight_vector",
    "sklearn.utils._random",
    "sklearn.utils._logistic_sigmoid",
    "sklearn.utils._readonly_array_wrapper",
    "sklearn.utils._typedefs",
    "sklearn.utils._heap",
    "sklearn.utils._sorting",
    "sklearn.utils._vector_sentinel",
    "sklearn.utils._isfinite",
    "sklearn.svm._newrand",
    "sklearn._isotonic",
)

# For some commands, use setuptools
SETUPTOOLS_COMMANDS = {
    "develop",
    "release",
    "bdist_egg",
    "bdist_rpm",
    "bdist_wininst",
    "install_egg_info",
    "build_sphinx",
    "egg_info",
    "easy_install",
    "upload",
    "bdist_wheel",
    "--single-version-externally-managed",
}
if SETUPTOOLS_COMMANDS.intersection(sys.argv):
    extra_setuptools_args = dict(
        zip_safe=False,  # the package can run out of an .egg file
        include_package_data=True,
        extras_require={
            key: min_deps.tag_to_packages[key]
            for key in ["examples", "docs", "tests", "benchmark"]
        },
    )
else:
    extra_setuptools_args = dict()


# Custom clean command to remove build artifacts


class CleanCommand(Clean):
    description = "Remove build artifacts from the source tree"

    def run(self):
        Clean.run(self)
        # Remove c files if we are not within a sdist package
        cwd = os.path.abspath(os.path.dirname(__file__))
        remove_c_files = not os.path.exists(os.path.join(cwd, "PKG-INFO"))
        if remove_c_files:
            print("Will remove generated .c files")
        if os.path.exists("build"):
            shutil.rmtree("build")
        for dirpath, dirnames, filenames in os.walk("sklearn"):
            for filename in filenames:
                if any(
                    filename.endswith(suffix)
                    for suffix in (".so", ".pyd", ".dll", ".pyc")
                ):
                    os.unlink(os.path.join(dirpath, filename))
                    continue
                extension = os.path.splitext(filename)[1]
                if remove_c_files and extension in [".c", ".cpp"]:
                    pyx_file = str.replace(filename, extension, ".pyx")
                    if os.path.exists(os.path.join(dirpath, pyx_file)):
                        os.unlink(os.path.join(dirpath, filename))
            for dirname in dirnames:
                if dirname == "__pycache__":
                    shutil.rmtree(os.path.join(dirpath, dirname))


cmdclass = {"clean": CleanCommand, "sdist": sdist}

# Custom build_ext command to set OpenMP compile flags depending on os and
# compiler. Also makes it possible to set the parallelism level via
# and environment variable (useful for the wheel building CI).
# build_ext has to be imported after setuptools
try:
    from numpy.distutils.command.build_ext import build_ext  # noqa

    class build_ext_subclass(build_ext):
        def finalize_options(self):
            super().finalize_options()
            if self.parallel is None:
                # Do not override self.parallel if already defined by
                # command-line flag (--parallel or -j)

                parallel = os.environ.get("SKLEARN_BUILD_PARALLEL")
                if parallel:
                    self.parallel = int(parallel)
            if self.parallel:
                print("setting parallel=%d " % self.parallel)

        def build_extensions(self):
            from sklearn._build_utils.openmp_helpers import get_openmp_flag

            for ext in self.extensions:
                if ext.name in USE_NEWEST_NUMPY_C_API:
                    print(
                        f"Use newest NumPy C API for extension {ext.name}"
                    )
                    ext.define_macros.append(DEFINE_MACRO_NUMPY_C_API)
                else:
                    print(f"Use old NumPy C API (version 1.7) for extension {ext.name}")

            if sklearn._OPENMP_SUPPORTED:
                openmp_flag = get_openmp_flag(self.compiler)

                for e in self.extensions:
                    e.extra_compile_args += openmp_flag
                    e.extra_link_args += openmp_flag

            build_ext.build_extensions(self)

    cmdclass["build_ext"] = build_ext_subclass

except ImportError:
    # Numpy should not be a dependency just to be able to introspect
    # that python 3.8 is required.
    pass


def configuration(parent_package="", top_path=None):
    if os.path.exists("MANIFEST"):
        os.remove("MANIFEST")

    from numpy.distutils.misc_util import Configuration
    from sklearn._build_utils import _check_cython_version

    config = Configuration(None, parent_package, top_path)

    # Avoid useless msg:
    # "Ignoring attempt to set 'name' (from ... "
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )

    # Cython is required by config.add_subpackage for templated extensions
    # that need the tempita sub-submodule. So check that we have the correct
    # version of Cython so as to be able to raise a more informative error
    # message from the start if it's not the case.
    _check_cython_version()

    config.add_subpackage("sklearn")

    return config


def check_package_status(package, min_version):
    """
    Returns a dictionary containing a boolean specifying whether given package
    is up-to-date, along with the version string (empty string if
    not installed).
    """
    package_status = {}
    try:
        module = importlib.import_module(package)
        package_version = module.__version__
        package_status["up_to_date"] = parse_version(package_version) >= parse_version(
            min_version
        )
        package_status["version"] = package_version
    except ImportError:
        traceback.print_exc()
        package_status["up_to_date"] = False
        package_status["version"] = ""

    req_str = "scikit-learn requires {} >= {}.\n".format(package, min_version)

    instructions = (
        "Installation instructions are available on the "
        "scikit-learn website: "
        "http://scikit-learn.org/stable/install.html\n"
    )

    if package_status["up_to_date"] is False:
        if package_status["version"]:
            raise ImportError(
                "Your installation of {} {} is out-of-date.\n{}{}".format(
                    package, package_status["version"], req_str, instructions
                )
            )
        else:
            raise ImportError(
                "{} is not installed.\n{}{}".format(package, req_str, instructions)
            )


def setup_package():

    # TODO: Require Python 3.8 for PyPy when PyPy3.8 is ready
    # https://github.com/conda-forge/conda-forge-pinning-feedstock/issues/2089
    if platform.python_implementation() == "PyPy":
        python_requires = ">=3.7"
        required_python_version = (3, 7)
    else:
        python_requires = ">=3.8"
        required_python_version = (3, 8)

    metadata = dict(
        name=DISTNAME,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        license=LICENSE,
        url=URL,
        download_url=DOWNLOAD_URL,
        project_urls=PROJECT_URLS,
        version=VERSION,
        long_description=LONG_DESCRIPTION,
        classifiers=[
            "Intended Audience :: Science/Research",
            "Intended Audience :: Developers",
            "License :: OSI Approved :: BSD License",
            "Programming Language :: C",
            "Programming Language :: Python",
            "Topic :: Software Development",
            "Topic :: Scientific/Engineering",
            "Development Status :: 5 - Production/Stable",
            "Operating System :: Microsoft :: Windows",
            "Operating System :: POSIX",
            "Operating System :: Unix",
            "Operating System :: MacOS",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: Implementation :: PyPy",
        ],
        cmdclass=cmdclass,
        python_requires=python_requires,
        install_requires=min_deps.tag_to_packages["install"],
        package_data={"": ["*.pxd"]},
        **extra_setuptools_args,
    )

    commands = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    if all(
        command in ("egg_info", "dist_info", "clean", "check") for command in commands
    ):
        # These actions are required to succeed without Numpy for example when
        # pip is used to install Scikit-learn when Numpy is not yet present in
        # the system.

        # These commands use setup from setuptools
        from setuptools import setup

        metadata["version"] = VERSION
        metadata["packages"] = ["sklearn"]
    else:
        if sys.version_info < required_python_version:
            required_version = "%d.%d" % required_python_version
            raise RuntimeError(
                "Scikit-learn requires Python %s or later. The current"
                " Python version is %s installed in %s."
                % (required_version, platform.python_version(), sys.executable)
            )

        check_package_status("numpy", min_deps.NUMPY_MIN_VERSION)

        check_package_status("scipy", min_deps.SCIPY_MIN_VERSION)

        # These commands require the setup from numpy.distutils because they
        # may use numpy.distutils compiler classes.
        from numpy.distutils.core import setup

        # Monkeypatches CCompiler.spawn to prevent random wheel build errors on Windows
        # The build errors on Windows was because msvccompiler spawn was not threadsafe
        # This fixed can be removed when we build with numpy >= 1.22.2 on Windows.
        # https://github.com/pypa/distutils/issues/5
        # https://github.com/scikit-learn/scikit-learn/issues/22310
        # https://github.com/numpy/numpy/pull/20640
        from numpy.distutils.ccompiler import replace_method
        from distutils.ccompiler import CCompiler
        from sklearn.externals._numpy_compiler_patch import CCompiler_spawn

        replace_method(CCompiler, "spawn", CCompiler_spawn)

        metadata["configuration"] = configuration

    setup(**metadata)


if __name__ == "__main__":
    setup_package()
