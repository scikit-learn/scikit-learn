#! /usr/bin/env python
#
# Copyright (C) 2007-2009 Cournapeau David <cournape@gmail.com>
#               2010 Fabian Pedregosa <fabian.pedregosa@inria.fr>
# License: 3-clause BSD

import sys
import os
from os.path import join
import platform
import shutil

from setuptools import Command, Extension, setup
from setuptools.command.build_ext import build_ext

import traceback
import importlib

try:
    import builtins
except ImportError:
    # Python 2 compat: just to be able to declare that Python >=3.8 is needed.
    import __builtin__ as builtins

# This is a bit (!) hackish: we are setting a global variable so that the main
# sklearn __init__ can detect if it is being loaded by the setup routine, to
# avoid attempting to load components that aren't built yet.
# TODO: can this be simplified or removed since the switch to setuptools
# away from numpy.distutils?
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
from sklearn._build_utils import _check_cython_version  # noqa
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
    "sklearn.cluster._dbscan_inner",
    "sklearn.cluster._hierarchical_fast",
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
    "sklearn.linear_model._sag_fast",
    "sklearn.linear_model._sgd_fast",
    "sklearn.manifold._barnes_hut_tsne",
    "sklearn.manifold._utils",
    "sklearn.metrics.cluster._expected_mutual_info_fast",
    "sklearn.metrics._pairwise_distances_reduction._datasets_pair",
    "sklearn.metrics._pairwise_distances_reduction._middle_term_computer",
    "sklearn.metrics._pairwise_distances_reduction._base",
    "sklearn.metrics._pairwise_distances_reduction._argkmin",
    "sklearn.metrics._pairwise_distances_reduction._radius_neighbors",
    "sklearn.metrics._pairwise_fast",
    "sklearn.neighbors._ball_tree",
    "sklearn.neighbors._kd_tree",
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
    "sklearn.utils.murmurhash",
    "sklearn.svm._newrand",
    "sklearn._isotonic",
)


# Custom clean command to remove build artifacts


class CleanCommand(Command):
    description = "Remove build artifacts from the source tree"

    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
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


# Custom build_ext command to set OpenMP compile flags depending on os and
# compiler. Also makes it possible to set the parallelism level via
# and environment variable (useful for the wheel building CI).
# build_ext has to be imported after setuptools


class build_ext_subclass(build_ext):
    def finalize_options(self):
        build_ext.finalize_options(self)
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
                print(f"Using newest NumPy C API for extension {ext.name}")
                ext.define_macros.append(DEFINE_MACRO_NUMPY_C_API)
            else:
                print(f"Using old NumPy C API (version 1.7) for extension {ext.name}")

        if sklearn._OPENMP_SUPPORTED:
            openmp_flag = get_openmp_flag(self.compiler)

            for e in self.extensions:
                e.extra_compile_args += openmp_flag
                e.extra_link_args += openmp_flag

        build_ext.build_extensions(self)

    def run(self):
        # Specifying `build_clib` allows running `python setup.py develop`
        # fully from a fresh clone.
        self.run_command("build_clib")
        build_ext.run(self)


cmdclass = {
    "clean": CleanCommand,
    "build_ext": build_ext_subclass,
}


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


extension_config = {
    "__check_build": [
        {"sources": ["_check_build.pyx"]},
    ],
    "": [
        {"sources": ["_isotonic.pyx"], "include_np": True},
    ],
    "_loss": [
        {"sources": ["_loss.pyx.tp"], "include_np": True},
    ],
    "cluster": [
        {"sources": ["_dbscan_inner.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_hierarchical_fast.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_k_means_common.pyx"], "include_np": True},
        {"sources": ["_k_means_lloyd.pyx"], "include_np": True},
        {"sources": ["_k_means_elkan.pyx"], "include_np": True},
        {"sources": ["_k_means_minibatch.pyx"], "include_np": True},
    ],
    "datasets": [
        {
            "sources": ["_svmlight_format_fast.pyx"],
            "include_np": True,
            "compile_for_pypy": False,
        }
    ],
    "decomposition": [
        {"sources": ["_online_lda_fast.pyx"], "include_np": True},
        {"sources": ["_cdnmf_fast.pyx"], "include_np": True},
    ],
    "ensemble": [
        {"sources": ["_gradient_boosting.pyx"], "include_np": True},
    ],
    "ensemble._hist_gradient_boosting": [
        {"sources": ["_gradient_boosting.pyx"], "include_np": True},
        {"sources": ["histogram.pyx"], "include_np": True},
        {"sources": ["splitting.pyx"], "include_np": True},
        {"sources": ["_binning.pyx"], "include_np": True},
        {"sources": ["_predictor.pyx"], "include_np": True},
        {"sources": ["_bitset.pyx"], "include_np": True},
        {"sources": ["common.pyx"], "include_np": True},
        {"sources": ["utils.pyx"], "include_np": True},
    ],
    "feature_extraction": [
        {"sources": ["_hashing_fast.pyx"], "language": "c++", "include_np": True},
    ],
    "linear_model": [
        {"sources": ["_cd_fast.pyx"], "include_np": True},
        {"sources": ["_sgd_fast.pyx"], "include_np": True},
        {"sources": ["_sag_fast.pyx.tp"], "include_np": True},
    ],
    "manifold": [
        {"sources": ["_utils.pyx"], "include_np": True},
        {"sources": ["_barnes_hut_tsne.pyx"], "include_np": True},
    ],
    "metrics": [
        {"sources": ["_pairwise_fast.pyx"], "include_np": True},
        {
            "sources": ["_dist_metrics.pyx.tp", "_dist_metrics.pxd.tp"],
            "include_np": True,
        },
    ],
    "metrics.cluster": [
        {"sources": ["_expected_mutual_info_fast.pyx"], "include_np": True},
    ],
    "metrics._pairwise_distances_reduction": [
        {
            "sources": ["_datasets_pair.pyx.tp", "_datasets_pair.pxd.tp"],
            "language": "c++",
            "include_np": True,
            "extra_compile_args": ["-std=c++11"],
        },
        {
            "sources": ["_middle_term_computer.pyx.tp", "_middle_term_computer.pxd.tp"],
            "language": "c++",
            "include_np": True,
            "extra_compile_args": ["-std=c++11"],
        },
        {
            "sources": ["_base.pyx.tp", "_base.pxd.tp"],
            "language": "c++",
            "include_np": True,
            "extra_compile_args": ["-std=c++11"],
        },
        {
            "sources": ["_argkmin.pyx.tp", "_argkmin.pxd.tp"],
            "language": "c++",
            "include_np": True,
            "extra_compile_args": ["-std=c++11"],
        },
        {
            "sources": ["_radius_neighbors.pyx.tp", "_radius_neighbors.pxd.tp"],
            "language": "c++",
            "include_np": True,
            "extra_compile_args": ["-std=c++11"],
        },
    ],
    "preprocessing": [
        {"sources": ["_csr_polynomial_expansion.pyx"], "include_np": True},
    ],
    "neighbors": [
        {"sources": ["_ball_tree.pyx"], "include_np": True},
        {"sources": ["_kd_tree.pyx"], "include_np": True},
        {"sources": ["_partition_nodes.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_quad_tree.pyx"], "include_np": True},
    ],
    "svm": [
        {
            "sources": ["_newrand.pyx"],
            "include_np": True,
            "include_dirs": [join("src", "newrand")],
            "language": "c++",
            # Use C++11 random number generator fix
            "extra_compile_args": ["-std=c++11"],
        },
        {
            "sources": ["_libsvm.pyx"],
            "depends": [
                join("src", "libsvm", "libsvm_helper.c"),
                join("src", "libsvm", "libsvm_template.cpp"),
                join("src", "libsvm", "svm.cpp"),
                join("src", "libsvm", "svm.h"),
                join("src", "newrand", "newrand.h"),
            ],
            "include_dirs": [
                join("src", "libsvm"),
                join("src", "newrand"),
            ],
            "libraries": ["libsvm-skl"],
            "extra_link_args": ["-lstdc++"],
            "include_np": True,
        },
        {
            "sources": ["_liblinear.pyx"],
            "libraries": ["liblinear-skl"],
            "include_dirs": [
                join("src", "liblinear"),
                join("src", "newrand"),
                join("..", "utils"),
            ],
            "include_np": True,
            "depends": [
                join("src", "liblinear", "tron.h"),
                join("src", "liblinear", "linear.h"),
                join("src", "liblinear", "liblinear_helper.c"),
                join("src", "newrand", "newrand.h"),
            ],
            "extra_link_args": ["-lstdc++"],
        },
        {
            "sources": ["_libsvm_sparse.pyx"],
            "libraries": ["libsvm-skl"],
            "include_dirs": [
                join("src", "libsvm"),
                join("src", "newrand"),
            ],
            "include_np": True,
            "depends": [
                join("src", "libsvm", "svm.h"),
                join("src", "newrand", "newrand.h"),
                join("src", "libsvm", "libsvm_sparse_helper.c"),
            ],
            "extra_link_args": ["-lstdc++"],
        },
    ],
    "tree": [
        {"sources": ["_tree.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_splitter.pyx"], "include_np": True},
        {"sources": ["_criterion.pyx"], "include_np": True},
        {"sources": ["_utils.pyx"], "include_np": True},
    ],
    "utils": [
        {"sources": ["sparsefuncs_fast.pyx"], "include_np": True},
        {"sources": ["_cython_blas.pyx"]},
        {"sources": ["arrayfuncs.pyx"], "include_np": True},
        {
            "sources": ["murmurhash.pyx", join("src", "MurmurHash3.cpp")],
            "include_dirs": ["src"],
            "include_np": True,
        },
        {"sources": ["_fast_dict.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_fast_dict.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_openmp_helpers.pyx"]},
        {"sources": ["_seq_dataset.pyx.tp", "_seq_dataset.pxd.tp"], "include_np": True},
        {
            "sources": ["_weight_vector.pyx.tp", "_weight_vector.pxd.tp"],
            "include_np": True,
        },
        {"sources": ["_random.pyx"], "include_np": True},
        {"sources": ["_logistic_sigmoid.pyx"], "include_np": True},
        {"sources": ["_readonly_array_wrapper.pyx"], "include_np": True},
        {"sources": ["_typedefs.pyx"], "include_np": True},
        {"sources": ["_heap.pyx"], "include_np": True},
        {"sources": ["_sorting.pyx"], "include_np": True},
        {"sources": ["_vector_sentinel.pyx"], "language": "c++", "include_np": True},
        {"sources": ["_isfinite.pyx"]},
    ],
}

# Paths in `libraries` must be relative to the root directory because `libraries` is
# passed directly to `setup`
libraries = [
    (
        "libsvm-skl",
        {
            "sources": [
                join("sklearn", "svm", "src", "libsvm", "libsvm_template.cpp"),
            ],
            "depends": [
                join("sklearn", "svm", "src", "libsvm", "svm.cpp"),
                join("sklearn", "svm", "src", "libsvm", "svm.h"),
                join("sklearn", "svm", "src", "newrand", "newrand.h"),
            ],
            # Use C++11 to use the random number generator fix
            "extra_compiler_args": ["-std=c++11"],
            "extra_link_args": ["-lstdc++"],
        },
    ),
    (
        "liblinear-skl",
        {
            "sources": [
                join("sklearn", "svm", "src", "liblinear", "linear.cpp"),
                join("sklearn", "svm", "src", "liblinear", "tron.cpp"),
            ],
            "depends": [
                join("sklearn", "svm", "src", "liblinear", "linear.h"),
                join("sklearn", "svm", "src", "liblinear", "tron.h"),
                join("sklearn", "svm", "src", "newrand", "newrand.h"),
            ],
            # Use C++11 to use the random number generator fix
            "extra_compiler_args": ["-std=c++11"],
            "extra_link_args": ["-lstdc++"],
        },
    ),
]


def configure_extension_modules():
    # Skip cythonization as we do not want to include the generated
    # C/C++ files in the release tarballs as they are not necessarily
    # forward compatible with future versions of Python for instance.
    if "sdist" in sys.argv or "--help" in sys.argv:
        return []

    from sklearn._build_utils import cythonize_extensions
    from sklearn._build_utils import gen_from_templates
    import numpy

    is_pypy = platform.python_implementation() == "PyPy"
    np_include = numpy.get_include()

    optimization_level = "O2"
    if os.name == "posix":
        default_extra_compile_args = [f"-{optimization_level}"]
        default_libraries = ["m"]
    else:
        default_extra_compile_args = [f"/{optimization_level}"]
        default_libraries = []

    build_with_debug_symbols = (
        os.environ.get("SKLEARN_BUILD_ENABLE_DEBUG_SYMBOLS", "0") != "0"
    )
    if os.name == "posix":
        if build_with_debug_symbols:
            default_extra_compile_args.append("-g")
        else:
            # Setting -g0 will strip symbols, reducing the binary size of extensions
            default_extra_compile_args.append("-g0")

    cython_exts = []
    for submodule, extensions in extension_config.items():
        submodule_parts = submodule.split(".")
        parent_dir = join("sklearn", *submodule_parts)
        for extension in extensions:
            if is_pypy and not extension.get("compile_for_pypy", True):
                continue

            # Generate files with Tempita
            tempita_sources = []
            sources = []
            for source in extension["sources"]:
                source = join(parent_dir, source)
                new_source_path, path_ext = os.path.splitext(source)

                if path_ext != ".tp":
                    sources.append(source)
                    continue

                # `source` is a Tempita file
                tempita_sources.append(source)

                # Do not include pxd files that were generated by tempita
                if os.path.splitext(new_source_path)[-1] == ".pxd":
                    continue
                sources.append(new_source_path)

            gen_from_templates(tempita_sources)

            # By convention, our extensions always use the name of the first source
            source_name = os.path.splitext(os.path.basename(sources[0]))[0]
            if submodule:
                name_parts = ["sklearn", submodule, source_name]
            else:
                name_parts = ["sklearn", source_name]
            name = ".".join(name_parts)

            # Make paths start from the root directory
            include_dirs = [
                join(parent_dir, include_dir)
                for include_dir in extension.get("include_dirs", [])
            ]
            if extension.get("include_np", False):
                include_dirs.append(np_include)

            depends = [
                join(parent_dir, depend) for depend in extension.get("depends", [])
            ]

            extra_compile_args = (
                extension.get("extra_compile_args", []) + default_extra_compile_args
            )
            libraries_ext = extension.get("libraries", []) + default_libraries

            new_ext = Extension(
                name=name,
                sources=sources,
                language=extension.get("language", None),
                include_dirs=include_dirs,
                libraries=libraries_ext,
                depends=depends,
                extra_link_args=extension.get("extra_link_args", None),
                extra_compile_args=extra_compile_args,
            )
            cython_exts.append(new_ext)

    return cythonize_extensions(cython_exts)


def setup_package():
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
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: Implementation :: PyPy",
        ],
        cmdclass=cmdclass,
        python_requires=python_requires,
        install_requires=min_deps.tag_to_packages["install"],
        package_data={"": ["*.csv", "*.gz", "*.txt", "*.pxd", "*.rst", "*.jpg"]},
        zip_safe=False,  # the package can run out of an .egg file
        extras_require={
            key: min_deps.tag_to_packages[key]
            for key in ["examples", "docs", "tests", "benchmark"]
        },
    )

    commands = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    if not all(
        command in ("egg_info", "dist_info", "clean", "check") for command in commands
    ):
        if sys.version_info < required_python_version:
            required_version = "%d.%d" % required_python_version
            raise RuntimeError(
                "Scikit-learn requires Python %s or later. The current"
                " Python version is %s installed in %s."
                % (required_version, platform.python_version(), sys.executable)
            )

        check_package_status("numpy", min_deps.NUMPY_MIN_VERSION)
        check_package_status("scipy", min_deps.SCIPY_MIN_VERSION)

        _check_cython_version()
        metadata["ext_modules"] = configure_extension_modules()
        metadata["libraries"] = libraries
    setup(**metadata)


if __name__ == "__main__":
    setup_package()
