"""All minimum dependencies for scikit-learn."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import argparse
from collections import defaultdict

# scipy and cython should by in sync with pyproject.toml
NUMPY_MIN_VERSION = "1.22.0"
SCIPY_MIN_VERSION = "1.8.0"
JOBLIB_MIN_VERSION = "1.2.0"
THREADPOOLCTL_MIN_VERSION = "3.1.0"
PYTEST_MIN_VERSION = "7.1.2"
CYTHON_MIN_VERSION = "3.0.10"


# 'build' and 'install' is included to have structured metadata for CI.
# It will NOT be included in setup's extras_require
# The values are (version_spec, comma separated tags)
dependent_packages = {
    "numpy": (NUMPY_MIN_VERSION, "build, install"),
    "scipy": (SCIPY_MIN_VERSION, "build, install"),
    "joblib": (JOBLIB_MIN_VERSION, "install"),
    "threadpoolctl": (THREADPOOLCTL_MIN_VERSION, "install"),
    "cython": (CYTHON_MIN_VERSION, "build"),
    "meson-python": ("0.17.1", "build"),
    "matplotlib": ("3.5.0", "benchmark, docs, examples, tests"),
    "scikit-image": ("0.19.0", "docs, examples, tests"),
    "pandas": ("1.4.0", "benchmark, docs, examples, tests"),
    "seaborn": ("0.9.0", "docs, examples"),
    "memory_profiler": ("0.57.0", "benchmark, docs"),
    "pytest": (PYTEST_MIN_VERSION, "tests"),
    "pytest-cov": ("2.9.0", "tests"),
    "ruff": ("0.11.7", "tests"),
    "mypy": ("1.15", "tests"),
    "pyamg": ("4.2.1", "tests"),
    "polars": ("0.20.30", "docs, tests"),
    "pyarrow": ("12.0.0", "tests"),
    "sphinx": ("7.3.7", "docs"),
    "sphinx-copybutton": ("0.5.2", "docs"),
    "sphinx-gallery": ("0.17.1", "docs"),
    "numpydoc": ("1.2.0", "docs, tests"),
    "Pillow": ("8.4.0", "docs"),
    "pooch": ("1.6.0", "docs, examples, tests"),
    "sphinx-prompt": ("1.4.0", "docs"),
    "sphinxext-opengraph": ("0.9.1", "docs"),
    "plotly": ("5.14.0", "docs, examples"),
    "sphinxcontrib-sass": ("0.3.4", "docs"),
    "sphinx-remove-toctrees": ("1.0.0.post1", "docs"),
    "sphinx-design": ("0.6.0", "docs"),
    "pydata-sphinx-theme": ("0.15.3", "docs"),
    "towncrier": ("24.8.0", "docs"),
    # XXX: Pin conda-lock to the latest released version (needs manual update
    # from time to time)
    "conda-lock": ("3.0.1", "maintenance"),
}


# create inverse mapping for setuptools
tag_to_packages: dict = defaultdict(list)
for package, (min_version, extras) in dependent_packages.items():
    for extra in extras.split(", "):
        tag_to_packages[extra].append("{}>={}".format(package, min_version))


# Used by CI to get the min dependencies
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get min dependencies for a package")

    parser.add_argument("package", choices=dependent_packages)
    args = parser.parse_args()
    min_version = dependent_packages[args.package][0]
    print(min_version)
