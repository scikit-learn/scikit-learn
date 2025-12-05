"""All minimum dependencies for scikit-learn."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import argparse
import re
from collections import defaultdict

# scipy and cython should be in sync with pyproject.toml
NUMPY_MIN_VERSION = "1.24.1"
SCIPY_MIN_VERSION = "1.10.0"
JOBLIB_MIN_VERSION = "1.3.0"
THREADPOOLCTL_MIN_VERSION = "3.2.0"
PYTEST_MIN_VERSION = "7.1.2"
CYTHON_MIN_VERSION = "3.1.2"

# Allowed tags for sanity checking
VALID_TAGS = {
    "build", "install", "benchmark", "docs",
    "examples", "tests", "maintenance"
}

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
    "matplotlib": ("3.6.1", "benchmark, docs, examples, tests"),
    "scikit-image": ("0.22.0", "docs, examples"),
    "pandas": ("1.5.0", "benchmark, docs, examples, tests"),
    "seaborn": ("0.13.0", "docs, examples"),
    "memory_profiler": ("0.57.0", "benchmark, docs"),
    "pytest": (PYTEST_MIN_VERSION, "tests"),
    "pytest-cov": ("2.9.0", "tests"),
    "ruff": ("0.11.7", "tests"),
    "mypy": ("1.15", "tests"),
    "pyamg": ("5.0.0", "tests"),
    "polars": ("0.20.30", "docs, tests"),
    "pyarrow": ("12.0.0", "tests"),
    "sphinx": ("7.3.7", "docs"),
    "sphinx-copybutton": ("0.5.2", "docs"),
    "sphinx-gallery": ("0.17.1", "docs"),
    "numpydoc": ("1.2.0", "docs, tests"),
    "Pillow": ("10.1.0", "docs"),
    "pooch": ("1.8.0", "docs, examples, tests"),
    "sphinx-prompt": ("1.4.0", "docs"),
    "sphinxext-opengraph": ("0.9.1", "docs"),
    "plotly": ("5.18.0", "docs, examples"),
    "sphinxcontrib-sass": ("0.3.4", "docs"),
    "sphinx-remove-toctrees": ("1.0.0.post1", "docs"),
    "sphinx-design": ("0.6.0", "docs"),
    "pydata-sphinx-theme": ("0.15.3", "docs"),
    "towncrier": ("24.8.0", "docs"),
    "conda-lock": ("3.0.1", "maintenance"),
}

# version sanity checking
_VERSION_PATTERN = re.compile(r"^\d+(\.\d+){1,2}(\.\w+)?$")

def _validate_version_format():
    """Ensure all dependency versions use valid semantic versioning."""
    for pkg, (version, _) in dependent_packages.items():
        if not _VERSION_PATTERN.match(version):
            raise ValueError(f"Invalid version format for {pkg}: {version}")

# tag sanity checking
def _validate_tags():
    for pkg, (_, tags) in dependent_packages.items():
        for tag in tags.split(", "):
            if tag not in VALID_TAGS:
                raise ValueError(f"Unknown tag '{tag}' in {pkg}")

_validate_version_format()
_validate_tags()

# Inverse mapping for setuptools
tag_to_packages: dict = defaultdict(list)
for package, (min_version, extras) in dependent_packages.items():
    for extra in extras.split(", "):
        tag_to_packages[extra].append(f"{package}>={min_version}")

# reusable helper
def get_min_version(package: str) -> str:
    return dependent_packages[package][0]


# Used by CI to get the min dependencies
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get min dependencies for a package")
    parser.add_argument("package", choices=dependent_packages)
    args = parser.parse_args()
    print(get_min_version(args.package))
