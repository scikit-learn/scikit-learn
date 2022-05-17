"""All minimum dependencies for scikit-learn."""
import platform
import argparse


# scipy and cython should by in sync with pyproject.toml

# NumPy version should match oldest-supported-numpy for the minimum supported
# Python version.
# see: https://github.com/scipy/oldest-supported-numpy/blob/main/setup.cfg
if platform.python_implementation() == "PyPy":
    NUMPY_MIN_VERSION = "1.19.2"
else:
    NUMPY_MIN_VERSION = "1.17.3"

SCIPY_MIN_VERSION = "1.3.2"
JOBLIB_MIN_VERSION = "1.0.0"
THREADPOOLCTL_MIN_VERSION = "2.0.0"
PYTEST_MIN_VERSION = "5.0.1"
CYTHON_MIN_VERSION = "0.29.24"


# 'build' and 'install' is included to have structured metadata for CI.
# It will NOT be included in setup's extras_require
# The values are (version_spec, comma separated tags)
dependent_packages = {
    "numpy": (NUMPY_MIN_VERSION, "build, install"),
    "scipy": (SCIPY_MIN_VERSION, "build, install"),
    "joblib": (JOBLIB_MIN_VERSION, "install"),
    "threadpoolctl": (THREADPOOLCTL_MIN_VERSION, "install"),
    "cython": (CYTHON_MIN_VERSION, "build"),
    "matplotlib": ("3.1.2", "benchmark, docs, examples, tests"),
    "scikit-image": ("0.14.5", "docs, examples, tests"),
    "pandas": ("1.0.5", "benchmark, docs, examples, tests"),
    "seaborn": ("0.9.0", "docs, examples"),
    "memory_profiler": ("0.57.0", "benchmark, docs"),
    "pytest": (PYTEST_MIN_VERSION, "tests"),
    "pytest-cov": ("2.9.0", "tests"),
    "flake8": ("3.8.2", "tests"),
    "black": ("22.3.0", "tests"),
    "mypy": ("0.770", "tests"),
    "pyamg": ("4.0.0", "tests"),
    "sphinx": ("4.0.1", "docs"),
    "sphinx-gallery": ("0.7.0", "docs"),
    "numpydoc": ("1.2.0", "docs, tests"),
    "Pillow": ("7.1.2", "docs"),
    "sphinx-prompt": ("1.3.0", "docs"),
    "sphinxext-opengraph": ("0.4.2", "docs"),
}


# create inverse mapping for setuptools
tag_to_packages: dict = {
    extra: []
    for extra in ["build", "install", "docs", "examples", "tests", "benchmark"]
}
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
