import os, sys
import pytest
import warnings
import shutil
import subprocess

try:
    import cffi
except ImportError:
    cffi = None

if sys.flags.optimize > 1:
    # no docstrings present to inspect when PYTHONOPTIMIZE/Py_OptimizeFlag > 1
    # cffi cannot succeed
    cffi = None

try:
    with warnings.catch_warnings(record=True) as w:
        # numba issue gh-4733
        warnings.filterwarnings('always', '', DeprecationWarning)
        import numba
except ImportError:
    numba = None

try:
    import cython
    from Cython.Compiler.Version import version as cython_version
except ImportError:
    cython = None
else:
    from distutils.version import LooseVersion
    # Cython 0.29.14 is required for Python 3.8 and there are
    # other fixes in the 0.29 series that are needed even for earlier
    # Python versions.
    # Note: keep in sync with the one in pyproject.toml
    required_version = LooseVersion('0.29.14')
    if LooseVersion(cython_version) < required_version:
        # too old or wrong cython, skip the test
        cython = None

@pytest.mark.skipif(cython is None, reason="requires cython")
@pytest.mark.slow
def test_cython(tmp_path):
    examples = os.path.join(os.path.dirname(__file__), '..', '_examples')
    # CPython 3.5 and below does not handle __fspath__ well: see bpo-26027
    shutil.copytree(examples, str(tmp_path / '_examples'))
    subprocess.check_call([sys.executable, 'setup.py', 'build'],
                          cwd=str(tmp_path / '_examples' / 'cython'))

@pytest.mark.skipif(numba is None or cffi is None,
                    reason="requires numba and cffi")
def test_numba():
    from numpy.random._examples.numba import extending

@pytest.mark.skipif(cffi is None, reason="requires cffi")
def test_cffi():
    from numpy.random._examples.cffi import extending
