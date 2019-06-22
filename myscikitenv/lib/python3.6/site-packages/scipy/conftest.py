# Pytest customization
from __future__ import division, absolute_import, print_function

import os
import pytest
import warnings

from distutils.version import LooseVersion
from scipy._lib._fpumode import get_fpu_mode
from scipy._lib._testutils import FPUModeChangeWarning


def pytest_runtest_setup(item):
    if LooseVersion(pytest.__version__) >= LooseVersion("3.6.0"):
        mark = item.get_closest_marker("xslow")
    else:
        mark = item.get_marker("xslow")
    if mark is not None:
        try:
            v = int(os.environ.get('SCIPY_XSLOW', '0'))
        except ValueError:
            v = False
        if not v:
            pytest.skip("very slow test; set environment variable SCIPY_XSLOW=1 to run it")


@pytest.fixture(scope="function", autouse=True)
def check_fpu_mode(request):
    """
    Check FPU mode was not changed during the test.
    """
    old_mode = get_fpu_mode()
    yield
    new_mode = get_fpu_mode()

    if old_mode != new_mode:
        warnings.warn("FPU mode changed from {0:#x} to {1:#x} during "
                      "the test".format(old_mode, new_mode),
                      category=FPUModeChangeWarning, stacklevel=0)
