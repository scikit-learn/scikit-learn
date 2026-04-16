import multiprocessing
import subprocess
import sys

import pytest

from .test_public_api import PUBLIC_MODULES

# Regression tests for gh-6793.
# Check that all modules are importable in a new Python process.
# This is not necessarily true if there are import cycles present.


def _check_single_module(module):
    pid = subprocess.Popen([sys.executable, '-X', 'faulthandler', '-c',
                            f'import {module}'])
    assert pid.wait() == 0, f'Failed to import {module}'


@pytest.mark.fail_slow(40)
@pytest.mark.slow
def test_public_modules_importable_2():
    # Ensure we use max 6 processes, to limit peak resource usage (memory, file handles)
    # on resource-constrained systems (e.g., RISC-V - see gh-24163).
    with multiprocessing.Pool(processes=6) as pool:
        pool.map(_check_single_module, PUBLIC_MODULES)
