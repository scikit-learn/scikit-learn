"""
Entrypoint for testing from the top-level namespace
"""
import os
import sys

PKG = os.path.dirname(os.path.dirname(__file__))


def test(extra_args=None):
    try:
        import pytest
    except ImportError:
        raise ImportError("Need pytest>=3.0 to run tests")
    cmd = ['--skip-slow', '--skip-network']
    if extra_args:
        if not isinstance(extra_args, list):
            extra_args = [extra_args]
        cmd = extra_args
    cmd += [PKG]
    print("running: pytest {}".format(' '.join(cmd)))
    sys.exit(pytest.main(cmd))


__all__ = ['test']
