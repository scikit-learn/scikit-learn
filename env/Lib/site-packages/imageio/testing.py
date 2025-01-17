# -*- coding: utf-8 -*-
# Distributed under the (new) BSD License. See LICENSE.txt for more info.

""" Functionality used for testing. This code itself is not covered in tests.
"""

import os
import sys
import pytest

# Get root dir
THIS_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = THIS_DIR
for i in range(9):
    ROOT_DIR = os.path.dirname(ROOT_DIR)
    if os.path.isfile(os.path.join(ROOT_DIR, ".gitignore")):
        break


# Functions to use from invoke tasks


def test_unit(cov_report="term"):
    """Run all unit tests. Returns exit code."""
    orig_dir = os.getcwd()
    os.chdir(ROOT_DIR)
    try:
        _clear_imageio()
        _enable_faulthandler()
        return pytest.main(
            [
                "-v",
                "--cov",
                "imageio",
                "--cov-config",
                ".coveragerc",
                "--cov-report",
                cov_report,
                "tests",
            ]
        )
    finally:
        os.chdir(orig_dir)
        import imageio

        print("Tests were performed on", str(imageio))


# Requirements


def _enable_faulthandler():
    """Enable faulthandler (if we can), so that we get tracebacks
    on segfaults.
    """
    try:
        import faulthandler

        faulthandler.enable()
        print("Faulthandler enabled")
    except Exception:
        print("Could not enable faulthandler")


def _clear_imageio():
    # Remove ourselves from sys.modules to force an import
    for key in list(sys.modules.keys()):
        if key.startswith("imageio"):
            del sys.modules[key]
