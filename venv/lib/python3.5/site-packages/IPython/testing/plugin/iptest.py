#!/usr/bin/env python
"""Nose-based test runner.
"""

from nose.core import main
from nose.plugins.builtin import plugins
from nose.plugins.doctests import Doctest

from . import ipdoctest
from .ipdoctest import IPDocTestRunner

if __name__ == '__main__':
    print('WARNING: this code is incomplete!')
    print()

    pp = [x() for x in plugins]  # activate all builtin plugins first
    main(testRunner=IPDocTestRunner(),
         plugins=pp+[ipdoctest.IPythonDoctest(),Doctest()])
