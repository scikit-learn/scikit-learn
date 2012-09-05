"""
Fixtures for the tests.

This module gets loaded by test discovery scanners (such as nose) in
their collection scan.
"""

import platform
import os
if platform.system() == 'Windows':
    # During the tests, under Windows, we don't want any multiprocessing
    os.environ['JOBLIB_MULTIPROCESSING'] = '0'
