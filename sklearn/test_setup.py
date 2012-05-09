"""
Fixtures for the tests.

This module gets loaded by test discovery scanners (such as nose) in
their collection scan.
"""

import os
if os.name.startswith('win') or os.name.startswith('nt'):
    # During the tests, under Windows, we don't want any multiprocessing
    os.environ['JOBLIB_MULTIPROCESSING'] = '0'
