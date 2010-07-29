"""
Helper for testing.
"""

import sys
import warnings

def warnings_to_stdout():
    """ Redirect all warnings to stdout.
    """
    showwarning_orig = warnings.showwarning
    def showwarning(msg, cat, fname, lno, file=None):
        showwarning_orig(msg, cat, fname.split('/')[-1], 0, sys.stdout)
    warnings.showwarning = showwarning
    #warnings.simplefilter('always')


