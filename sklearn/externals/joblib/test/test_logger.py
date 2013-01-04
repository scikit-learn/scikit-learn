"""
Test the logger module.
"""

# Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
# Copyright (c) 2009 Gael Varoquaux
# License: BSD Style, 3 clauses.

import shutil
import os
import sys
import io
from tempfile import mkdtemp
import re

from ..logger import PrintTime

try:
    # Python 2/Python 3 compat
    unicode('str')
except NameError:
    unicode = lambda s: s

###############################################################################
# Test fixtures
env = dict()


def setup():
    """ Test setup.
    """
    cachedir = mkdtemp()
    if os.path.exists(cachedir):
        shutil.rmtree(cachedir)
    env['dir'] = cachedir


def teardown():
    """ Test teardown.
    """
    #return True
    shutil.rmtree(env['dir'])


###############################################################################
# Tests
def test_print_time():
    # A simple smoke test for PrintTime.
    try:
        orig_stderr = sys.stderr
        sys.stderr = io.StringIO()
        print_time = PrintTime(logfile=os.path.join(env['dir'], 'test.log'))
        print_time(unicode('Foo'))
        # Create a second time, to smoke test log rotation.
        print_time = PrintTime(logfile=os.path.join(env['dir'], 'test.log'))
        print_time(unicode('Foo'))
        # And a third time
        print_time = PrintTime(logfile=os.path.join(env['dir'], 'test.log'))
        print_time(unicode('Foo'))
        printed_text = sys.stderr.getvalue()
        # Use regexps to be robust to time variations
        match = r"Foo: 0\..s, 0\.0min\nFoo: 0\..s, 0.0min\nFoo: " + \
                r".\..s, 0.0min\n"
        if not re.match(match, printed_text):
            raise AssertionError('Excepted %s, got %s' %
                                    (match, printed_text))
    finally:
        sys.stderr = orig_stderr
